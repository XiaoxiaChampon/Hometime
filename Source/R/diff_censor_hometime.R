############################################################
# Copyright 2023 Xiaoxia Champon

# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the “Software”), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:

# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.
############################################################
# Purpose: Selecting the Optimal Approach to Estimate Marginal Mean Effect for Hometime
# Author:  Xiaoxia Champon, special thanks to Laine Thomas, Sean O'Brien
# Date: 04/27/202
##############################################################


##############################################################

# load required libraries
library(MASS)
library(survival)
library(sandwich)
library(rms)
library(Hmisc)
library(splines)
library(quantreg)
library(data.table)
options(warn = -1)

# ---- For: parallelization ----
# For: foreach loop
library(foreach)

run_parallel <- TRUE
time_elapsed <- list()
if(run_parallel)
{
    print("RUNNING PARALLEL")
    
    # For: makeCluster
    library(doParallel)
    
    # For: %dorng% or registerDoRNG for reproducable parallel random number generation
    library(doRNG)
    
    if(exists("initialized_parallel") && initialized_parallel == TRUE)
    {
        parallel::stopCluster(cl = my.cluster)
    }
    n.cores <- parallel::detectCores() - 1
    my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
    doParallel::registerDoParallel(cl = my.cluster)
    cat("Parellel Registered: ", foreach::getDoParRegistered(), "\n")
    initialized_parallel <- TRUE
    

}



###########################################################
#Function to generate hoemtime data based on different scenarios
#Input:
#B: scalar, the number of repetition
#n: scalar, the number of total sample size
# censor=0, no right censor, censor=1, yes ;
# equalsize=0, unequal allocation between treatment and control, equalsize=1, equal allocation
#effect: scalar, treatment effect, usually 0 means no effect, 1 means has some effect
#effect.d: scalar, effect of treatment on death, usually 0 means no effect of death, 1 means has some effect
#diff_censor: scalar, 1 means different censoring rate for treatment and control group, 0, otherwise


#Output: 3D array, n*5*B, n is the number of observation, 
#5 is the column of the data (outcome, group, outcome.t , outcome.b, event), 
#B is the number of replications
generate_home_time_scenario = function(B, n, censor, effect, effect.d, equalsize,diff_censor) {
 
    # 
    home_time_data=array(0,c(n,4,B))
    #same censoring rate
    
    if (censor==1 && diff_censor==0){censor_p=c(0)}
    
    
    # 
    #different censoring rate for treatment and non-treatment
    if (censor==1 && diff_censor==1){censor_p=array(0,c(2,2,B))}
    
    for (j in 1:B) {
        #j=1
        if (censor == 0) {
            admixture <-
                c(rep(3, n / 4), rep(0, n / 4), rep(0, n / 4), rep(-1, n / 4))
            
            frailty <- rnorm(n, admixture, .5)
            
            
            if (equalsize == 1) {
                trt <- rbinom(n, 1, .5)
            }
            if (equalsize == 0) {
                trt <- rbinom(n, 1, .6)
            }
            
            cum.ao <- rep(0, n)
            prev.dead <- rep(0, n)
            
            for (t in 1:365) {
                if (t > 1) {
                    prev.dead <- dead
                }
                
                lin <- -7 + .00015 * t + frailty - effect.d * (t < 14) * trt
                #' *The probability of the death from the linear predictor*
                p.d <- exp(lin) / (1 + exp(lin))
                dead <- ((rbinom(n, 1, p.d) + prev.dead) > 0)
                #' *The probability of out of the hospital from the linear predictor*
                lin <- 2 - 1.5 * frailty + effect * trt
                p.o <- exp(lin) / (1 + exp(lin))
                out <- rbinom(n, 1, p.o)
                
                alive.out <- (1 - dead) * out
                
                cum.ao <- cum.ao + alive.out
                if (t == 365)
                    slice.365 <- cum.ao
            }
            data <- data.frame(outcome = slice.365, group = trt)
            data$outcome.t <- slice.365
            
            # htevent 1: non censored, 0: censored
            data$htevent <- rep(1, n)
        }
        
        
        if (censor == 1) {
            admixture <-
                c(rep(3, n / 4), rep(0, n / 4), rep(0, n / 4), rep(-1, n / 4))
            frailty <- rnorm(n, admixture, .5)
            
            if (equalsize == 1) {
                trt <- rbinom(n, 1, .5)
            }
            if (equalsize == 0) {
                trt <- rbinom(n, 1, .6)
            }
            cum.ao <- rep(0, n)
            prev.dead <- rep(0, n)
            
            prob_dead_plot_data=matrix(0,nrow=400,ncol=n)
            prob_out_plot_data=matrix(0,nrow=400,ncol=n)
            for (t in 1:400) {
                if (t > 1) {
                    prev.dead <- dead
                }
                
                lin <- -7 + .00015 * t + frailty - effect.d * (t < 14) * trt
                p.d <- exp(lin) / (1 + exp(lin))
                
                dead <- ((rbinom(n, 1, p.d) + prev.dead) > 0)
                #######
                prob_dead_plot_data[t,]=p.d
                ########
                lin <- 2 - 1.5 * frailty + effect * trt
                p.o <- exp(lin) / (1 + exp(lin))
                
                #####
                prob_out_plot_data[t,]=p.o
                #####
                out <- rbinom(n, 1, p.o)
                
                alive.out <- (1 - dead) * out
                
                cum.ao <- cum.ao + alive.out
                if (t == 365)
                    slice.365 <- cum.ao
                #' *still need a event indicator for Cox model, even if not censored*
                
            }
            # create censoring indicator
            data <- data.frame(outcome = slice.365, group = trt)
            data$outcome.t <- slice.365
            
            #####################Aug 30, uniform censoring
            if (diff_censor==0){
                uniform_cencoring=runif(n, 180, 365)
                for (subject_index in 1:n){
                    ##add the proportion for the censoring, 70% for treatment and 35% for non-treatment
                    
                    if (data$outcome.t[subject_index]<uniform_cencoring[subject_index]){
                        data$htevent[subject_index] <- 1
                    }
                    
                    if (data$outcome.t[subject_index]>=uniform_cencoring[subject_index]){
                        data$htevent[subject_index] <- 0
                        data$outcome.t[subject_index]=uniform_cencoring[subject_index]
                    }
                    censor_p[j]=1-sum(data$htevent)/n
                }
            }
            
            ###########################################################################
            #Nov 15, 2023 change the percentage of censoring based on treatment and non-treatment
            #########################################################################
            if (diff_censor==1) {
                data_nontrt=data[data$group==0,]
                data_trt=data[data$group==1,]
                n_trt=dim(data_trt)[1]
                n_nontrt=dim(data_nontrt)[1]
                ######treatment group
                ##censoring
                #uniform_cencoring_trt=runif(n_trt, 180, 365)
                for (subject_index in 1:n_trt){
                    uniform_cencoring_trt=runif(n_trt, 240, 365)
                    ##add the proportion for the censoring, 70% for treatment and 50% for non-treatment
                    
                    if (data_trt$outcome.t[subject_index]<uniform_cencoring_trt[subject_index]){
                        data_trt$htevent[subject_index] <- 1
                    }
                    
                    if (data_trt$outcome.t[subject_index]>=uniform_cencoring_trt[subject_index]){
                        data_trt$htevent[subject_index] <- 0
                        data_trt$outcome.t[subject_index]=uniform_cencoring_trt[subject_index]
                    }
                    censor_prop_trt=1-sum(data_trt$htevent)/n
                    
                }
                
                
             
                
                
                ###non-treatment group
                #uniform_cencoring_nontrt=runif(n_nontrt, 180, 365)
                for (subject_index in 1:n_nontrt){
                    uniform_cencoring_nontrt=runif(n_nontrt, 360,365 )
                    ##add the proportion for the censoring, 70% for treatment and 50% for non-treatment
                    
                    if (data_nontrt$outcome.t[subject_index]<uniform_cencoring_nontrt[subject_index]){
                        data_nontrt$htevent[subject_index] <- 1
                    }
                    
                    if (data_nontrt$outcome.t[subject_index]>=uniform_cencoring_nontrt[subject_index]){
                        data_nontrt$htevent[subject_index] <- 0
                        data_nontrt$outcome.t[subject_index]=uniform_cencoring_nontrt[subject_index]
                    }
                    censor_prop_nontrt=1-sum(data_nontrt$htevent)/n
                }
                data=rbind(data_trt,data_nontrt)
                censor_p[1,,j]=c(1,censor_prop_trt)
                censor_p[2,,j]=c(0,censor_prop_nontrt)
            }
            

        }
      
     
       
        
        home_time_data[,,j]=as.matrix(data)
    }
    return(list("home_time_data"=home_time_data,"censor_p"=censor_p))
}




#Function to find the results for simulated home time data
#Input : 3D array, n*5*B hometime data, n-observations, B: replications
#Output: average values (parameter, standard error for parameter and p-value) for B replications

home_time_simulation_results=function(home_time_scenario_data){
    num_replica=dim(home_time_scenario_data)[3]
    power <- matrix(nrow = num_replica, ncol = 6)
    param.est <- matrix(nrow = num_replica, ncol = 6)
    n=dim(home_time_scenario_data)[1]
    
    foreach_out <- foreach (j = 1:num_replica, .combine = cbind, .init = NULL) %dorng% {
        source("functions_hometime.R")
        home_time_data_rep=data.frame(home_time_scenario_data[,,j])
        colnames(home_time_data_rep)=c("outcome", "group","outcome.t","htevent" )
        param.est_j=home_time_regressions(home_time_data_rep)$param.est
        power_j=home_time_regressions(home_time_data_rep)$power
        return(c(param.est_j, power_j))
    }
    
    param.est=t(foreach_out[1:6,])
    power=t(foreach_out[7:12,])
    a.param <- apply(param.est, MARGIN = 2, FUN = mean)
    a.paramsd <- apply(param.est, MARGIN = 2, FUN = sd) / sqrt(n)
    a.power <- apply(power, MARGIN = 2, FUN = mean)
   
    result <- cbind(a.param, a.paramsd, a.power)
  
    rownames(result) <- c("lin", "cox", "med", "nb", "poi", "temp")
 
    colnames(result) <-c("est", "se", "power")
    return(result)
}



#Function to produce the table for manuscript
#Input : B:scalar- the number of replications, n: scalar-the number of observations
#Output: average values (parameter, standard error for parameter and p-value) for B replications
# censor=0, no right censor, censor=1, yes ;
# equalsize=0, unequal allocation between treatment and control, equalsize=1, equal allocation
#effect: scalar, treatment effect, usually 0 means no effect, 1 means has some effect
#effect.d: scalar, effect of treatment on death, usually 0 means no effect of death, 1 means has some effect

#Output: 6 rows (one per model), 3 columns (parameter, standard error, p-value)
home_time_table=function(B, n, censor, effect, effect.d, equalsize,diff_censor){
    home_time_data=generate_home_time_scenario(B, n, censor, effect, effect.d, equalsize,diff_censor)
    result_table=home_time_simulation_results(home_time_data$home_time_data)
    return(result_table)
}



##################################################################################################
####################################################################################################
#censored balanced type I
set.seed(123)
uncensor_balance_500 <- home_time_table(B = 5000, n = 500, censor = 1, effect = 0, effect.d = 0, equalsize = 1,diff_censor=1)


set.seed(123)
uncensor_balance_1000 <- home_time_table(B = 5000, n = 1000, censor = 1, effect = 0, effect.d = 0, equalsize = 1,diff_censor=1)



set.seed(123)
uncensor_balance_5000 <- home_time_table(B = 5000, n = 5000, censor = 1, effect = 0, effect.d = 0, equalsize = 1,diff_censor=1)

##################################################################################
#####power censor balanced
set.seed(123)
uncensor_balance_500 <- home_time_table(B = 5000, n = 500, censor = 1, effect = 1, effect.d = 0, equalsize = 1,diff_censor=1)


set.seed(123)
uncensor_balance_1000 <- home_time_table(B = 5000, n = 1000, censor = 1, effect = 1, effect.d = 0, equalsize = 1,diff_censor=1)



set.seed(123)
uncensor_balance_5000 <- home_time_table(B = 5000, n = 5000, censor = 1, effect = 1, effect.d = 0, equalsize = 1,diff_censor=1)


power_censor_balance=cbind(uncensor_balance_500,uncensor_balance_1000,uncensor_balance_5000)
save(power_censor_balance,file="power_censor_balance.RData")

###############################################################################
#####type I censor unbalanced
set.seed(123)
uncensor_balance_500 <- home_time_table(B = 5000, n = 500, censor = 1, effect = 0, effect.d = 0, equalsize = 0,diff_censor=1)



set.seed(123)
uncensor_balance_1000 <- home_time_table(B = 5000, n = 1000, censor = 1, effect = 0, effect.d = 0, equalsize = 0,diff_censor=1)




set.seed(123)
uncensor_balance_5000 <- home_time_table(B = 5000, n = 5000, censor = 1, effect = 0, effect.d = 0, equalsize = 0,diff_censor=1)






load("/Users/xzhao17/Desktop/hometime-refs/diff_censor_data/typeI_censor_unbalance.RData")
xtable(typeI_censor_unbalance)
##################################################################################
#####power censor unbalanced
set.seed(123)
uncensor_balance_500 <- home_time_table(B = 5000, n = 500, censor = 1, effect = 1, effect.d = 0, equalsize = 0,diff_censor=1)



set.seed(123)
uncensor_balance_1000 <- home_time_table(B = 5000, n = 1000, censor = 1, effect = 1, effect.d = 0, equalsize = 0,diff_censor=1)




set.seed(123)
uncensor_balance_5000 <- home_time_table(B = 5000, n = 5000, censor = 1, effect = 1, effect.d = 0, equalsize = 0,diff_censor=1)



power_censor_unbalance=cbind(uncensor_balance_500,uncensor_balance_1000,uncensor_balance_5000)
save(power_censor_unbalance,file="power_censor_unbalance.RData")



if(run_parallel)
{
    parallel::stopCluster(cl = my.cluster)
    initialized_parallel <- FALSE
}
