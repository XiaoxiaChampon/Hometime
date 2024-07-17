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
#setwd("~/Documents/GitHub/Hometime/Source/R")
# load required libraries
library(MASS)
library(survival)
library(sandwich)
library(rms)
library(Hmisc)
library(splines)
library(quantreg)
library(data.table)
#options(warn = -1)

###########
 library(optparse)
# #B, n, censor, effect, equalsize,diff_censor,censorbig=1
# # Define options

option_list <- list(
        make_option(c("-j", "--jobid"), type="integer", default=123,
                    help="Job Index", metavar="JOBID"),
        make_option(c("-n", "--numcpus"), type="integer", default=16,
                    help="Num CPUs", metavar="NUMCPUS")
       
    )


# # Create parser and parse options
parser <- OptionParser(option_list=option_list)
options <- parse_args(parser)
# 
options_jobid <- options$jobid
options_numcpus <- options$numcpus

# 
 # options_jobid <- 1
 # options_numcpus <- 8

# 
# # Use the options
cat("Job Idx:", options_jobid, "\n")
cat("Num CPUs:", options_numcpus, "\n")


#################################################
#################################################



###########
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
    # n.cores <- parallel::detectCores()
    n.cores <- options_numcpus
    my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
    doParallel::registerDoParallel(cl = my.cluster)
    cat("Parellel Registered: ", foreach::getDoParRegistered(), " (num cores=", n.cores, ")\n")
    initialized_parallel <- TRUE
    
    # registerDoRNG(123) # ///<<<< THIS CREATES THE ERROR FOR FADPClust !!!
}

ensure_dir_exist <- function(directory_path){
    # Check if the directory exists
    if(!dir.exists(directory_path)) {
        # Directory doesn't exist, so create it
        dir.create(directory_path, recursive = TRUE)
        cat("Directory created:", directory_path, "\n")
    } else {
        cat("Directory already exists:", directory_path, "\n")
    }
}

# scenario_folder = "table_365_noextra"
# ensure_dir_exist(scenario_folder)

final_table_folder = paste0("365dayscensor1")
ensure_dir_exist(final_table_folder)

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
#Output: 3D array, n*5*B, n is the number of observation, 
#5 is the column of the data (outcome, group, outcome.t , outcome.b, event), 
#B is the number of replications
generate_home_time_scenario = function(B, n, censor, effect,  equalsize, diff_censor, censorbig, option_365) {
    # B=3
    # n=1000
    # censor=1
    # effect=0
    # equalsize=0
    # diff_censor=0
    # censorbig=0
    
    home_time_data=array(0,c(n,4,B))
    #same censoring rate
    if (censor==0 && diff_censor==0){censor_p=c(0)}
    if (censor==0 && diff_censor==1){censor_p=c(0)}
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
                #for (t in 1:90) {
                if (t > 1) {
                    prev.dead <- dead
                }
                
                lin <- -7 + .00015 * t + frailty - effect * (t < 14) * trt
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
                    #if (t == 90)
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
            
            # prob_dead_plot_data=matrix(0,nrow=400,ncol=n)
            # prob_out_plot_data=matrix(0,nrow=400,ncol=n)
            #boundary_date=120
            boundary_date=400
            prob_dead_plot_data=matrix(0,nrow=boundary_date,ncol=n)
            prob_out_plot_data=matrix(0,nrow=boundary_date,ncol=n)
            for (t in 1:boundary_date) {
                if (t > 1) {
                    prev.dead <- dead
                }
                
                lin <- -7 + .00015 * t + frailty - effect * (t < 14) * trt
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
                    #if (t == 90)
                    slice.365 <- cum.ao
                #' *still need a event indicator for Cox model, even if not censored*
                
            }
            # create censoring indicator
            data <- data.frame(outcome = slice.365, group = trt)
            data$outcome.t <- slice.365
            
            #####################Aug 30, uniform censoring
            if (diff_censor==0){
                uniform_cencoring=runif(n, 180, 365) #77%
                #uniform_cencoring=runif(n, 45, 90) #77%
                for (subject_index in 1:n){
                    ##add the proportion for the censoring, 70% for treatment and 35% for non-treatment
                    
                    if (data$outcome.t[subject_index]<uniform_cencoring[subject_index]){
                        data$htevent[subject_index] <- 1
                    }
                    
                    if (data$outcome.t[subject_index]>=uniform_cencoring[subject_index]){
                        data$htevent[subject_index] <- 0
                        data$outcome.t[subject_index]=uniform_cencoring[subject_index]
                    }
                    if (option_365==1){
                        if(data$outcome.t[subject_index]>=365){
                            data$htevent[subject_index] <- 0
                            data$outcome.t[subject_index]=365
                        }
                    }
                    
                    
                    censor_p[j]=1-sum(data$htevent)/n
                }
            }
            
            ###########################################################################
            #Nov 15, 2023 change the percentage of censoring based on treatment and non-treatment
            #########################################################################
            #0.1, 0.05, censor_prop_trt [1] 0.018 > censor_prop_nontrt [1] 0.004
            #0.4, 0.2,  0.09, 0.015,  0.3, 0.15., 0.046, 0.01
            #0.8 , 0.4 , 20%, 4%
            if (diff_censor==1) {
                data_nontrt=data[data$group==0,]
                data_trt=data[data$group==1,]
                
                if (censorbig==1){
                    n_trt=sample(1:dim(data_trt)[1],round(dim(data_trt)[1]*0.8))
                    n_nontrt=sample(1:dim(data_nontrt)[1],round(dim(data_nontrt)[1]*0.4))
                }
                
                if (censorbig==0){
                    n_trt=sample(1:dim(data_trt)[1],round(dim(data_trt)[1]*0.4))
                    n_nontrt=sample(1:dim(data_nontrt)[1],round(dim(data_nontrt)[1]*0.2))
                }
                
                if (censorbig==3){
                    n_trt=sample(1:dim(data_trt)[1],round(dim(data_trt)[1]*0.4))
                    n_nontrt=sample(1:dim(data_nontrt)[1],round(dim(data_nontrt)[1]*0.8))
                }
                
                if (censorbig==4){
                    n_trt=sample(1:dim(data_trt)[1],round(dim(data_trt)[1]*0.2))
                    n_nontrt=sample(1:dim(data_nontrt)[1],round(dim(data_nontrt)[1]*0.4))
                }
                
                #non censored hevent 1
                data$htevent=1
                ######treatment group
                ##censoring
                #uniform_cencoring_trt=runif(n_trt, 180, 365)
                for (subject_index in n_trt){
                    #uniform_cencoring_trt=runif(n_trt, 240, 365)
                    #10%
                    uniform_cencoring_trt=runif(length(n_trt), 1, 365)
                    #uniform_cencoring_trt=runif(length(n_trt), 1, 90)
                    ##add the proportion for the censoring, 70% for treatment and 50% for non-treatment
                    index_union = 1
                    if (data$outcome.t[subject_index]>uniform_cencoring_trt[index_union]){
                        data$htevent[subject_index] <- 0
                        data$outcome.t[subject_index]=uniform_cencoring_trt[index_union]
                    }
                    if (option_365==1){
                        if(data$outcome.t[subject_index]>=365){
                            data$htevent[subject_index] <- 0
                            data$outcome.t[subject_index]=365
                        }
                    }
                   
                    index_union=index_union+1
                    censor_prop_trt=(length(n_trt)-sum(data[n_trt,]$htevent))/n
                    
                }
                
                
                ###non-treatment group
                #uniform_cencoring_nontrt=runif(n_nontrt, 180, 365)
                for (subject_index in n_nontrt){
                    #5%
                    #uniform_cencoring_nontrt=runif(n_nontrt, 360,365 )
                    uniform_cencoring_nontrt=runif(length(n_nontrt), 1,365 )
                    #uniform_cencoring_nontrt=runif(length(n_nontrt), 1,90 )
                    ##add the proportion for the censoring, 70% for treatment and 50% for non-treatment
                    index_union = 1
                    if (data$outcome.t[subject_index]>uniform_cencoring_nontrt[index_union]){
                        data$htevent[subject_index] <- 0
                        data$outcome.t[subject_index]=uniform_cencoring_nontrt[index_union]
                    }
                    if (option_365==1){
                        if(data$outcome.t[subject_index]>=365){
                            data$htevent[subject_index] <- 0
                            data$outcome.t[subject_index]=365
                        }
                    }
                    index_union=index_union+1
                    censor_prop_nontrt=(length(n_nontrt)-sum(data[n_nontrt,]$htevent))/n
                }
                censor_p[1,,j]=c(1,censor_prop_trt)
                censor_p[2,,j]=c(0,censor_prop_nontrt)
            }
            
            
        }
        
        home_time_data[,,j]=as.matrix(data)
    }
    return(list("home_time_data"=home_time_data,"censor_p"=censor_p))
}

# 
# n500_diff_censor_no_effect=generate_home_time_scenario(B=2, n=500, censor=1, effect=0, 
#                                                     equalsize=0, diff_censor=1,censorbig=1)
# 
# start_time=Sys.time()
# n1000_same_censor_effect=generate_home_time_scenario(B=2, n=1000, censor=1, effect=1,
#                                   equalsize=0, diff_censor=0,censorbig=1)
# end_time=Sys.time()
# time_used=end_time-start_time
# time_used
# #save(n1000_same_censor_effect,file="n1000_same_censor_effect.RData")
# 
# start_time=Sys.time()
# n500_same_censor_effect=generate_home_time_scenario(B=5000, n=500, censor=1, effect=1, 
#                                                         equalsize=0, diff_censor=0,censorbig=1)
# end_time=Sys.time()
# time_used=end_time-start_time
# time_used
# 
# #save(n500_same_censor_effect,file="n500_same_censor_effect.RData")
# 
# start_time=Sys.time()
# n1000_diff_censor_effect=generate_home_time_scenario(B=5000, n=1000, censor=1, effect=1, 
#                                                         equalsize=0, diff_censor=1,censorbig=1)
# end_time=Sys.time()
# time_used=end_time-start_time
# time_used
# 
# #save(n1000_diff_censor_effect,file="n1000_diff_censor_effect.RData")
# 
# start_time=Sys.time()
# n500_diff_censor_effect=generate_home_time_scenario(B=5000, n=500, censor=1, effect=1, 
#                                                        equalsize=0, diff_censor=1,censorbig=1)
# end_time=Sys.time()
# time_used=end_time-start_time
# time_used
# save(n500_diff_censor_effect,n1000_diff_censor_effect,
#      n1000_same_censor_effect,n500_same_censor_effect,
#      n1000_same_censor_no_effect,n500_same_censor_no_effect,
#      n1000_diff_censor_no_effect,n500_diff_censor_no_effect,
#      file="n500_n1000_censor_effect.RData")

# apply(abcs$censor_p,c(1,2),mean)[1,2] #0.077
# apply(abcs$censor_p,c(1,2),mean)[2,2] #0.0178

#Function to find the results for simulated home time data
#Input : 3D array, n*5*B hometime data, n-observations, B: replications
#Output: average values (parameter, standard error for parameter and p-value) for B replications

home_time_simulation_results=function(home_time_scenario_data){
    num_replica=dim(home_time_scenario_data)[3]
    power <- matrix(nrow = num_replica, ncol = 6)
    param.est <- matrix(nrow = num_replica, ncol = 6)
    se.est <- matrix(nrow = num_replica, ncol = 6)
    n=dim(home_time_scenario_data)[1]
    
    foreach_out <- foreach (j = 1:num_replica, .combine = cbind, .init = NULL) %dorng% {
        source("./Source/R/functions_hometime.R")
        home_time_data_rep=data.frame(home_time_scenario_data[,,j])
        colnames(home_time_data_rep)=c("outcome", "group","outcome.t","htevent" )
        result_hometime=home_time_regressions(home_time_data_rep)
        param.est_j=result_hometime$param.est
        power_j=result_hometime$power
        se.est_j=result_hometime$se_model
        return(c(param.est_j, power_j,se.est_j))
    }
    
    param.est=t(foreach_out[1:6,])
    power=t(foreach_out[7:12,])
    se_model=t(foreach_out[13:18,])
    a.param <- apply(param.est, MARGIN = 2, FUN = mean)
    a.paramsd <- apply(param.est, MARGIN = 2, FUN = sd) 
    a.paramsdmodel <- apply(se_model, MARGIN = 2, FUN = mean) 
    a.power <- apply(power, MARGIN = 2, FUN = mean)
    a.powerse<- sqrt((a.power*(1-a.power))/num_replica)
    result <- cbind(a.param,a.paramsdmodel, a.paramsd, a.power, a.powerse)
    
    rownames(result) <- c("lin",  "med",  "poi","nb","cox", "temp")
    
    colnames(result) <-c("est", "semodel","se", "power","powerse")
    return(result)
}



#Function to produce the table for manuscript
#Input : B:scalar- the number of replications, n: scalar-the number of observations
#Output: average values (parameter, standard error for parameter and p-value) for B replications
# censor=0, no right censor, censor=1, yes ;
# equalsize=0, unequal allocation between treatment and control, equalsize=1, equal allocation
#effect: scalar, treatment effect, usually 0 means no effect, 1 means has some effect
#effect.d: scalar, effect of treatment on death, usually 0 means no effect of death, 1 means has some effect

#Output: 6 rows (one per model), 4 columns (parameter, standard error model, standard error, p-value)
home_time_table=function(B, n, censor, effect, equalsize,diff_censor,censorbig, option_365){
    home_time_data=generate_home_time_scenario(B, n, censor, effect, equalsize,diff_censor,censorbig,option_365)
    result_table=home_time_simulation_results(home_time_data$home_time_data)
    return(result_table)
}
#################################################################################################
##################################################################################################
#uncensored balanced type I
# set.seed(123)
# uncensor_balance_500 <- home_time_table(B = 5000, n = 500, censor = 0, effect = 0,  equalsize = 1,diff_censor=1)

begin_exp_time <- Sys.time()
    
#censor, effect,  equalsize, diff_censor, censorbig
get_file_name <- function(censor, effect,  equalsize, diff_censor, censorbig, option_365){
    part1 <- "TypeI"
    if(effect == 1){
        part1 <- "Power"
    }
    part2 <- "Uncensor"
    if(censor == 1){
        part2 <- "UniformCensor"
        if(diff_censor == 1) {
            part2 <- "SmallDiffCensor"
            if(censorbig == 1){
                part2 <- "BigDiffCensor"
            }else if(censorbig == 3){
                part2 <- "ReverseBigDiffCensor"
            }else if(censorbig == 4){
                part2 <- "ReverseSmallDiffCensor"
            }
        }
    }
    part3 <- "Balance"
    if(equalsize == 0){
        part3 <- "Unbalance"
    }
    part4 <- "NoCensor365"
    if(option_365 == 1){
        part4 <- "Censor365"
    }
    return(paste(sep="_",part1, part2, part3, part4))
}

run_exp_set <- function(censor, effect,  equalsize, diff_censor, censorbig, option_365){
    
    file_name <- get_file_name(censor, effect,  equalsize, diff_censor, censorbig, option_365)
    
    begin_exp_set_time <- Sys.time()
    set.seed(123)
    uncensor_balance_500 <- home_time_table(B = 5000, n =500, 
                                            censor = censor, 
                                            effect = effect,  
                                            equalsize = equalsize, 
                                            diff_censor = diff_censor, 
                                            censorbig = censorbig,
                                            option_365 = option_365)
    set.seed(123)
    uncensor_balance_1000 <- home_time_table(B = 5000, n = 1000, 
                                             censor = censor, 
                                             effect = effect,  
                                             equalsize = equalsize, 
                                             diff_censor = diff_censor, 
                                             censorbig = censorbig,
                                             option_365 = option_365)
    set.seed(123)
    uncensor_balance_5000 <- home_time_table(B = 5000, n = 5000, 
                                             censor = censor, 
                                             effect = effect,  
                                             equalsize = equalsize, 
                                             diff_censor = diff_censor, 
                                             censorbig = censorbig,
                                             option_365 = option_365)
    end_exp_set_time <- Sys.time()
    
    cat("\n",
        " ", file_name, " took : ", capture.output(end_exp_set_time - begin_exp_set_time), 
        "\n")
    
    typeI_censor_balance=cbind(uncensor_balance_500, uncensor_balance_1000, uncensor_balance_5000)
    save(typeI_censor_balance,file=paste0(final_table_folder, "/",file_name,".RData"))
}

# run_exp_set(0,0,1,1,1,1)
# run_exp_set(0,1,1,1,1,1)
# 
# run_exp_set(0,0,0,1,1,1)
# run_exp_set(0,1,0,1,1,1)
# 
# run_exp_set(1,0,1,0,1,1)
# run_exp_set(1,1,1,0,1,1)
# 
# run_exp_set(1,0,0,0,1,1)
# run_exp_set(1,1,0,0,1,1)
# 
# run_exp_set(1,0,1,1,1,1)
# run_exp_set(1,1,1,1,1,1)
# 
# run_exp_set(1,0,0,1,1,1)
# run_exp_set(1,1,0,1,1,1)

run_exp_set(1,0,1,1,0,1)
run_exp_set(1,1,1,1,0,1)

run_exp_set(1,0,0,1,0,1)
run_exp_set(1,1,0,1,0,1)

run_exp_set(1,0,1,1,3,1)
run_exp_set(1,1,1,1,3,1)

run_exp_set(1,0,0,1,3,1)
run_exp_set(1,1,0,1,3,1)

run_exp_set(1,0,1,1,4,1)
run_exp_set(1,1,1,1,4,1)

run_exp_set(1,0,0,1,4,1)
run_exp_set(1,1,0,1,4,1)

end_exp_time <- Sys.time()

cat("\n================\n",
    "\t All Experiments took : ", capture.output(end_exp_time - begin_exp_time), 
    "\n================\n")


if(run_parallel)
{
    parallel::stopCluster(cl = my.cluster)
    initialized_parallel <- FALSE
}
