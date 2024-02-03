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

############################################################
# Function to use Method of Zhan and Schaubel.
# Thanks to Jianghao Li for translating from article into code and extending
# Thanks to Sean O'Brien to provide the code
#######
# function to estimate
estimate.daoh <- function(AOU, U) {
  ### AOU (N x M). 1 row per patient. 1 column per day. Entry indicates whether patient was
  ###              (A) alive and (O) out-of-hospital and (U) uncensored on a given day
  ### U   (N x M). Entry indicates whether patient was (U) uncensored on a given day.
  ### Calculation assumes censoring times in absence of death are known e.g. administrative censoring only.
  ### For ISCHEMIA, we censored patients who died on their maximum last potential f/u (=6/30/2019).
  ### In simulation code below, I think we do something even fudgier.
  sumU <- colSums(U)
  PI <- colSums(AOU) / sumU
  D <- t((t(AOU) - t(U) * PI) / sumU)
  cum.est <- cumsum(PI)
  cum.se <- sqrt(colSums(t(apply(D, 1, cumsum)) ^ 2))
  data.frame(est = last(cum.est), se = last(cum.se))
}


#####
# Function to summarize results: parameter estimate, standard error and typeI error rate/Power
# Input: result from the simulation
#Out put: parameter estimate, standard error and typeI error rate/Power
# analysis_final_result <- function(result) {
#   n <- dim(result$power)[1]
#   a.param <- apply(result$param.est, MARGIN = 2, FUN = mean)
#   a.paramsd <- apply(result$param.est, MARGIN = 2, FUN = sd) / sqrt(n)
#   a.power <- apply(result$power, MARGIN = 2, FUN = mean)
#   result <- cbind(a.param, a.paramsd, a.power)
#   rownames(result) <- c("lin", "cox", "med", "nb", "poi", "temp")
#   colnames(result) <- c("est", "se", "power")
#   return(result)
# }


###########################################################
#Function to generate hoemtime data based on different scenarios
#Input:
#B: scalar, the number of repetition
#n: scalar, the number of total sample size
# censor=0, no right censor, censor=1, yes ;
# equalsize=0, unequal allocation between treatment and control, equalsize=1, equal allocation
#effect: scalar, treatment effect, usually 0 means no effect, 1 means has some effect
#effect.d: scalar, effect of treatment on death, usually 0 means no effect of death, 1 means has some effect

#Output: 3D array, n*5*B, n is the number of observation, 
#5 is the column of the data (outcome, group, outcome.t , outcome.b, event), 
#B is the number of replications
generate_home_time_scenario = function(B, n, censor, effect, effect.d, equalsize,diff_censor) {
  # B=3
  # n=1000
  # censor=1
  # effect=0
  # effect.d=0
  # equalsize=1
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
        uniform_cencoring_trt=runif(n_trt, 180, 365)
        for (subject_index in 1:n_trt){
          ##add the proportion for the censoring, 70% for treatment and 35% for non-treatment

          if (data_trt$outcome.t[subject_index]<uniform_cencoring_trt[subject_index]){
            data_trt$htevent[subject_index] <- 1
          }

          if (data_trt$outcome.t[subject_index]>=uniform_cencoring_trt[subject_index]){
            data_trt$htevent[subject_index] <- 0
            data_trt$outcome.t[subject_index]=uniform_cencoring_trt[subject_index]
          }
          censor_prop_trt=1-sum(data_trt$htevent)/n

        }

        
        #trt 0.7, non-trt 0.35
        # repeat{
        #   uniform_cencoring_trt=runif(n_trt, 180, 365)
        #   for (subject_index in 1:n_trt){
        #     ##add the proportion for the censoring, 70% for treatment and 35% for non-treatment
        #     
        #     if (data_trt$outcome.t[subject_index]<uniform_cencoring_trt[subject_index]){
        #       data_trt$htevent[subject_index] <- 1
        #     }
        #     
        #     if (data_trt$outcome.t[subject_index]>=uniform_cencoring_trt[subject_index]){
        #       data_trt$htevent[subject_index] <- 0
        #       data_trt$outcome.t[subject_index]=uniform_cencoring_trt[subject_index]
        #     }
        #     censor_prop_trt=1-sum(data_trt$htevent)/n
        #     
        #   }
        #   if (censor_prop_trt<0.8 && censor_prop_trt>0.7){break}
        # }
        
        
        ###non-treatment group
        uniform_cencoring_nontrt=runif(n_nontrt, 180, 365)
        for (subject_index in 1:n_nontrt){
          ##add the proportion for the censoring, 70% for treatment and 35% for non-treatment
          
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
     
      ###############################################################################
      ###############################################################################
      #############
      
      # if (data$outcome.t >= 365) {
      #   data$htevent <-
      #     0 #I think for censoring we need to force them to have an event at t=365?
      # } else {
      #   data$htevent <- 1
      # }
      # data$outcome.t[data$outcome.t > 365] <-
      #   365   #sets any HT beyond 365 to 365
    }
    #############figure1 in the manuscript
    # library(ggplot2)
    # data_admixture=data.frame(admixture,as.factor(trt))
    # colnames(data_admixture)=c("Admixture","Treatment")
    # admix_hist<-ggplot(data_admixture, aes(x=admixture, fill=Treatment, color=Treatment)) +
    #   geom_histogram(position="dodge", alpha=0.5)+
    #   theme(text=element_text(size=20))
    # 
    # data_frailty=data.frame(frailty,as.factor(trt))
    # colnames(data_frailty)=c("frailty","Treatment")
    # frailty_box=ggplot(data_frailty,aes(Treatment,frailty))+geom_boxplot(aes(fill=Treatment))+
    #   theme(text = element_text(size = 20))
    # 
    # library(gridExtra)
    # grid.arrange(admix_hist,frailty_box,  ncol = 2)
    ###################################
    
    ###Fig2 in the manuscript 1000*400
    ####
    # pd_plot_data=data.frame(t(prob_dead_plot_data),trt)
    # pd_plot_data_trt=pd_plot_data[pd_plot_data$trt==1,-401]
    # pd_plot_data_notrt=pd_plot_data[pd_plot_data$trt==0,-401]
    # mean_pd=apply(t(as.matrix(pd_plot_data_trt)),1,mean)
    # mean_pd_notrt=apply(t(as.matrix(pd_plot_data_notrt)),1,mean)
    # ###
    # prob_out_data=data.frame(t(prob_out_plot_data),trt)
    # pout_data_trt=prob_out_data[prob_out_data$trt==1,-401]
    # pout_data_notrt=prob_out_data[prob_out_data$trt==0,-401]
    # mean_pout=apply(t(as.matrix(pout_data_trt)),1,mean)
    # mean_pout_notrt=apply(t(as.matrix(pout_data_notrt)),1,mean)
    ###
    
    ###jsut the mean
    # par(mfrow=c(1,2))
    # plot(1:400,mean_pd,col="green",type = "l",ylim=c(0,0.01),xlab="Day",ylab="Probability of Death",cex.lab=2, cex.axis=2,lty=1, lwd=1)
    # lines(1:400,mean_pd_notrt,col="red",type = "l",cex.lab=2, cex.axis=2,lty=1, lwd=1)
    # legend("topleft",legend = c("Treatment","No Treatment"),col=c("green","red"),lty=c(1,1),lwd=3)
    # 
    # 
    # plot(1:400,mean_pout,col="green",type = "l",ylim=c(0,0.7),xlab="Day",ylab="Probability of Out",cex.lab=2, cex.axis=2,lty=1, lwd=1)
    # lines(1:400,mean_pout_notrt,col="red",type = "l",cex.lab=2, cex.axis=2,lty=1, lwd=1)
    # #legend("topleft",legend = c("Treatment","No Treatment"),col=c("green","red"),lty=c(1,1),lwd=3)
    # 
    # ####all individuals
    # par(mfrow=c(1,2))
    # matplot(1:400,t(as.matrix(pout_data_trt)),col="green",type = "l",ylim=c(0,0.1),xlab="Day",ylab="Probability of Death",cex.lab=2, cex.axis=2,lty=1, lwd=1)
    # matlines(1:400,t(as.matrix(pout_data_notrt)),col="red",type = "l",cex.lab=2, cex.axis=2,lty=1, lwd=1)
    # #legend("topleft",legend = c("Treatment","No Treatment"),col=c("green","red"),lty=c(1,1),lwd=3)
    # 
    # matplot(1:400,t(as.matrix(pd_plot_data_trt)),col="green",type = "l",ylim=c(0,0.1),xlab="Day",ylab="Probability of Out",cex.lab=2, cex.axis=2,lty=1, lwd=1)
    # matlines(1:400,t(as.matrix(pd_plot_data_notrt)),col="red",type = "l",cex.lab=2, cex.axis=2,lty=1, lwd=1)
    # legend("topleft",legend = c("Treatment","No Treatment"),col=c("green","red"),lty=c(1,1),lwd=3)
    
    
    home_time_data[,,j]=as.matrix(data)
  }
  return(list("home_time_data"=home_time_data,"censor_p"=censor_p))
}
#sample data abc
#generate_home_time_scenario = function(B, n, censor, effect, effect.d, equalsize)
abc=generate_home_time_scenario(5000, 1000, 1, 0, 0, 1,0)
#mean(abc$censor_p) #0.488444
# different censoring
#abcdiff=generate_home_time_scenario(5000, 1000, 1, 0, 0, 1,1)
#colnames(abcdiff$censor_p)=c("proportion","treatment")
#trt_p
#apply(abcdiff$censor_p,c(1,2),mean)[1,1]
#non_trt_p
#apply(abcdiff$censor_p,c(1,2),mean)[2,1]

#sample data abc
#generate_home_time_scenario = function(B, n, censor, effect, effect.d, equalsize,diff_censor=0)
#abc=generate_home_time_scenario(5000, 1000, 1, 0, 0, 1)
#mean(abc$censor_p) #0.488444



#Function to test 6 different models
#Input: 2D array, n*5, n: the number of the subjects, 
#                      columns:outcome group outcome.t  outcome.b event
#Output: list including two items: parameter (6 elements, one for each model) and 
#                                  power (6 elements, one for each model)
home_time_regressions = function(home_time_data) {
  #home_time_data=data.frame(abc[,,2])
  #colnames(home_time_data)=c("outcome", "group","outcome.t","outcome.b" ,"htevent"  )
  #colnames(home_time_data)=c("outcome", "group","outcome.t","htevent"  )
  
  power <- c(0)
  param.est <- c(0)
  # analysis options
  ###############################################
  # t-test, linear regression
  m1 <- lm(outcome.t ~ group, data = home_time_data)
  co1 <- summary(m1)$coefficients[2, 1]
  p1 <- summary(m1)$coefficients[2, 4] < 0.05
  ##################################################
  # COX MODEL - htevent 1: non-censored, htevent 0: censored
  
  
  m2.2 <- coxph(Surv(home_time_data$outcome.t, home_time_data$htevent) ~ home_time_data$group)
  co2new <- unname(m2.2$coefficients)
  p2new <- (summary(m2.2))$coef[5] < 0.05
  ###############################################
  # median regression
  m3 <- rq(outcome.t ~ group, data = home_time_data)
  co3 <- summary(m3, se = "ker")[[3]][2, 1]
  p3 <- summary(m3, se = "ker")[[3]][2, 4] < 0.05
  ###############################################
  
  # negative binomial
  
  m4 <-
    summary(glm(outcome.t ~ group, family = negative.binomial(.3), data =  home_time_data))
  p4 <- m4$coefficients[2, 4] < 0.05
  co4 <- m4$coefficients[2, 1]
  
  
  # Poisson
  # ROBUST VARIANCE IS EXTREMELY IMPORTANT
  ###############################################
  m5 <- glm(outcome.t ~ group, family = poisson, data =  home_time_data)
  sand_vcov <- sandwich(m5)
  sand_se <- sqrt(diag(sand_vcov))
  robust_z <- m5$coef / sand_se
  robust_p <- 2 * pnorm(-abs(robust_z), 0, 1)
  p5 <- robust_p[2] < 0.05
  co5 <- m5$coef[2]
  
  # m5.5 <- summary(m5)
  # p5new <- m5.5$coef[2, 4] < 0.05
  # co5new <- m5.5$coef[2, 1]
  ################################################
  # temporal process regression
  
  data0 <-  home_time_data[ home_time_data$group == 0,]$outcome.t
  data1 <-  home_time_data[ home_time_data$group == 1,]$outcome.t
  # find maximum time for each group
  max0 <- max(data0)
  max1 <- max(data1)
  nsub0 <- length(data0)
  nsub1 <- length(data1)
  # consider htevent=1 as censored
  # assign all values uncensored
  U0 <- matrix(rep(1, nsub0 * max0), nrow = nsub0, ncol = max0)
  U1 <- matrix(rep(1, nsub1 * max1), nrow = nsub1, ncol = max1)
  data0f <- U0
  data1f <- U1
  U00 <- U0
  U11 <- U1
  # for each subject
  # if htevent=0 (censored), any value from that row in U should be 0 after that day
  for (i in 1:nsub0) {
    data0f[i,][data0[i]:max0] <- 0
    if (home_time_data[home_time_data$group == 0,]$htevent[i] == 0) {
      U00[i,][data0[i]:max0] <- 0
    }
  }
  
  for (i in 1:nsub1) {
    data1f[i,][data1[i]:max1] <- 0
    if (home_time_data[home_time_data$group == 1,]$htevent[i] == 0) {
      U11[i,][data1[i]:max1] <- 0
    }
  }
  
  tmp1 <- estimate.daoh(data0f, U00)
  tmp2 <- estimate.daoh(data1f, U11)
  diff <- tmp1$est - tmp2$est
  var.diff <- (tmp1$se ^ 2 + tmp2$se ^ 2)
  pval <- 1 - pchisq(diff ^ 2 / var.diff, df = 1)
  co6 <- diff
  p6 <- pval < 0.05
  ################################################
  #linear regression, cox, median, nb, possion, temporal process
  param.est <- c(co1, co2new, co3, co4, co5, co6)
  names(param.est)=c("lr","cox","median","nb","poission","tp")
  power<- c(p1, p2new, p3, p4, p5, p6)
  names(power)=c("lr","cox","median","nb","poission","tp")
  return(list("param.est "=param.est ,"power"=power))
}

#Function to find the results for simulated home time data
#Input : 3D array, n*5*B hometime data, n-observations, B: replications
#Output: average values (parameter, standard error for parameter and p-value) for B replications
home_time_simulation_results=function(home_time_scenario_data){
  #home_time_scenario_data=abc
  num_replica=dim(home_time_scenario_data)[3]
  power <- matrix(nrow = num_replica, ncol = 6)
  param.est <- matrix(nrow = num_replica, ncol = 6)
  n=dim(home_time_scenario_data)[1]
  for (j in 1:num_replica ){
    #j=1
    home_time_data_rep=data.frame(home_time_scenario_data[,,j])
    # colnames(home_time_data_rep)=c("outcome", "group","outcome.t",
    #                                "outcome.b" ,"htevent" )
    colnames(home_time_data_rep)=c("outcome", "group","outcome.t","htevent" )
    param.est[j, ]=home_time_regressions(home_time_data_rep)$param.est
    power[j, ]=home_time_regressions(home_time_data_rep)$power
  }
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
home_time_table=function(B, n, censor, effect, effect.d, equalsize){
  home_time_data=generate_home_time_scenario(B, n, censor, effect, effect.d, equalsize)
  result_table=home_time_simulation_results(home_time_data$home_time_data)
  return(result_table)
}


#example to create the hometime table for sample size=1000, 5 replications
#home_time_table(3, 1000, 0, 0, 0, 1)
set.seed(123)
uncensor_balance_500 <- home_time_table(B = 5000, n = 500, censor = 0, effect = 0,
                                        effect.d = 0, equalsize = 1)
set.seed(123)
uncensor_balance_1000  <- home_time_table(B = 5000, n = 1000, censor = 0, effect = 0, 
                                          effect.d = 0, equalsize = 1)
# # Second setting, effect = 0, don't want power
#uncensor_balance_1000 <- analysisf(result1000)


set.seed(123)
uncensor_balance_5000 <- home_time_table(B = 5000, n = 5000, censor = 0, effect = 0, effect.d = 0, equalsize = 1)
# # Second setting, effect = 0, don't want power


typeI_uncensor_balance=cbind(uncensor_balance_500,uncensor_balance_1000,uncensor_balance_5000)
save(typeI_uncensor_balance,file="typeI_uncensor_balance.RData")
load("typeI_uncensor_balance.RData")
xtable(typeI_uncensor_balance)
##################################################################################
#####power uncensor balanced
set.seed(123)
uncensor_balance_500 <- home_time_table(B = 5000, n = 500, censor = 0, effect = 1, effect.d = 0, equalsize = 1)
# # Second setting, effect = 0, don't want power


set.seed(123)
uncensor_balance_1000 <- home_time_table(B = 5000, n = 1000, censor = 0, effect = 1, effect.d = 0, equalsize = 1)
# # Second setting, effect = 0, don't want power

set.seed(123)
uncensor_balance_5000 <- home_time_table(B = 5000, n = 5000, censor = 0, effect = 1, effect.d = 0, equalsize = 1)
# # Second setting, effect = 0, don't want power


power_uncensor_balance=cbind(uncensor_balance_500,uncensor_balance_1000,uncensor_balance_5000)
save(power_uncensor_balance,file="power_uncensor_balance.RData")

###############################################################################
#####type I uncensor unbalanced
set.seed(123)
uncensor_balance_500 <- home_time_table(B = 5000, n = 500, censor = 0, effect = 0, effect.d = 0, equalsize = 0)
# # Second setting, effect = 0, don't want power
#uncensor_balance_500 <- analysisf(result500)

set.seed(123)
uncensor_balance_1000 <- home_time_table(B = 5000, n = 1000, censor = 0, effect = 0, effect.d = 0, equalsize = 0)
# # Second setting, effect = 0, don't want power
#uncensor_balance_1000 <- analysisf(result1000)


set.seed(123)
uncensor_balance_5000 <- home_time_table(B = 5000, n = 5000, censor = 0, effect = 0, effect.d = 0, equalsize = 0)
# # Second setting, effect = 0, don't want power
#uncensor_balance_5000 <- analysisf(result5000)

typeI_uncensor_unbalance=cbind(uncensor_balance_500,uncensor_balance_1000,uncensor_balance_5000)
save(typeI_uncensor_unbalance,file="typeI_uncensor_unbalance.RData")

##################################################################################
#####power uncensor unbalanced
set.seed(123)
uncensor_balance_500 <- home_time_table(B = 5000, n = 500, censor = 0, effect = 1, effect.d = 0, equalsize = 0)
# # Second setting, effect = 0, don't want power
#uncensor_balance_500 <- analysisf(result500)

set.seed(123)
uncensor_balance_1000 <- home_time_table(B = 5000, n = 1000, censor = 0, effect = 1, effect.d = 0, equalsize = 0)
# # Second setting, effect = 0, don't want power
#uncensor_balance_1000 <- analysisf(result1000)


set.seed(123)
uncensor_balance_5000 <- home_time_table(B = 5000, n = 5000, censor = 0, effect = 1, effect.d = 0, equalsize = 0)
# # Second setting, effect = 0, don't want power
#uncensor_balance_5000 <- analysisf(result5000)

power_uncensor_unbalance=cbind(uncensor_balance_500,uncensor_balance_1000,uncensor_balance_5000)
save(power_uncensor_unbalance,file="power_uncensor_unbalance.RData")


##################################################################################################
####################################################################################################
#censored balanced type I
set.seed(123)
uncensor_balance_500 <- home_time_table(B = 5000, n = 500, censor = 1, effect = 0, effect.d = 0, equalsize = 1)
# # Second setting, effect = 0, don't want power
#uncensor_balance_500 <- analysisf(result500)

set.seed(123)
uncensor_balance_1000 <- home_time_table(B = 5000, n = 1000, censor = 1, effect = 0, effect.d = 0, equalsize = 1)
# # Second setting, effect = 0, don't want power
#uncensor_balance_1000 <- analysisf(result1000)


set.seed(123)
uncensor_balance_5000 <- home_time_table(B = 5000, n = 5000, censor = 1, effect = 0, effect.d = 0, equalsize = 1)
# # Second setting, effect = 0, don't want power
#uncensor_balance_5000 <- analysisf(result5000)

typeI_censor_balance=cbind(uncensor_balance_500,uncensor_balance_1000,uncensor_balance_5000)
save(typeI_censor_balance,file="typeI_censor_balance.RData")

##################################################################################
#####power censor balanced
set.seed(123)
uncensor_balance_500 <- home_time_table(B = 5000, n = 500, censor = 1, effect = 1, effect.d = 0, equalsize = 1)
# # Second setting, effect = 0, don't want power
#uncensor_balance_500 <- analysisf(result500)

set.seed(123)
uncensor_balance_1000 <- home_time_table(B = 5000, n = 1000, censor = 1, effect = 1, effect.d = 0, equalsize = 1)
# # Second setting, effect = 0, don't want power
#uncensor_balance_1000 <- analysisf(result1000)


set.seed(123)
uncensor_balance_5000 <- home_time_table(B = 5000, n = 5000, censor = 1, effect = 1, effect.d = 0, equalsize = 1)
# # Second setting, effect = 0, don't want power
#uncensor_balance_5000 <- analysisf(result5000)

power_censor_balance=cbind(uncensor_balance_500,uncensor_balance_1000,uncensor_balance_5000)
save(power_censor_balance,file="power_censor_balance.RData")

###############################################################################
#####type I censor unbalanced
set.seed(123)
uncensor_balance_500 <- home_time_table(B = 5000, n = 500, censor = 1, effect = 0, effect.d = 0, equalsize = 0)
# # Second setting, effect = 0, don't want power
#uncensor_balance_500 <- analysisf(result500)

set.seed(123)
uncensor_balance_1000 <- home_time_table(B = 5000, n = 1000, censor = 1, effect = 0, effect.d = 0, equalsize = 0)
# # Second setting, effect = 0, don't want power
#uncensor_balance_1000 <- analysisf(result1000)


set.seed(123)
uncensor_balance_5000 <- home_time_table(B = 5000, n = 5000, censor = 1, effect = 0, effect.d = 0, equalsize = 0)
# # Second setting, effect = 0, don't want power
#uncensor_balance_5000 <- analysisf(result5000)

typeI_censor_unbalance=cbind(uncensor_balance_500,uncensor_balance_1000,uncensor_balance_5000)
save(typeI_censor_unbalance,file="typeI_censor_unbalance.RData")

##################################################################################
#####power censor unbalanced
set.seed(123)
uncensor_balance_500 <- home_time_table(B = 5000, n = 500, censor = 1, effect = 1, effect.d = 0, equalsize = 0)
# # Second setting, effect = 0, don't want power
#uncensor_balance_500 <- analysisf(result500)

set.seed(123)
uncensor_balance_1000 <- home_time_table(B = 5000, n = 1000, censor = 1, effect = 1, effect.d = 0, equalsize = 0)
# # Second setting, effect = 0, don't want power
#uncensor_balance_1000 <- analysisf(result1000)


set.seed(123)
uncensor_balance_5000 <- home_time_table(B = 5000, n = 5000, censor = 1, effect = 1, effect.d = 0, equalsize = 0)
# # Second setting, effect = 0, don't want power
#uncensor_balance_5000 <- analysisf(result5000)

power_censor_unbalance=cbind(uncensor_balance_500,uncensor_balance_1000,uncensor_balance_5000)
save(power_censor_unbalance,file="power_censor_unbalance.RData")

#' 
#' ########big function-old version
#' simfun <- function(B, n, censor, effect, effect.d, equalsize) {
#' 
#'  
#'   power <- matrix(nrow = B, ncol = 6)
#'   param.est <- matrix(nrow = B, ncol = 6)
#' 
#' 
#'   for (j in 1:B) {
#'     if (censor == 0) {
#'     
#'       admixture <- c(rep(3, n / 4), rep(0, n / 4), rep(0, n / 4), rep(-1, n / 4))
#' #' *what is purpose of admixture?*
#' #' 
#'       frailty <- rnorm(n, admixture, .5)
#'       
#'       
#'       if (equalsize == 1) {
#'         trt <- rbinom(n, 1, .5)
#'       }
#'       if (equalsize == 0) {
#'         trt <- rbinom(n, 1, .6)
#'       }
#' 
#'       cum.ao <- rep(0, n)
#'       prev.dead <- rep(0, n)
#' 
#' #' *I would make the 365 below into a variable that is set in the function call as well*
#' #'  *Comment what this loop is doing*     
#'       for (t in 1:365) {
#'         if (t > 1) {
#'           prev.dead <- dead
#'         }
#'         #' *Is there a specific reason for -7 as the slope and 0.00015 as the slope? What about the frailty?*
#'         lin <- -7 + .00015 * t + frailty - effect.d * (t < 14) * trt
#' #' *The probability of the death from the linear predictor*
#' #' 
#'         p.d <- exp(lin) / (1 + exp(lin))
#'         dead <- ((rbinom(n, 1, p.d) + prev.dead) > 0)
#' #' *The probability of out of the hospital from the linear predictor*
#'         lin <- 2 - 1.5 * frailty + effect * trt
#'         p.o <- exp(lin) / (1 + exp(lin))
#'         out <- rbinom(n, 1, p.o)
#' 
#'         alive.out <- (1 - dead) * out
#' 
#'         cum.ao <- cum.ao + alive.out
#'         if (t == 365) slice.365 <- cum.ao
#'       }
#'       data <- data.frame(outcome = slice.365, group = trt)
#'       data$outcome.t <- slice.365
#'       data$outcome.b <- slice.365 / 365
#'       
#'       # htevent 1: non censored, 0: censored
#'       data$htevent <- rep(1, n)
#'     }
#' 
#' 
#'     if (censor == 1) {
#'       admixture <- c(rep(3, n / 4), rep(0, n / 4), rep(0, n / 4), rep(-1, n / 4))
#'       frailty <- rnorm(n, admixture, .5)
#' 
#'       if (equalsize == 1) {
#'         trt <- rbinom(n, 1, .5)
#'       }
#'       if (equalsize == 0) {
#'         trt <- rbinom(n, 1, .6)
#'       }
#'       cum.ao <- rep(0, n)
#'       prev.dead <- rep(0, n)
#'       for (t in 1:400) {
#'         if (t > 1) {
#'           prev.dead <- dead
#'         }
#' 
#'         lin <- -7 + .00015 * t + frailty - effect.d * (t < 14) * trt
#'         p.d <- exp(lin) / (1 + exp(lin))
#'         dead <- ((rbinom(n, 1, p.d) + prev.dead) > 0)
#' 
#'         lin <- 2 - 1.5 * frailty + effect * trt
#'         p.o <- exp(lin) / (1 + exp(lin))
#'         out <- rbinom(n, 1, p.o)
#' 
#'         alive.out <- (1 - dead) * out
#' 
#'         cum.ao <- cum.ao + alive.out
#'         if (t == 365) slice.365 <- cum.ao   
#' #' *still need a event indicator for Cox model, even if not censored*
#' 
#'       }
#'     # create censoring indicator 
#'       data <- data.frame(outcome = slice.365, group = trt)
#'       data$outcome.t <- slice.365
#'       if (data$outcome.t >= 365) {
#'         data$htevent <- 0 #I think for censoring we need to force them to have an event at t=365? 
#'       } else {
#'         data$htevent <- 1
#'       }
#'       data$outcome.t[data$outcome.t > 365] <- 365   #sets any HT beyond 365 to 365
#'     }
#' 
#' 
#' #' *Consider splitting your code into small functions*: 
#' #' *data generation*
#' #' *model running and storage of results*
#' #' *summarization of results* 
#'     
#'     # analysis options
#'     ###############################################
#'     # t-test, linear regression
#'     m1 <- lm(outcome.t ~ group, data = data)
#'     co1 <- summary(m1)$coefficients[2, 1]
#'     p1 <- summary(m1)$coefficients[2, 4] < 0.05
#'     ##################################################
#'     # COX MODEL - htevent 1: non-censored, htevent 0: censored
#' 
#' 
#'     m2.2 <- coxph(Surv(data$outcome.t, data$htevent) ~ data$group)
#'     co2new <- unname(m2.2$coefficients)
#'     p2new <- (summary(m2.2))$coef[5] < 0.05
#'     ###############################################
#'     # median regression
#'     m3 <- rq(outcome.t ~ group, data = data)
#'     co3 <- summary(m3, se = "ker")[[3]][2, 1]
#'     p3 <- summary(m3, se = "ker")[[3]][2, 4] < 0.05
#'     ###############################################
#' 
#'     # negative binomial
#' 
#'     m4 <- summary(glm(outcome.t ~ group, family = negative.binomial(.3), data = data))
#'     p4 <- m4$coefficients[2, 4] < 0.05
#'     co4 <- m4$coefficients[2, 1]
#' 
#' 
#'     # Poisson
#'     # ROBUST VARIANCE IS EXTREMELY IMPORTANT
#'     ###############################################
#'     m5 <- glm(outcome.t ~ group, family = poisson, data = data)
#'     sand_vcov <- sandwich(m5)
#'     sand_se <- sqrt(diag(sand_vcov))
#'     robust_z <- m5$coef / sand_se
#'     robust_p <- 2 * pnorm(-abs(robust_z), 0, 1)
#'     p5 <- robust_p[2] < 0.05
#'     co5 <- m5$coef[2]
#' 
#'     m5.5 <- summary(m5)
#'     p5new <- m5.5$coef[2, 4] < 0.05
#'     co5new <- m5.5$coef[2, 1]
#'     ################################################
#'     # temporal process regression
#' 
#'     data0 <- data[data$group == 0, ]$outcome.t
#'     data1 <- data[data$group == 1, ]$outcome.t
#'     # find maximum time for each group
#'     max0 <- max(data0)
#'     max1 <- max(data1)
#'     nsub0 <- length(data0)
#'     nsub1 <- length(data1)
#'     # consider htevent=1 as censored
#'     # assign all values uncensored
#'     U0 <- matrix(rep(1, nsub0 * max0), nrow = nsub0, ncol = max0)
#'     U1 <- matrix(rep(1, nsub1 * max1), nrow = nsub1, ncol = max1)
#'     data0f <- U0
#'     data1f <- U1
#'     U00 <- U0
#'     U11 <- U1
#'     # for each subject
#'     # if htevent=0 (censored), any value from that row in U should be 0 after that day
#'     for (i in 1:nsub0) {
#'       data0f[i, ][data0[i]:max0] <- 0
#'       if (data[data$group == 0, ]$htevent[i] == 0) {
#'         U00[i, ][data0[i]:max0] <- 0
#'       }
#'     }
#' 
#'     for (i in 1:nsub1) {
#'       data1f[i, ][data1[i]:max1] <- 0
#'       if (data[data$group == 1, ]$htevent[i] == 0) {
#'         U11[i, ][data1[i]:max1] <- 0
#'       }
#'     }
#' 
#'     tmp1 <- estimate.daoh(data0f, U00)
#'     tmp2 <- estimate.daoh(data1f, U11)
#'     diff <- tmp1$est - tmp2$est
#'     var.diff <- (tmp1$se^2 + tmp2$se^2)
#'     pval <- 1 - pchisq(diff^2 / var.diff, df = 1)
#'     co6 <- diff
#'     p6 <- pval < 0.05
#'     ################################################
#' 
#'     param.est[j, ] <- c(co1, co2new, co3, co4, co5new, co6)
#'     power[j, ] <- c(p1, p2new, p3, p4, p5new, p6)
#'   }
#'   output <- list(power = power, param.est = param.est)
#'   return(output)
#' }
#' 
#' 
#' 




#############################################################################################
# simfunc <- function(B, n, censor, effect, effect.d,equalsize)
#####type I uncensor balanced
#result500 <- simfun(B = 5000, n = 500, censor = 0, effect = 0, effect.d = 0, equalsize = 1)
# # Second setting, effect = 0, don't want power
#home_time_table=function(B, n, censor, effect, effect.d, equalsize)
#uncensor_balance_500 <- analysisf(result500)
# uncensor_balance_500 <- home_time_table(B = 5000, n = 500, censor = 0, effect = 0,
#                                         effect.d = 0, equalsize = 1)
# 
# result1000 <- simfun(B = 5000, n = 1000, censor = 0, effect = 0, effect.d = 0, equalsize = 1)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_1000 <- analysisf(result1000)
# 
# 
# 
# result5000 <- simfun(B = 5000, n = 5000, censor = 0, effect = 0, effect.d = 0, equalsize = 1)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_5000 <- analysisf(result5000)
# 
# typeI_uncensor_balance=cbind(uncensor_balance_500,uncensor_balance_1000,uncensor_balance_5000)
# save(typeI_uncensor_balance,file="typeI_uncensor_balance.RData")
# load("typeI_uncensor_balance.RData")
# xtable(typeI_uncensor_balance)
# % latex table generated in R 4.2.1 by xtable 1.8-4 package
# % Mon Jun 26 11:13:33 2023
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrrr}
# \hline
# & a.param & a.paramsd & a.power & a.param & a.paramsd & a.power & a.param & a.paramsd & a.power \\ 
# \hline
# lin & 0.03 & 0.18 & 0.05 & 0.20 & 0.13 & 0.05 & -0.03 & 0.06 & 0.05 \\ 
# cox & 0.00 & 0.00 & 0.05 & -0.00 & 0.00 & 0.06 & -0.00 & 0.00 & 0.05 \\ 
# med & 0.04 & 0.22 & 0.03 & 0.31 & 0.15 & 0.05 & -0.03 & 0.07 & 0.11 \\ 
# nb & 0.00 & 0.00 & 0.05 & 0.00 & 0.00 & 0.05 & -0.00 & 0.00 & 0.05 \\ 
# poi & 0.00 & 0.00 & 0.84 & 0.00 & 0.00 & 0.84 & -0.00 & 0.00 & 0.84 \\ 
# temp & -0.03 & 0.18 & 0.05 & -0.20 & 0.13 & 0.05 & 0.03 & 0.06 & 0.05 \\ 
# \hline
# \end{tabular}
# \end{table}
##################################################################################
#####power uncensor balanced
# result500 <- simfun(B = 5000, n = 500, censor = 0, effect = 1, effect.d = 0, equalsize = 1)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_500 <- analysisf(result500)
# 
# 
# result1000 <- simfun(B = 5000, n = 1000, censor = 0, effect = 1, effect.d = 0, equalsize = 1)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_1000 <- analysisf(result1000)
# 
# 
# 
# result5000 <- simfun(B = 5000, n = 5000, censor = 0, effect = 1, effect.d = 0, equalsize = 1)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_5000 <- analysisf(result5000)
# 
# power_uncensor_balance=cbind(uncensor_balance_500,uncensor_balance_1000,uncensor_balance_5000)
# save(power_uncensor_balance,file="power_uncensor_balance.RData")
# load("power_uncensor_balance.RData")
# % latex table generated in R 4.2.1 by xtable 1.8-4 package
# % Mon Jun 26 11:12:41 2023
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrrr}
# \hline
# & a.param & a.paramsd & a.power & a.param & a.paramsd & a.power & a.param & a.paramsd & a.power \\ 
# \hline
# lin & 15.92 & 0.19 & 0.22 & 15.68 & 0.13 & 0.38 & 15.79 & 0.06 & 0.96 \\ 
# cox & -0.47 & 0.00 & 1.00 & -0.47 & 0.00 & 1.00 & -0.47 & 0.00 & 1.00 \\ 
# med & 35.82 & 0.19 & 0.54 & 36.21 & 0.12 & 0.97 & 36.46 & 0.05 & 1.00 \\ 
# nb & 0.07 & 0.00 & 0.21 & 0.07 & 0.00 & 0.38 & 0.07 & 0.00 & 0.96 \\ 
# poi & 0.07 & 0.00 & 0.92 & 0.07 & 0.00 & 0.96 & 0.07 & 0.00 & 1.00 \\ 
# temp & -15.90 & 0.19 & 0.22 & -15.65 & 0.13 & 0.38 & -15.77 & 0.06 & 0.96 \\ 
# \hline
# \end{tabular}
# \end{table}
###############################################################################
#####type I uncensor unbalanced
# result500 <- simfun(B = 5000, n = 500, censor = 0, effect = 0, effect.d = 0, equalsize = 0)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_500 <- analysisf(result500)
# 
# 
# result1000 <- simfun(B = 5000, n = 1000, censor = 0, effect = 0, effect.d = 0, equalsize = 0)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_1000 <- analysisf(result1000)
# 
# 
# 
# result5000 <- simfun(B = 5000, n = 5000, censor = 0, effect = 0, effect.d = 0, equalsize = 0)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_5000 <- analysisf(result5000)
# 
# typeI_uncensor_unbalance=cbind(uncensor_balance_500,uncensor_balance_1000,uncensor_balance_5000)
# save(typeI_uncensor_unbalance,file="typeI_uncensor_unbalance.RData")
# load("typeI_uncensor_unbalance.RData")
# xtable(typeI_uncensor_unbalance)
# % latex table generated in R 4.2.1 by xtable 1.8-4 package
# % Mon Jun 26 11:14:28 2023
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrrr}
# \hline
# & a.param & a.paramsd & a.power & a.param & a.paramsd & a.power & a.param & a.paramsd & a.power \\ 
# \hline
# lin & -0.12 & 0.19 & 0.05 & 0.42 & 0.13 & 0.05 & -0.00 & 0.06 & 0.05 \\ 
# cox & -0.00 & 0.00 & 0.05 & -0.00 & 0.00 & 0.05 & -0.00 & 0.00 & 0.05 \\ 
# med & 0.79 & 0.23 & 0.03 & 0.82 & 0.16 & 0.05 & 0.06 & 0.07 & 0.12 \\ 
# nb & -0.00 & 0.00 & 0.05 & 0.00 & 0.00 & 0.05 & 0.00 & 0.00 & 0.05 \\ 
# poi & -0.00 & 0.00 & 0.84 & 0.00 & 0.00 & 0.84 & 0.00 & 0.00 & 0.84 \\ 
# temp & 0.12 & 0.19 & 0.05 & -0.42 & 0.13 & 0.05 & 0.00 & 0.06 & 0.05 \\ 
# \hline
# \end{tabular}
# \end{table}
##################################################################################
# #####power uncensor unbalanced
# result500 <- simfun(B = 5000, n = 500, censor = 0, effect = 1, effect.d = 0, equalsize = 0)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_500 <- analysisf(result500)
# 
# 
# result1000 <- simfun(B = 5000, n = 1000, censor = 0, effect = 1, effect.d = 0, equalsize = 0)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_1000 <- analysisf(result1000)
# 
# 
# 
# result5000 <- simfun(B = 5000, n = 5000, censor = 0, effect = 1, effect.d = 0, equalsize = 0)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_5000 <- analysisf(result5000)
# 
# power_uncensor_unbalance=cbind(uncensor_balance_500,uncensor_balance_1000,uncensor_balance_5000)
# save(power_uncensor_unbalance,file="power_uncensor_unbalance.RData")
# load("power_uncensor_unbalance.RData")
# xtable(power_uncensor_unbalance)
# % latex table generated in R 4.2.1 by xtable 1.8-4 package
# % Mon Jun 26 11:15:19 2023
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrrr}
# \hline
# & a.param & a.paramsd & a.power & a.param & a.paramsd & a.power & a.param & a.paramsd & a.power \\ 
# \hline
# lin & 16.05 & 0.19 & 0.21 & 15.77 & 0.14 & 0.37 & 15.80 & 0.06 & 0.96 \\ 
# cox & -0.47 & 0.00 & 1.00 & -0.47 & 0.00 & 1.00 & -0.47 & 0.00 & 1.00 \\ 
# med & 37.01 & 0.20 & 0.52 & 36.58 & 0.14 & 0.97 & 36.48 & 0.06 & 1.00 \\ 
# nb & 0.07 & 0.00 & 0.22 & 0.07 & 0.00 & 0.38 & 0.07 & 0.00 & 0.96 \\ 
# poi & 0.07 & 0.00 & 0.92 & 0.07 & 0.00 & 0.96 & 0.07 & 0.00 & 1.00 \\ 
# temp & -16.02 & 0.19 & 0.22 & -15.74 & 0.14 & 0.38 & -15.77 & 0.06 & 0.96 \\ 
# \hline
# \end{tabular}
# \end{table}


##################################################################################################
####################################################################################################
#censored balanced type I
# result500 <- simfun(B = 5000, n = 500, censor = 1, effect = 0, effect.d = 0, equalsize = 1)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_500 <- analysisf(result500)
# 
# 
# result1000 <- simfun(B = 5000, n = 1000, censor = 1, effect = 0, effect.d = 0, equalsize = 1)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_1000 <- analysisf(result1000)
# 
# 
# 
# result5000 <- simfun(B = 5000, n = 5000, censor = 1, effect = 0, effect.d = 0, equalsize = 1)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_5000 <- analysisf(result5000)
# 
# typeI_censor_balance=cbind(uncensor_balance_500,uncensor_balance_1000,uncensor_balance_5000)
# save(typeI_censor_balance,file="typeI_censor_balance.RData")
# # load("typeI_censor_balance.RData")
# # xtable(typeI_censor_balance)
# # % latex table generated in R 4.2.1 by xtable 1.8-4 package
# # % Mon Jun 26 11:17:05 2023
# # \begin{table}[ht]
# # \centering
# # \begin{tabular}{rrrrrrrrrr}
# # \hline
# # & a.param & a.paramsd & a.power & a.param & a.paramsd & a.power & a.param & a.paramsd & a.power \\ 
# # \hline
# # lin & 0.08 & 0.19 & 0.05 & -0.09 & 0.13 & 0.05 & 0.00 & 0.06 & 0.05 \\ 
# # cox & -0.00 & 0.00 & 0.05 & 0.00 & 0.00 & 0.05 & -0.00 & 0.00 & 0.05 \\ 
# # med & 0.19 & 0.22 & 0.03 & 0.01 & 0.15 & 0.04 & 0.04 & 0.07 & 0.11 \\ 
# # nb & 0.00 & 0.00 & 0.05 & -0.00 & 0.00 & 0.05 & 0.00 & 0.00 & 0.05 \\ 
# # poi & 0.00 & 0.00 & 0.85 & -0.00 & 0.00 & 0.85 & 0.00 & 0.00 & 0.84 \\ 
# # temp & -0.08 & 0.19 & 0.06 & 0.09 & 0.13 & 0.05 & -0.00 & 0.06 & 0.05 \\ 
# # \hline
# # \end{tabular}
# # \end{table}
# ##################################################################################
# #####power censor balanced
# result500 <- simfun(B = 5000, n = 500, censor = 1, effect = 1, effect.d = 0, equalsize = 1)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_500 <- analysisf(result500)
# 
# 
# result1000 <- simfun(B = 5000, n = 1000, censor = 1, effect = 1, effect.d = 0, equalsize = 1)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_1000 <- analysisf(result1000)
# 
# 
# 
# result5000 <- simfun(B = 5000, n = 5000, censor = 1, effect = 1, effect.d = 0, equalsize = 1)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_5000 <- analysisf(result5000)
# 
# power_censor_balance=cbind(uncensor_balance_500,uncensor_balance_1000,uncensor_balance_5000)
# save(power_censor_balance,file="power_censor_balance.RData")
# 
# # load("power_censor_balance.RData")
# # xtable(power_censor_balance)
# # % latex table generated in R 4.2.1 by xtable 1.8-4 package
# # % Mon Jun 26 11:17:45 2023
# # \begin{table}[ht]
# # \centering
# # \begin{tabular}{rrrrrrrrrr}
# # \hline
# # & a.param & a.paramsd & a.power & a.param & a.paramsd & a.power & a.param & a.paramsd & a.power \\ 
# # \hline
# # lin & 15.98 & 0.19 & 0.22 & 15.81 & 0.13 & 0.39 & 15.77 & 0.06 & 0.97 \\ 
# # cox & -0.47 & 0.00 & 1.00 & -0.47 & 0.00 & 1.00 & -0.47 & 0.00 & 1.00 \\ 
# # med & 36.06 & 0.19 & 0.54 & 36.40 & 0.12 & 0.97 & 36.42 & 0.05 & 1.00 \\ 
# # nb & 0.07 & 0.00 & 0.22 & 0.07 & 0.00 & 0.39 & 0.07 & 0.00 & 0.97 \\ 
# # poi & 0.07 & 0.00 & 0.93 & 0.07 & 0.00 & 0.96 & 0.07 & 0.00 & 1.00 \\ 
# # temp & -15.95 & 0.19 & 0.22 & -15.78 & 0.13 & 0.39 & -15.75 & 0.06 & 0.97 \\ 
# # \hline
# # \end{tabular}
# # \end{table}
# ###############################################################################
# #####type I censor unbalanced
# result500 <- simfun(B = 5000, n = 500, censor = 1, effect = 0, effect.d = 0, equalsize = 0)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_500 <- analysisf(result500)
# 
# 
# result1000 <- simfun(B = 5000, n = 1000, censor = 1, effect = 0, effect.d = 0, equalsize = 0)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_1000 <- analysisf(result1000)
# 
# 
# 
# result5000 <- simfun(B = 5000, n = 5000, censor = 1, effect = 0, effect.d = 0, equalsize = 0)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_5000 <- analysisf(result5000)
# 
# typeI_censor_unbalance=cbind(uncensor_balance_500,uncensor_balance_1000,uncensor_balance_5000)
# save(typeI_censor_unbalance,file="typeI_censor_unbalance.RData")
# 
# # load("typeI_censor_unbalance.RData")
# # xtable(typeI_censor_unbalance)
# # 
# # % latex table generated in R 4.2.1 by xtable 1.8-4 package
# # % Mon Jun 26 11:18:34 2023
# # \begin{table}[ht]
# # \centering
# # \begin{tabular}{rrrrrrrrrr}
# # \hline
# # & a.param & a.paramsd & a.power & a.param & a.paramsd & a.power & a.param & a.paramsd & a.power \\ 
# # \hline
# # lin & -0.18 & 0.19 & 0.05 & 0.08 & 0.14 & 0.05 & 0.06 & 0.06 & 0.05 \\ 
# # cox & -0.00 & 0.00 & 0.05 & -0.00 & 0.00 & 0.05 & 0.00 & 0.00 & 0.05 \\ 
# # med & 0.47 & 0.23 & 0.03 & 0.37 & 0.16 & 0.04 & 0.04 & 0.07 & 0.12 \\ 
# # nb & -0.00 & 0.00 & 0.05 & 0.00 & 0.00 & 0.05 & 0.00 & 0.00 & 0.05 \\ 
# # poi & -0.00 & 0.00 & 0.84 & 0.00 & 0.00 & 0.84 & 0.00 & 0.00 & 0.84 \\ 
# # temp & 0.18 & 0.19 & 0.05 & -0.08 & 0.14 & 0.05 & -0.06 & 0.06 & 0.05 \\ 
# # \hline
# # \end{tabular}
# # \end{table}
# ##################################################################################
# #####power censor unbalanced
# result500 <- simfun(B = 5000, n = 500, censor = 1, effect = 1, effect.d = 0, equalsize = 0)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_500 <- analysisf(result500)
# 
# 
# result1000 <- simfun(B = 5000, n = 1000, censor = 1, effect = 1, effect.d = 0, equalsize = 0)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_1000 <- analysisf(result1000)
# 
# 
# 
# result5000 <- simfun(B = 5000, n = 5000, censor = 1, effect = 1, effect.d = 0, equalsize = 0)
# # # Second setting, effect = 0, don't want power
# uncensor_balance_5000 <- analysisf(result5000)
# 
# power_censor_unbalance=cbind(uncensor_balance_500,uncensor_balance_1000,uncensor_balance_5000)
# save(power_censor_unbalance,file="power_censor_unbalance.RData")

# load("power_censor_unbalance.RData")
# xtable(power_censor_unbalance)
# % latex table generated in R 4.2.1 by xtable 1.8-4 package
# % Mon Jun 26 11:19:19 2023
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrrr}
# \hline
# & a.param & a.paramsd & a.power & a.param & a.paramsd & a.power & a.param & a.paramsd & a.power \\ 
# \hline
# lin & 15.62 & 0.19 & 0.21 & 15.88 & 0.13 & 0.36 & 15.70 & 0.06 & 0.95 \\ 
# cox & -0.47 & 0.00 & 1.00 & -0.47 & 0.00 & 1.00 & -0.47 & 0.00 & 1.00 \\ 
# med & 36.49 & 0.20 & 0.49 & 36.81 & 0.13 & 0.97 & 36.41 & 0.06 & 1.00 \\ 
# nb & 0.07 & 0.00 & 0.21 & 0.07 & 0.00 & 0.37 & 0.07 & 0.00 & 0.96 \\ 
# poi & 0.07 & 0.00 & 0.92 & 0.07 & 0.00 & 0.96 & 0.07 & 0.00 & 1.00 \\ 
# temp & -15.60 & 0.19 & 0.21 & -15.85 & 0.13 & 0.37 & -15.67 & 0.06 & 0.96 \\ 
# \hline
# \end{tabular}
# \end{table}
##################################################################################
###################################################################################
