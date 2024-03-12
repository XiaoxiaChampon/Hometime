
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




#Function to test 6 different models
#Input: 2D array, n*5, n: the number of the subjects, 
#                      columns:outcome group outcome.t  outcome.b event
#Output: list including two items: parameter (6 elements, one for each model) and 
#                                  power (6 elements, one for each model)
home_time_regressions = function(home_time_data) {
    
    
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
    
    
    m2.2 <- coxph(Surv(time=outcome.t, event=htevent) ~ group,data = home_time_data)
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
    param.est <- c(co1,  co3, co5,co4, co2new, co6)
    names(param.est)=c("lr","median","poission","nb","cox","tp")
    power<- c(p1,  p3, p5, p4,p2new, p6)
    names(power)=c("lr","median","poission","nb","cox","tp")
    return(list("param.est "=param.est ,"power"=power))
}
