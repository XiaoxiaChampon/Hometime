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
# Purpose: Get Histogram for Hometime
# Author:  Xiaoxia Champon, special thanks to Laine Thomas, Sean O'Brien
# Date: 05/28/2024
##############################################################
source("./Source/R/diff_censor_hometime.R")

#get sample data
#B, n, censor, effect,  equalsize,diff_censor,censorbig

n1000_no_censor_no_effect=generate_home_time_scenario(B=2, n=1000, censor=0, effect=0,
                                                        equalsize=0, diff_censor=0,censorbig=1)

n1000_same_censor_no_effect=generate_home_time_scenario(B=2, n=1000, censor=1, effect=0,
                                                     equalsize=0, diff_censor=0,censorbig=1)

n1000_diff_censor_no_effect=generate_home_time_scenario(B=2, n=1000, censor=1, effect=0,
                                                     equalsize=0, diff_censor=1,censorbig=1)
n1000_diff_censor_small_no_effect=generate_home_time_scenario(B=2, n=1000, censor=1, effect=0,
                                                        equalsize=0, diff_censor=1,censorbig=0)

par(mfrow=c(2,2))
hist(n1000_no_censor_no_effect$home_time_data[,,1][,3],main="1000 sub, no censor",xlab="Hometime")
hist(n1000_same_censor_no_effect$home_time_data[,,1][,3],main="1000 sub, same censor",xlab="Hometime")
hist(n1000_diff_censor_no_effect$home_time_data[,,1][,3],main="1000 sub, diff censor,  big",xlab="Hometime")
hist(n1000_diff_censor_small_no_effect$home_time_data[,,1][,3],main="1000 sub, diff censor small",xlab="Hometime")

n1000_no_censor_effect=generate_home_time_scenario(B=2, n=1000, censor=0, effect=1,
                                                      equalsize=0, diff_censor=0,censorbig=1)
#
n1000_same_censor_effect=generate_home_time_scenario(B=2, n=1000, censor=1, effect=1,
                                                        equalsize=0, diff_censor=0,censorbig=1)

n1000_diff_censor_effect=generate_home_time_scenario(B=2, n=1000, censor=1, effect=1,
                                                        equalsize=0, diff_censor=1,censorbig=1)
n1000_diff_censor_small_effect=generate_home_time_scenario(B=2, n=1000, censor=1, effect=1,
                                                              equalsize=0, diff_censor=1,censorbig=0)
par(mfrow=c(2,2))
hist(n1000_no_censor_effect$home_time_data[,,2][,3],main="1000 sub, no censor, effect",xlab="Hometime")
hist(n1000_same_censor_effect$home_time_data[,,2][,3],main="1000 sub, same censor, effect",xlab="Hometime")
hist(n1000_diff_censor_effect$home_time_data[,,2][,3],main="1000 sub, diff censor,  big, effect",xlab="Hometime")
hist(n1000_diff_censor_small_effect$home_time_data[,,2][,3],main="1000 sub, diff censor small, effect",xlab="Hometime")


