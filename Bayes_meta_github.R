#################
# Bayesian meta-analysis with misclassification adjustment
# Citation: Burstyn I, Gustafson P, Boon D, Goodman J, Goldstein ND. Bayesian meta-analysis of effect estimates for ever-use of talc on ovarian cancer risk with adjustment for exposure misclassification. Manuscript in preparation.
# 5/30/26 -- Neal Goldstein (code based on script from Igor Burstyn)
#################


### FUNCTIONS ###

library("rje")
library("rjags") #install JAGS 4.3.2 from here: https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Mac%20OS%20X/
library("MCMCvis")
#install.packages("forestplot")
library("forestplot")
#install.packages("epiR")
library("epiR")


### SETUP ###

#scenario = 1 #hard code for scenario 1...8

scenario_title = paste("Priors of Scenario ", 1:8, sep="") #labels for output
output_dir <- "./talc/" #for publication files

burnin.iter = 50000 #burn in to discard
samples.iter = 500000 #coda samples

set.seed(777) #reproducibility

#automate analysis and output for all scenarios
for (scenario in 1:8) {
  
  cat("\n\n************** ","Scenario: ",scenario," of 8 **************\n",sep="")
  
  
  ### PRIOR DERIVATION FOLLOWING WU ET AL. https://pubmed.ncbi.nlm.nih.gov/41872762/ ###
  
  #code is from supplement 1
  
  #Format for each row = expert best guess, 2.5%ile, 97.5%ile
  #note that case control is differential while cohort is non-differential (but code written to allow differential if priors available)
  #from table A2.1
  
  #expert 1
  SN_case1=c(0.96, 0.94, 0.98)
  SP_case1=c(0.85, 0.80, 0.90) 
  SN_control1=c(.9, .85, .95) 
  SP_control1=c(.9, .85, .95)
  
  SN_cohort1=c(0.9, 0.85, 0.95)
  SP_cohort1=c(0.9, 0.85, 0.95)
  
  #expert 2
  SN_case2=c(0.75, 0.65, 0.85)
  SP_case2=c(0.90, 0.85, 0.95) 
  SN_control2=c(.65, .6, .7) 
  SP_control2=c(.9, .85, .95)
  
  SN_cohort2=c(0.8, 0.75, 0.85)
  SP_cohort2=c(0.9, 0.85, 0.95)
  
  #expert 3
  SN_case3=c(0.90, 0.85, 0.95)
  SP_case3=c(0.85, 0.8, 0.9) 
  SN_control3=c(.8, .75, .85) 
  SP_control3=c(.9, .85, .95)
  
  SN_cohort3=c(0.8, 0.75, 0.85)
  SP_cohort3=c(0.8, 0.78, 0.82)
  
  #combine
  #average best, lowest and highest of expert guesses  for sensitivity in cases
  SN_case=rbind(SN_case1, SN_case2, SN_case3)
  SN_case_mean=mean(SN_case[,1])
  SN_case_low=mean(SN_case[,2])
  SN_case_high=mean(SN_case[,3])
  c(SN_case_mean, SN_case_low, SN_case_high)
  
  SN_cohort=rbind(SN_cohort1, SN_cohort2, SN_cohort3)
  SN_cohort_mean=mean(SN_cohort[,1])
  SN_cohort_low=mean(SN_cohort[,2])
  SN_cohort_high=mean(SN_cohort[,3])
  c(SN_cohort_mean, SN_cohort_low, SN_cohort_high)
  
  #average best, lowest and highest of expert guesses for specificity in cases
  SP_case=rbind(SP_case1, SP_case2, SP_case3)
  SP_case_mean=mean(SP_case[,1])
  SP_case_low=mean(SP_case[,2])
  SP_case_high=mean(SP_case[,3])
  c(SP_case_mean, SP_case_low, SP_case_high)
  
  SP_cohort=rbind(SP_cohort1, SP_cohort2, SP_cohort3)
  SP_cohort_mean=mean(SP_cohort[,1])
  SP_cohort_low=mean(SP_cohort[,2])
  SP_cohort_high=mean(SP_cohort[,3])
  c(SP_cohort_mean, SP_cohort_low, SP_cohort_high)
  
  #average best, lowest and highest of expert guesses  for sensitivity in controls
  SN_control=rbind(SN_control1, SN_control2, SN_control3)
  SN_control_mean=mean(SN_control[,1])
  SN_control_low=mean(SN_control[,2])
  SN_control_high=mean(SN_control[,3])
  c(SN_control_mean, SN_control_low, SN_control_high)
  
  SN_cohort=rbind(SN_cohort1, SN_cohort2, SN_cohort3)
  SN_cohort_mean=mean(SN_cohort[,1])
  SN_cohort_low=mean(SN_cohort[,2])
  SN_cohort_high=mean(SN_cohort[,3])
  c(SN_cohort_mean, SN_cohort_low, SN_cohort_high)
  
  #average best, lowest and highest of expert guesses  for specificity in controls
  SP_control=rbind(SP_control1, SP_control2, SP_control3)
  SP_control_mean=mean(SP_control[,1])
  SP_control_low=mean(SP_control[,2])
  SP_control_high=mean(SP_control[,3])
  c(SP_control_mean, SP_control_low, SP_control_high)
  
  SP_cohort=rbind(SP_cohort1, SP_cohort2, SP_cohort3)
  SP_cohort_mean=mean(SP_cohort[,1])
  SP_cohort_low=mean(SP_cohort[,2])
  SP_cohort_high=mean(SP_cohort[,3])
  c(SP_cohort_mean, SP_cohort_low, SP_cohort_high)
  
  #select scenario for calculating beta shape parameters
  if (scenario==1) {
    
    #average of best guesses (default)
    
    #get parameters of beta distribution of SN cases
    SN1_cc=epi.betabuster(mode=SN_case_mean, conf=0.975, imsure="greater than", x=SN_case_low, conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SN1_cohort=epi.betabuster(mode=SN_cohort_mean, conf=0.975, imsure="greater than", x=SN_cohort_low, conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
    #get parameters of beta distribution of SP cases
    SP1_cc=epi.betabuster(mode=SP_case_mean, conf=0.975, imsure="greater than", x=SP_case_low, conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SP1_cohort=epi.betabuster(mode=SP_cohort_mean, conf=0.975, imsure="greater than", x=SP_cohort_low, conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
    #get parameters of beta distribution of SN controls
    SN0_cc=epi.betabuster(mode=SN_control_mean, conf=0.975, imsure="greater than", x=SN_control_low, conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SN0_cohort=epi.betabuster(mode=SN_cohort_mean, conf=0.975, imsure="greater than", x=SN_cohort_low, conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
    #get parameters of beta distribution of SP controls
    SP0_cc=epi.betabuster(mode=SP_control_mean, conf=0.975, imsure="greater than", x=SP_control_low, conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SP0_cohort=epi.betabuster(mode=SP_cohort_mean, conf=0.975, imsure="greater than", x=SP_cohort_low, conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
  } else if (scenario==2) {
    
    #prior sample size = 10
    prior_SS = 10
    
    SN0_cc = list("shape1"=NA,"shape2"=NA)
    SP0_cc = list("shape1"=NA,"shape2"=NA)
    SN1_cc = list("shape1"=NA,"shape2"=NA)
    SP1_cc = list("shape1"=NA,"shape2"=NA)
    
    SN0_cohort = list("shape1"=NA,"shape2"=NA)
    SP0_cohort = list("shape1"=NA,"shape2"=NA)
    SN1_cohort = list("shape1"=NA,"shape2"=NA)
    SP1_cohort = list("shape1"=NA,"shape2"=NA)
    
    SN0_cc$shape1 = SN_control_mean*prior_SS 
    SN0_cc$shape2 = prior_SS - SN0_cc$shape1
    
    SP0_cc$shape1 = SP_control_mean*prior_SS
    SP0_cc$shape2 = prior_SS - SP0_cc$shape1
    
    SN1_cc$shape1 = SN_case_mean*prior_SS 
    SN1_cc$shape2 = prior_SS - SN1_cc$shape1
    
    SP1_cc$shape1 = SP_case_mean*prior_SS 
    SP1_cc$shape2 = prior_SS - SP1_cc$shape1
    
    SN0_cohort$shape1 = SN_cohort_mean*prior_SS 
    SN0_cohort$shape2 = prior_SS - SN0_cohort$shape1
    
    SP0_cohort$shape1 = SP_cohort_mean*prior_SS
    SP0_cohort$shape2 = prior_SS - SP0_cohort$shape1
    
    SN1_cohort$shape1 = SN_cohort_mean*prior_SS 
    SN1_cohort$shape2 = prior_SS - SN1_cohort$shape1
    
    SP1_cohort$shape1 = SP_cohort_mean*prior_SS 
    SP1_cohort$shape2 = prior_SS - SP1_cohort$shape1
    
    rm(prior_SS)
    
  } else if (scenario==3) {
    
    #prior sample size = 40
    prior_SS = 40
    
    SN0_cc = list("shape1"=NA,"shape2"=NA)
    SP0_cc = list("shape1"=NA,"shape2"=NA)
    SN1_cc = list("shape1"=NA,"shape2"=NA)
    SP1_cc = list("shape1"=NA,"shape2"=NA)
    
    SN0_cohort = list("shape1"=NA,"shape2"=NA)
    SP0_cohort = list("shape1"=NA,"shape2"=NA)
    SN1_cohort = list("shape1"=NA,"shape2"=NA)
    SP1_cohort = list("shape1"=NA,"shape2"=NA)
    
    SN0_cc$shape1 = SN_control_mean*prior_SS 
    SN0_cc$shape2 = prior_SS - SN0_cc$shape1
    
    SP0_cc$shape1 = SP_control_mean*prior_SS
    SP0_cc$shape2 = prior_SS - SP0_cc$shape1
    
    SN1_cc$shape1 = SN_case_mean*prior_SS 
    SN1_cc$shape2 = prior_SS - SN1_cc$shape1
    
    SP1_cc$shape1 = SP_case_mean*prior_SS 
    SP1_cc$shape2 = prior_SS - SP1_cc$shape1
    
    SN0_cohort$shape1 = SN_cohort_mean*prior_SS 
    SN0_cohort$shape2 = prior_SS - SN0_cohort$shape1
    
    SP0_cohort$shape1 = SP_cohort_mean*prior_SS
    SP0_cohort$shape2 = prior_SS - SP0_cohort$shape1
    
    SN1_cohort$shape1 = SN_cohort_mean*prior_SS 
    SN1_cohort$shape2 = prior_SS - SN1_cohort$shape1
    
    SP1_cohort$shape1 = SP_cohort_mean*prior_SS 
    SP1_cohort$shape2 = prior_SS - SP1_cohort$shape1
    
    rm(prior_SS)
    
  } else if (scenario==4) {
    
    #prior sample size = 400
    prior_SS = 400
    
    SN0_cc = list("shape1"=NA,"shape2"=NA)
    SP0_cc = list("shape1"=NA,"shape2"=NA)
    SN1_cc = list("shape1"=NA,"shape2"=NA)
    SP1_cc = list("shape1"=NA,"shape2"=NA)
    
    SN0_cohort = list("shape1"=NA,"shape2"=NA)
    SP0_cohort = list("shape1"=NA,"shape2"=NA)
    SN1_cohort = list("shape1"=NA,"shape2"=NA)
    SP1_cohort = list("shape1"=NA,"shape2"=NA)
    
    SN0_cc$shape1 = SN_control_mean*prior_SS 
    SN0_cc$shape2 = prior_SS - SN0_cc$shape1
    
    SP0_cc$shape1 = SP_control_mean*prior_SS
    SP0_cc$shape2 = prior_SS - SP0_cc$shape1
    
    SN1_cc$shape1 = SN_case_mean*prior_SS 
    SN1_cc$shape2 = prior_SS - SN1_cc$shape1
    
    SP1_cc$shape1 = SP_case_mean*prior_SS 
    SP1_cc$shape2 = prior_SS - SP1_cc$shape1
    
    SN0_cohort$shape1 = SN_cohort_mean*prior_SS 
    SN0_cohort$shape2 = prior_SS - SN0_cohort$shape1
    
    SP0_cohort$shape1 = SP_cohort_mean*prior_SS
    SP0_cohort$shape2 = prior_SS - SP0_cohort$shape1
    
    SN1_cohort$shape1 = SN_cohort_mean*prior_SS 
    SN1_cohort$shape2 = prior_SS - SN1_cohort$shape1
    
    SP1_cohort$shape1 = SP_cohort_mean*prior_SS 
    SP1_cohort$shape2 = prior_SS - SP1_cohort$shape1
    
    rm(prior_SS)
    
  } else if (scenario==5) {
    
    #expert 1 only
    expert=1
    
    #get parameters of beta distribution of SN cases
    SN1_cc=epi.betabuster(mode=SN_case[expert,1], conf=0.975, imsure="greater than", x=SN_case[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SN1_cohort=epi.betabuster(mode=SN_cohort[expert,1], conf=0.975, imsure="greater than", x=SN_cohort[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
    #get parameters of beta distribution of SP cases
    SP1_cc=epi.betabuster(mode=SP_case[expert,1], conf=0.975, imsure="greater than", x=SP_case[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SP1_cohort=epi.betabuster(mode=SP_cohort[expert,1], conf=0.975, imsure="greater than", x=SP_cohort[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
    #get parameters of beta distribution of SN controls
    SN0_cc=epi.betabuster(mode=SN_control[expert,1], conf=0.975, imsure="greater than", x=SN_control[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SN0_cohort=epi.betabuster(mode=SN_cohort[expert,1], conf=0.975, imsure="greater than", x=SN_cohort[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
    #get parameters of beta distribution of SP controls
    SP0_cc=epi.betabuster(mode=SP_control[expert,1], conf=0.975, imsure="greater than", x=SP_control[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SP0_cohort=epi.betabuster(mode=SP_cohort[expert,1], conf=0.975, imsure="greater than", x=SP_cohort[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
    rm(expert)
    
  } else if (scenario==6) {
    
    #expert 2 only
    expert=2
    
    #get parameters of beta distribution of SN cases
    SN1_cc=epi.betabuster(mode=SN_case[expert,1], conf=0.975, imsure="greater than", x=SN_case[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SN1_cohort=epi.betabuster(mode=SN_cohort[expert,1], conf=0.975, imsure="greater than", x=SN_cohort[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
    #get parameters of beta distribution of SP cases
    SP1_cc=epi.betabuster(mode=SP_case[expert,1], conf=0.975, imsure="greater than", x=SP_case[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SP1_cohort=epi.betabuster(mode=SP_cohort[expert,1], conf=0.975, imsure="greater than", x=SP_cohort[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
    #get parameters of beta distribution of SN controls
    SN0_cc=epi.betabuster(mode=SN_control[expert,1], conf=0.975, imsure="greater than", x=SN_control[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SN0_cohort=epi.betabuster(mode=SN_cohort[expert,1], conf=0.975, imsure="greater than", x=SN_cohort[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
    #get parameters of beta distribution of SP controls
    SP0_cc=epi.betabuster(mode=SP_control[expert,1], conf=0.975, imsure="greater than", x=SP_control[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SP0_cohort=epi.betabuster(mode=SP_cohort[expert,1], conf=0.975, imsure="greater than", x=SP_cohort[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
    rm(expert)
    
  } else if (scenario==7) {
    
    #expert 3 only
    expert=3
    
    #get parameters of beta distribution of SN cases
    SN1_cc=epi.betabuster(mode=SN_case[expert,1], conf=0.975, imsure="greater than", x=SN_case[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SN1_cohort=epi.betabuster(mode=SN_cohort[expert,1], conf=0.975, imsure="greater than", x=SN_cohort[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
    #get parameters of beta distribution of SP cases
    SP1_cc=epi.betabuster(mode=SP_case[expert,1], conf=0.975, imsure="greater than", x=SP_case[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SP1_cohort=epi.betabuster(mode=SP_cohort[expert,1], conf=0.975, imsure="greater than", x=SP_cohort[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
    #get parameters of beta distribution of SN controls
    SN0_cc=epi.betabuster(mode=SN_control[expert,1], conf=0.975, imsure="greater than", x=SN_control[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SN0_cohort=epi.betabuster(mode=SN_cohort[expert,1], conf=0.975, imsure="greater than", x=SN_cohort[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
    #get parameters of beta distribution of SP controls
    SP0_cc=epi.betabuster(mode=SP_control[expert,1], conf=0.975, imsure="greater than", x=SP_control[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SP0_cohort=epi.betabuster(mode=SP_cohort[expert,1], conf=0.975, imsure="greater than", x=SP_cohort[expert,2], conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
    rm(expert)
    
  } else if (scenario==8) {
    
    #mean of highest value for 95% prior density
    
    #get parameters of beta distribution of SN cases
    SN1_cc=epi.betabuster(mode=SN_case_mean, conf=0.975, imsure="less than", x=SN_case_high, conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SN1_cohort=epi.betabuster(mode=SN_cohort_mean, conf=0.975, imsure="less than", x=SN_cohort_high, conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
    #get parameters of beta distribution of SP cases
    SP1_cc=epi.betabuster(mode=SP_case_mean, conf=0.975, imsure="less than", x=SP_case_high, conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SP1_cohort=epi.betabuster(mode=SP_cohort_mean, conf=0.975, imsure="less than", x=SP_cohort_high, conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
    #get parameters of beta distribution of SN controls
    SN0_cc=epi.betabuster(mode=SN_control_mean, conf=0.975, imsure="less than", x=SN_control_high, conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SN0_cohort=epi.betabuster(mode=SN_cohort_mean, conf=0.975, imsure="less than", x=SN_cohort_high, conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
    #get parameters of beta distribution of SP controls
    SP0_cc=epi.betabuster(mode=SP_control_mean, conf=0.975, imsure="less than", x=SP_control_high, conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    SP0_cohort=epi.betabuster(mode=SP_cohort_mean, conf=0.975, imsure="less than", x=SP_cohort_high, conf.level = 0.95, max.shape1 = 1000, step = 0.001)
    
  } else {
    
    stop("unknown scenario")
    
  }
  
  SN1_cc$shape1
  SN1_cc$shape2
  SN1_cohort$shape1
  SN1_cohort$shape2
  
  SP1_cc$shape1
  SP1_cc$shape2
  SP1_cohort$shape1
  SP1_cohort$shape2
  
  SN0_cc$shape1
  SN0_cc$shape2
  SN0_cohort$shape1
  SN0_cohort$shape2
  
  SP0_cc$shape1
  SP0_cc$shape2
  SP0_cohort$shape1
  SP0_cohort$shape2
  
  
  ### READ DATA ###
  
  #notes for data setup 
  
  #study specific 2x2 tables from IARC Table A2.2 (append .m for each study 1...m)
  #x0.m = exposed controls
  #x1.m = exposed cases
  #n0.m = total controls
  #n1.m = total cases
  
  #classifier accuracy, we vary this per study (append .m for each study 1...m)
  #parameters of priors on SN and SP as Beta(shape1=a,shape2=b)
  #allows for differential exposure misclassification among referents (SN0/SP0) and cases (SN1/SP1)
  #a.sn0.m = sensitivity (shape1), controls
  #b.sn0.m = sensitivity (shape2), controls
  #a.sp0.m = specificity (shape1), controls
  #b.sp0.m = specificity (shape2), controls
  #a.sn1.m = sensitivity (shape1), cases
  #b.sn1.m = sensitivity (shape2), cases
  #a.sp1.m = specificity (shape1), cases
  #b.sp1.m = specificity (shape2), cases
  
  #Wu et al priors: this is per study
  a.sn0.1=SN0_cohort$shape1; b.sn0.1=SN0_cohort$shape2
  a.sp0.1=SP0_cohort$shape1; b.sp0.1=SP0_cohort$shape2
  a.sn1.1=SN1_cohort$shape1; b.sn1.1=SN1_cohort$shape2
  a.sp1.1=SP1_cohort$shape1; b.sp1.1=SP1_cohort$shape2
  
  a.sn0.2=SN0_cohort$shape1; b.sn0.2=SN0_cohort$shape2
  a.sp0.2=SP0_cohort$shape1; b.sp0.2=SP0_cohort$shape2
  a.sn1.2=SN1_cohort$shape1; b.sn1.2=SN1_cohort$shape2
  a.sp1.2=SP1_cohort$shape1; b.sp1.2=SP1_cohort$shape2
  
  a.sn0.3=SN0_cohort$shape1; b.sn0.3=SN0_cohort$shape2
  a.sp0.3=SP0_cohort$shape1; b.sp0.3=SP0_cohort$shape2
  a.sn1.3=SN1_cohort$shape1; b.sn1.3=SN1_cohort$shape2
  a.sp1.3=SP1_cohort$shape1; b.sp1.3=SP1_cohort$shape2
  
  a.sn0.4=SN0_cohort$shape1; b.sn0.4=SN0_cohort$shape2
  a.sp0.4=SP0_cohort$shape1; b.sp0.4=SP0_cohort$shape2
  a.sn1.4=SN1_cohort$shape1; b.sn1.4=SN1_cohort$shape2
  a.sp1.4=SP1_cohort$shape1; b.sp1.4=SP1_cohort$shape2
  
  a.sn0.5=SN0_cc$shape1; b.sn0.5=SN0_cc$shape2
  a.sp0.5=SP0_cc$shape1; b.sp0.5=SP0_cc$shape2
  a.sn1.5=SN1_cc$shape1; b.sn1.5=SN1_cc$shape2
  a.sp1.5=SP1_cc$shape1; b.sp1.5=SP1_cc$shape2
  
  a.sn0.6=SN0_cc$shape1; b.sn0.6=SN0_cc$shape2
  a.sp0.6=SP0_cc$shape1; b.sp0.6=SP0_cc$shape2
  a.sn1.6=SN1_cc$shape1; b.sn1.6=SN1_cc$shape2
  a.sp1.6=SP1_cc$shape1; b.sp1.6=SP1_cc$shape2
  
  a.sn0.7=SN0_cc$shape1; b.sn0.7=SN0_cc$shape2
  a.sp0.7=SP0_cc$shape1; b.sp0.7=SP0_cc$shape2
  a.sn1.7=SN1_cc$shape1; b.sn1.7=SN1_cc$shape2
  a.sp1.7=SP1_cc$shape1; b.sp1.7=SP1_cc$shape2
  
  a.sn0.8=SN0_cc$shape1; b.sn0.8=SN0_cc$shape2
  a.sp0.8=SP0_cc$shape1; b.sp0.8=SP0_cc$shape2
  a.sn1.8=SN1_cc$shape1; b.sn1.8=SN1_cc$shape2
  a.sp1.8=SP1_cc$shape1; b.sp1.8=SP1_cc$shape2
  
  a.sn0.9=SN0_cc$shape1; b.sn0.9=SN0_cc$shape2
  a.sp0.9=SP0_cc$shape1; b.sp0.9=SP0_cc$shape2
  a.sn1.9=SN1_cc$shape1; b.sn1.9=SN1_cc$shape2
  a.sp1.9=SP1_cc$shape1; b.sp1.9=SP1_cc$shape2
  
  a.sn0.10=SN0_cc$shape1; b.sn0.10=SN0_cc$shape2
  a.sp0.10=SP0_cc$shape1; b.sp0.10=SP0_cc$shape2
  a.sn1.10=SN1_cc$shape1; b.sn1.10=SN1_cc$shape2
  a.sp1.10=SP1_cc$shape1; b.sp1.10=SP1_cc$shape2
  
  a.sn0.11=SN0_cc$shape1; b.sn0.11=SN0_cc$shape2
  a.sp0.11=SP0_cc$shape1; b.sp0.11=SP0_cc$shape2
  a.sn1.11=SN1_cc$shape1; b.sn1.11=SN1_cc$shape2
  a.sp1.11=SP1_cc$shape1; b.sp1.11=SP1_cc$shape2
  
  a.sn0.12=SN0_cc$shape1; b.sn0.12=SN0_cc$shape2
  a.sp0.12=SP0_cc$shape1; b.sp0.12=SP0_cc$shape2
  a.sn1.12=SN1_cc$shape1; b.sn1.12=SN1_cc$shape2
  a.sp1.12=SP1_cc$shape1; b.sp1.12=SP1_cc$shape2
  
  a.sn0.13=SN0_cc$shape1; b.sn0.13=SN0_cc$shape2
  a.sp0.13=SP0_cc$shape1; b.sp0.13=SP0_cc$shape2
  a.sn1.13=SN1_cc$shape1; b.sn1.13=SN1_cc$shape2
  a.sp1.13=SP1_cc$shape1; b.sp1.13=SP1_cc$shape2
  
  a.sn0.14=SN0_cc$shape1; b.sn0.14=SN0_cc$shape2
  a.sp0.14=SP0_cc$shape1; b.sp0.14=SP0_cc$shape2
  a.sn1.14=SN1_cc$shape1; b.sn1.14=SN1_cc$shape2
  a.sp1.14=SP1_cc$shape1; b.sp1.14=SP1_cc$shape2
  
  a.sn0.15=SN0_cc$shape1; b.sn0.15=SN0_cc$shape2
  a.sp0.15=SP0_cc$shape1; b.sp0.15=SP0_cc$shape2
  a.sn1.15=SN1_cc$shape1; b.sn1.15=SN1_cc$shape2
  a.sp1.15=SP1_cc$shape1; b.sp1.15=SP1_cc$shape2
  
  #and uniform prior on exposure probability (r) Beta(1,1)
  aa=1
  bb=1
  
  #study 1...15 (1...4 are cohort, 5...15 case-control)
  #IARCin Table A2.2 gives fractional counts of cases and controls, not as integers.
  ##This is not allowed when specifying binomial distribution of events for given probability and sample size
  ##Therefore, we rounded up the counts to the nearest integer
  IARC_data = list(
    x0.1=32413, x1.1=514, n0.1=79055, n1.1=1224,
    x0.2=15721, x1.2=18, n0.2=60464, n1.2=76,
    x0.3=10852, x1.3=64, n0.3=40193, n1.3=219,
    x0.4=37558, x1.4=363, n0.4=70865, n1.4=649,
    x0.5=658, x1.5=705, n0.5=963, n1.5=1005,
    x0.6=297, x1.6=272, n0.6=1841, n1.6=1565,
    x0.7=112, x1.7=74, n0.7=601, n1.7=400,
    x0.8=316, x1.8=194, n0.8=1305, n1.8=633,
    x0.9=122, x1.9=195, n0.9=513, n1.9=664,
    x0.10=636, x1.10=755, n0.10=1875, n1.10=1884,
    x0.11=200, x1.11=197, n0.11=564, n1.11=449,
    x0.12=170, x1.12=208, n0.12=664, n1.12=643,
    x0.13=202, x1.13=119, n0.13=596, n1.13=315,
    x0.14=15, x1.14=14, n0.14=80, n1.14=44,
    x0.15=75, x1.15=53, n0.15=421, n1.15=233,
    a.sn0.1=a.sn0.1, b.sn0.1=b.sn0.1,
    a.sp0.1=a.sp0.1, b.sp0.1=b.sp0.1,
    a.sn1.1=a.sn1.1, b.sn1.1=b.sn1.1,
    a.sp1.1=a.sp1.1, b.sp1.1=b.sp1.1,
    a.sn0.2=a.sn0.2, b.sn0.2=b.sn0.2,
    a.sp0.2=a.sp0.2, b.sp0.2=b.sp0.2,
    a.sn1.2=a.sn1.2, b.sn1.2=b.sn1.2,
    a.sp1.2=a.sp1.2, b.sp1.2=b.sp1.2,
    a.sn0.3=a.sn0.3, b.sn0.3=b.sn0.3,
    a.sp0.3=a.sp0.3, b.sp0.3=b.sp0.3,
    a.sn1.3=a.sn1.3, b.sn1.3=b.sn1.3,
    a.sp1.3=a.sp1.3, b.sp1.3=b.sp1.3,
    a.sn0.4=a.sn0.4, b.sn0.4=b.sn0.4,
    a.sp0.4=a.sp0.4, b.sp0.4=b.sp0.4,
    a.sn1.4=a.sn1.4, b.sn1.4=b.sn1.4,
    a.sp1.4=a.sp1.4, b.sp1.4=b.sp1.4,
    a.sn0.5=a.sn0.5, b.sn0.5=b.sn0.5,
    a.sp0.5=a.sp0.5, b.sp0.5=b.sp0.5,
    a.sn1.5=a.sn1.5, b.sn1.5=b.sn1.5,
    a.sp1.5=a.sp1.5, b.sp1.5=b.sp1.5,
    a.sn0.6=a.sn0.6, b.sn0.6=b.sn0.6,
    a.sp0.6=a.sp0.6, b.sp0.6=b.sp0.6,
    a.sn1.6=a.sn1.6, b.sn1.6=b.sn1.6,
    a.sp1.6=a.sp1.6, b.sp1.6=b.sp1.6,
    a.sn0.7=a.sn0.7, b.sn0.7=b.sn0.7,
    a.sp0.7=a.sp0.7, b.sp0.7=b.sp0.7,
    a.sn1.7=a.sn1.7, b.sn1.7=b.sn1.7,
    a.sp1.7=a.sp1.7, b.sp1.7=b.sp1.7,
    a.sn0.8=a.sn0.8, b.sn0.8=b.sn0.8,
    a.sp0.8=a.sp0.8, b.sp0.8=b.sp0.8,
    a.sn1.8=a.sn1.8, b.sn1.8=b.sn1.8,
    a.sp1.8=a.sp1.8, b.sp1.8=b.sp1.8,
    a.sn0.9=a.sn0.9, b.sn0.9=b.sn0.9,
    a.sp0.9=a.sp0.9, b.sp0.9=b.sp0.9,
    a.sn1.9=a.sn1.9, b.sn1.9=b.sn1.9,
    a.sp1.9=a.sp1.9, b.sp1.9=b.sp1.9,
    a.sn0.10=a.sn0.10, b.sn0.10=b.sn0.10,
    a.sp0.10=a.sp0.10, b.sp0.10=b.sp0.10,
    a.sn1.10=a.sn1.10, b.sn1.10=b.sn1.10,
    a.sp1.10=a.sp1.10, b.sp1.10=b.sp1.10,
    a.sn0.11=a.sn0.11, b.sn0.11=b.sn0.11,
    a.sp0.11=a.sp0.11, b.sp0.11=b.sp0.11,
    a.sn1.11=a.sn1.11, b.sn1.11=b.sn1.11,
    a.sp1.11=a.sp1.11, b.sp1.11=b.sp1.11,
    a.sn0.12=a.sn0.12, b.sn0.12=b.sn0.12,
    a.sp0.12=a.sp0.12, b.sp0.12=b.sp0.12,
    a.sn1.12=a.sn1.12, b.sn1.12=b.sn1.12,
    a.sp1.12=a.sp1.12, b.sp1.12=b.sp1.12,
    a.sn0.13=a.sn0.13, b.sn0.13=b.sn0.13,
    a.sp0.13=a.sp0.13, b.sp0.13=b.sp0.13,
    a.sn1.13=a.sn1.13, b.sn1.13=b.sn1.13,
    a.sp1.13=a.sp1.13, b.sp1.13=b.sp1.13,
    a.sn0.14=a.sn0.14, b.sn0.14=b.sn0.14,
    a.sp0.14=a.sp0.14, b.sp0.14=b.sp0.14,
    a.sn1.14=a.sn1.14, b.sn1.14=b.sn1.14,
    a.sp1.14=a.sp1.14, b.sp1.14=b.sp1.14,
    a.sn0.15=a.sn0.15, b.sn0.15=b.sn0.15,
    a.sp0.15=a.sp0.15, b.sp0.15=b.sp0.15,
    a.sn1.15=a.sn1.15, b.sn1.15=b.sn1.15,
    a.sp1.15=a.sp1.15, b.sp1.15=b.sp1.15,
    aa=aa, bb=bb,
    nstudy=15)
  
  
  ### JAGS MODEL ###
  
  IARC_model <- "model{

#DATA

#study 1...15 (1...4 are cohort, 5...15 case-control)
x0.1 ~ dbin(p0.1, n0.1)
x1.1 ~ dbin(p1.1, n1.1)

x0.2 ~ dbin(p0.2, n0.2)
x1.2 ~ dbin(p1.2, n1.2)

x0.3 ~ dbin(p0.3, n0.3)
x1.3 ~ dbin(p1.3, n1.3)

x0.4 ~ dbin(p0.4, n0.4)
x1.4 ~ dbin(p1.4, n1.4)

x0.5 ~ dbin(p0.5, n0.5)
x1.5 ~ dbin(p1.5, n1.5)

x0.6 ~ dbin(p0.6, n0.6)
x1.6 ~ dbin(p1.6, n1.6)

x0.7 ~ dbin(p0.7, n0.7)
x1.7 ~ dbin(p1.7, n1.7)

x0.8 ~ dbin(p0.8, n0.8)
x1.8 ~ dbin(p1.8, n1.8)

x0.9 ~ dbin(p0.9, n0.9)
x1.9 ~ dbin(p1.9, n1.9)

x0.10 ~ dbin(p0.10, n0.10)
x1.10 ~ dbin(p1.10, n1.10)

x0.11 ~ dbin(p0.11, n0.11)
x1.11 ~ dbin(p1.11, n1.11)

x0.12 ~ dbin(p0.12, n0.12)
x1.12 ~ dbin(p1.12, n1.12)

x0.13 ~ dbin(p0.13, n0.13)
x1.13 ~ dbin(p1.13, n1.13)

x0.14 ~ dbin(p0.14, n0.14)
x1.14 ~ dbin(p1.14, n1.14)

x0.15 ~ dbin(p0.15, n0.15)
x1.15 ~ dbin(p1.15, n1.15)

#MODELS

#study 1...15 (1...4 are cohort, 5...15 case-control)
p0.1 <- r0.1*SN0.1 + (1-r0.1)*(1-SP0.1)
p1.1 <- r1.1*SN1.1 + (1-r1.1)*(1-SP1.1)
r1.1 <- (OR.1*r0.1)/(1-r0.1+OR.1*r0.1)

p0.2 <- r0.2*SN0.2 + (1-r0.2)*(1-SP0.2)
p1.2 <- r1.2*SN1.2 + (1-r1.2)*(1-SP1.2)
r1.2 <- (OR.2*r0.2)/(1-r0.2+OR.2*r0.2)

p0.3 <- r0.3*SN0.3 + (1-r0.3)*(1-SP0.3)
p1.3 <- r1.3*SN1.3 + (1-r1.3)*(1-SP1.3)
r1.3 <- (OR.3*r0.3)/(1-r0.3+OR.3*r0.3)

p0.4 <- r0.4*SN0.4 + (1-r0.4)*(1-SP0.4)
p1.4 <- r1.4*SN1.4 + (1-r1.4)*(1-SP1.4)
r1.4 <- (OR.4*r0.4)/(1-r0.4+OR.4*r0.4)

p0.5 <- r0.5*SN0.5 + (1-r0.5)*(1-SP0.5)
p1.5 <- r1.5*SN1.5 + (1-r1.5)*(1-SP1.5)
r1.5 <- (OR.5*r0.5)/(1-r0.5+OR.5*r0.5)

p0.6 <- r0.6*SN0.6 + (1-r0.6)*(1-SP0.6)
p1.6 <- r1.6*SN1.6 + (1-r1.6)*(1-SP1.6)
r1.6 <- (OR.6*r0.6)/(1-r0.6+OR.6*r0.6)

p0.7 <- r0.7*SN0.7 + (1-r0.7)*(1-SP0.7)
p1.7 <- r1.7*SN1.7 + (1-r1.7)*(1-SP1.7)
r1.7 <- (OR.7*r0.7)/(1-r0.7+OR.7*r0.7)

p0.8 <- r0.8*SN0.8 + (1-r0.8)*(1-SP0.8)
p1.8 <- r1.8*SN1.8 + (1-r1.8)*(1-SP1.8)
r1.8 <- (OR.8*r0.8)/(1-r0.8+OR.8*r0.8)

p0.9 <- r0.9*SN0.9 + (1-r0.9)*(1-SP0.9)
p1.9 <- r1.9*SN1.9 + (1-r1.9)*(1-SP1.9)
r1.9 <- (OR.9*r0.9)/(1-r0.9+OR.9*r0.9)

p0.10 <- r0.10*SN0.10 + (1-r0.10)*(1-SP0.10)
p1.10 <- r1.10*SN1.10 + (1-r1.10)*(1-SP1.10)
r1.10 <- (OR.10*r0.10)/(1-r0.10+OR.10*r0.10)

p0.11 <- r0.11*SN0.11 + (1-r0.11)*(1-SP0.11)
p1.11 <- r1.11*SN1.11 + (1-r1.11)*(1-SP1.11)
r1.11 <- (OR.11*r0.11)/(1-r0.11+OR.11*r0.11)

p0.12 <- r0.12*SN0.12 + (1-r0.12)*(1-SP0.12)
p1.12 <- r1.12*SN1.12 + (1-r1.12)*(1-SP1.12)
r1.12 <- (OR.12*r0.12)/(1-r0.12+OR.12*r0.12)

p0.13 <- r0.13*SN0.13 + (1-r0.13)*(1-SP0.13)
p1.13 <- r1.13*SN1.13 + (1-r1.13)*(1-SP1.13)
r1.13 <- (OR.13*r0.13)/(1-r0.13+OR.13*r0.13)

p0.14 <- r0.14*SN0.14 + (1-r0.14)*(1-SP0.14)
p1.14 <- r1.14*SN1.14 + (1-r1.14)*(1-SP1.14)
r1.14 <- (OR.14*r0.14)/(1-r0.14+OR.14*r0.14)

p0.15 <- r0.15*SN0.15 + (1-r0.15)*(1-SP0.15)
p1.15 <- r1.15*SN1.15 + (1-r1.15)*(1-SP1.15)
r1.15 <- (OR.15*r0.15)/(1-r0.15+OR.15*r0.15)

#to adjust for confounding:
#RR.conf = RR_crude/RR_adjusted
#RR.conf*RR_adjusted = RR_crude
#CAN GET THIS RATIO FROM DATA BY RATIO OF DISTRIBUTIONS IN PRIOR
#LOG(RR_crude)~N(mc,vc) */* log(RR_adjusted)~N(ma,va)
#r1.1 <- (OR.1*RR.conf.1*r0.1)/(1-r0.1+OR.1*r0.1*RR.conf.1)

#PRIORS

#can make informative as per prior elucidation but flat initially as beta(1,1)
r0.1  ~ dbeta(aa,bb) 
r0.2  ~ dbeta(aa,bb)
r0.3  ~ dbeta(aa,bb)
r0.4  ~ dbeta(aa,bb) 
r0.5  ~ dbeta(aa,bb)
r0.6  ~ dbeta(aa,bb)
r0.7  ~ dbeta(aa,bb) 
r0.8  ~ dbeta(aa,bb)
r0.9  ~ dbeta(aa,bb)
r0.10  ~ dbeta(aa,bb) 
r0.11  ~ dbeta(aa,bb)
r0.12  ~ dbeta(aa,bb)
r0.13  ~ dbeta(aa,bb) 
r0.14  ~ dbeta(aa,bb)
r0.15  ~ dbeta(aa,bb)

#hyper-priors
##between-study variance

lor.1 ~ dnorm(mu,lambda)
lor.2 ~ dnorm(mu,lambda) 
lor.3 ~ dnorm(mu,lambda) 
lor.4 ~ dnorm(mu,lambda)
lor.5 ~ dnorm(mu,lambda) 
lor.6 ~ dnorm(mu,lambda) 
lor.7 ~ dnorm(mu,lambda)
lor.8 ~ dnorm(mu,lambda) 
lor.9 ~ dnorm(mu,lambda) 
lor.10 ~ dnorm(mu,lambda)
lor.11 ~ dnorm(mu,lambda) 
lor.12 ~ dnorm(mu,lambda) 
lor.13 ~ dnorm(mu,lambda)
lor.14 ~ dnorm(mu,lambda) 
lor.15 ~ dnorm(mu,lambda) 

#vague prior on precision, following: https://chjackson.github.io/openbugsdoc/Examples/Blockers.html
lambda ~ dgamma(0.1, 0.1) 

##within-study variance, following: https://chjackson.github.io/openbugsdoc/Examples/Blockers.html
#mu ~ dnorm(0, 0.001)
mu ~ dnorm(0,tau) 
tau ~ dgamma(0.1, 0.1) #same prior on precision withing studies as for between studies

#misclassification: priors differ per study 1...m
SN0.1 ~ dbeta(a.sn0.1, b.sn0.1)
SP0.1 ~ dbeta(a.sp0.1, b.sp0.1)
SN1.1 ~ dbeta(a.sn1.1, b.sn1.1)
SP1.1 ~ dbeta(a.sp1.1, b.sp1.1)

SN0.2 ~ dbeta(a.sn0.2, b.sn0.2)
SP0.2 ~ dbeta(a.sp0.2, b.sp0.2)
SN1.2 ~ dbeta(a.sn1.2, b.sn1.2)
SP1.2 ~ dbeta(a.sp1.2, b.sp1.2)

SN0.3 ~ dbeta(a.sn0.3, b.sn0.3)
SP0.3 ~ dbeta(a.sp0.3, b.sp0.3)
SN1.3 ~ dbeta(a.sn1.3, b.sn1.3)
SP1.3 ~ dbeta(a.sp1.3, b.sp1.3)

SN0.4 ~ dbeta(a.sn0.4, b.sn0.4)
SP0.4 ~ dbeta(a.sp0.4, b.sp0.4)
SN1.4 ~ dbeta(a.sn1.4, b.sn1.4)
SP1.4 ~ dbeta(a.sp1.4, b.sp1.4)

SN0.5 ~ dbeta(a.sn0.5, b.sn0.5)
SP0.5 ~ dbeta(a.sp0.5, b.sp0.5)
SN1.5 ~ dbeta(a.sn1.5, b.sn1.5)
SP1.5 ~ dbeta(a.sp1.5, b.sp1.5)

SN0.6 ~ dbeta(a.sn0.6, b.sn0.6)
SP0.6 ~ dbeta(a.sp0.6, b.sp0.6)
SN1.6 ~ dbeta(a.sn1.6, b.sn1.6)
SP1.6 ~ dbeta(a.sp1.6, b.sp1.6)

SN0.7 ~ dbeta(a.sn0.7, b.sn0.7)
SP0.7 ~ dbeta(a.sp0.7, b.sp0.7)
SN1.7 ~ dbeta(a.sn1.7, b.sn1.7)
SP1.7 ~ dbeta(a.sp1.7, b.sp1.7)

SN0.8 ~ dbeta(a.sn0.8, b.sn0.8)
SP0.8 ~ dbeta(a.sp0.8, b.sp0.8)
SN1.8 ~ dbeta(a.sn1.8, b.sn1.8)
SP1.8 ~ dbeta(a.sp1.8, b.sp1.8)

SN0.9 ~ dbeta(a.sn0.9, b.sn0.9)
SP0.9 ~ dbeta(a.sp0.9, b.sp0.9)
SN1.9 ~ dbeta(a.sn1.9, b.sn1.9)
SP1.9 ~ dbeta(a.sp1.9, b.sp1.9)

SN0.10 ~ dbeta(a.sn0.10, b.sn0.10)
SP0.10 ~ dbeta(a.sp0.10, b.sp0.10)
SN1.10 ~ dbeta(a.sn1.10, b.sn1.10)
SP1.10 ~ dbeta(a.sp1.10, b.sp1.10)

SN0.11 ~ dbeta(a.sn0.11, b.sn0.11)
SP0.11 ~ dbeta(a.sp0.11, b.sp0.11)
SN1.11 ~ dbeta(a.sn1.11, b.sn1.11)
SP1.11 ~ dbeta(a.sp1.11, b.sp1.11)

SN0.12 ~ dbeta(a.sn0.12, b.sn0.12)
SP0.12 ~ dbeta(a.sp0.12, b.sp0.12)
SN1.12 ~ dbeta(a.sn1.12, b.sn1.12)
SP1.12 ~ dbeta(a.sp1.12, b.sp1.12)

SN0.13 ~ dbeta(a.sn0.13, b.sn0.13)
SP0.13 ~ dbeta(a.sp0.13, b.sp0.13)
SN1.13 ~ dbeta(a.sn1.13, b.sn1.13)
SP1.13 ~ dbeta(a.sp1.13, b.sp1.13)

SN0.14 ~ dbeta(a.sn0.14, b.sn0.14)
SP0.14 ~ dbeta(a.sp0.14, b.sp0.14)
SN1.14 ~ dbeta(a.sn1.14, b.sn1.14)
SP1.14 ~ dbeta(a.sp1.14, b.sp1.14)

SN0.15 ~ dbeta(a.sn0.15, b.sn0.15)
SP0.15 ~ dbeta(a.sp0.15, b.sp0.15)
SN1.15 ~ dbeta(a.sn1.15, b.sn1.15)
SP1.15 ~ dbeta(a.sp1.15, b.sp1.15)

#calculations of parameters of interest
OR.1 <- exp(lor.1)
OR.2 <- exp(lor.2)
OR.3 <- exp(lor.3)
OR.4 <- exp(lor.4)
OR.5 <- exp(lor.5)
OR.6 <- exp(lor.6)
OR.7 <- exp(lor.7)
OR.8 <- exp(lor.8)
OR.9 <- exp(lor.9)
OR.10 <- exp(lor.10)
OR.11 <- exp(lor.11)
OR.12 <- exp(lor.12)
OR.13 <- exp(lor.13)
OR.14 <- exp(lor.14)
OR.15 <- exp(lor.15)
OR.gm <- exp((lor.1+lor.2+lor.3+lor.4+lor.5+lor.6+lor.7+lor.8+lor.9+lor.10+lor.11+lor.12+lor.13+lor.14+lor.15) / nstudy) #geometric mean
OR.meta <- exp(mu)

#cohort misclassification (studies 1...4)
cohort_misclass = abs(SN0.1-SN1.1) + abs(SP0.1-SP1.1) + abs(SN0.2-SN1.2) + abs(SP0.2-SP1.2) + abs(SN0.3-SN1.3) + abs(SP0.3-SP1.3) + abs(SN0.4-SN1.4) + abs(SP0.4-SP1.4)
}"
  
  
  ### RUN GIBBS SAMPLER ###
  
  mod <- jags.model(textConnection(IARC_model),
                    data=IARC_data,
                    n.chains=3)
  
  update(mod, burnin.iter) #burn-in
  
  meta_vars <- c("OR.1","OR.2","OR.3","OR.4","OR.5","OR.6","OR.7","OR.8","OR.9","OR.10","OR.11","OR.12","OR.13","OR.14","OR.15","OR.meta", 
                 "SN0.1", "SP0.1","SN1.1", "SP1.1","SN0.2", "SP0.2","SN1.2", "SP1.2","SN0.3", "SP0.3","SN1.3", "SP1.3","SN0.4", "SP0.4","SN1.4", "SP1.4","SN0.5", "SP0.5","SN1.5", "SP1.5","SN0.6", "SP0.6","SN1.6", "SP1.6","SN0.7", "SP0.7","SN1.7", "SP1.7","SN0.8", "SP0.8","SN1.8", "SP1.8","SN0.9", "SP0.9","SN1.9", "SP1.9","SN0.10", "SP0.10","SN1.10", "SP1.10","SN0.11", "SP0.11","SN1.11", "SP1.11","SN0.12", "SP0.12","SN1.12", "SP1.12","SN0.13", "SP0.13","SN1.13", "SP1.13","SN0.14", "SP0.14","SN1.14", "SP1.14","SN0.15", "SP0.15","SN1.15", "SP1.15",
                 "cohort_misclass")
  
  opt.JAGS <- coda.samples(mod, n.iter=samples.iter, variable.names=meta_vars)
  
  meta_results <- MCMCsummary(opt.JAGS)
  
  
  ### FOREST PLOT of ORs ###
  
  #to align with IARC Table A2.2
  study_labels = c("NHS-I (O'Brien et al., 2020)", "NHS-II (O'Brien et al., 2020)", "SIS (O'Brien et al., 2020)", "WHI-OS (O'Brien et al., 2020)", "AUS (Terry et al., 2013)", "DOV (Terry et al., 2013)", "HAW (Terry et al., 2013)", "HOP (Terry et al., 2013)", "NCO (Terry et al., 2013)", "NEC (Terry et al., 2013)", "SON (Terry et al., 2013)", "USC (Terry et al., 2013)", "AACES_B (Davis et al., 2021)", "CCCS_B (Davis et al., 2021)", "CCCS_W (Davis et al., 2021)")
  sort_order = c(2,9:16,3:8,17)
  
  #for base R vignette, see: https://web.archive.org/web/20201101045358/https://cran.r-project.org/web/packages/forestplot/vignettes/forestplot.html
  pdf(paste(output_dir,"Figure_forestplot_OR_",scenario_title[scenario],".pdf",sep=""),height=6,width=10,onefile=F) 
  #tiff(paste(output_dir,"Figure_forestplot_OR_",scenario_title[scenario],".tiff",sep=""),height=6,width=10,units='in',res=600) #uncomment if need high resolution TIFF
  fp <- forestplot(
    title=scenario_title[scenario],
    labeltext = c(study_labels,"Summary"),
    mean = meta_results$`50%`[sort_order],
    lower = meta_results$`2.5%`[sort_order],
    upper = meta_results$`97.5%`[sort_order],
    is.summary = c(rep(FALSE,15),T),
    zero = 1,
    xlab = "Odds Ratio",
    clip=c(0.5, 2.5),
    hrzl_lines = gpar(col="#444444"),
    col=fpColors(box="royalblue",line="darkblue", summary="royalblue")
  )
  print(fp) #workaround needed when embedded in a for loop
  dev.off()
  
  
  ### FOREST PLOT of SN/SP ###
  
  #to align with IARC Table A2.2
  study_labels = c("NHS-I (O'Brien et al., 2020)", "NHS-II (O'Brien et al., 2020)", "SIS (O'Brien et al., 2020)", "WHI-OS (O'Brien et al., 2020)", "AUS (Terry et al., 2013)", "DOV (Terry et al., 2013)", "HAW (Terry et al., 2013)", "HOP (Terry et al., 2013)", "NCO (Terry et al., 2013)", "NEC (Terry et al., 2013)", "SON (Terry et al., 2013)", "USC (Terry et al., 2013)", "AACES_B (Davis et al., 2021)", "CCCS_B (Davis et al., 2021)", "CCCS_W (Davis et al., 2021)")
  sort_order_sn0 = c(18,25:32,19:24)
  sort_order_sn1 = c(33,40:47,34:39)
  sort_order_sp0 = c(48,55:62,49:54)
  sort_order_sp1 = c(63,70:77,64:69)
  
  pdf(paste(output_dir,"Figure_forestplot_SNSP_",scenario_title[scenario],".pdf",sep=""),height=6,width=10,onefile=F) 
  #tiff("paste(output_dir,"Figure_forestplot_SNSP_",scenario_title[scenario],".tiff",sep=""),height=6,width=10,units='in',res=600) #uncomment if need high resolution TIFF
  fp <- forestplot(
    title=scenario_title[scenario],
    labeltext = c(study_labels),
    legend = c("SN0","SN1","SP0","SP1"),
    mean = cbind(meta_results$`50%`[sort_order_sn0],meta_results$`50%`[sort_order_sn1],meta_results$`50%`[sort_order_sp0],meta_results$`50%`[sort_order_sp1]),
    lower = cbind(meta_results$`2.5%`[sort_order_sn0],meta_results$`2.5%`[sort_order_sn1],meta_results$`2.5%`[sort_order_sp0],meta_results$`2.5%`[sort_order_sp1]),
    upper = cbind(meta_results$`97.5%`[sort_order_sn0],meta_results$`97.5%`[sort_order_sn1],meta_results$`97.5%`[sort_order_sp0],meta_results$`97.5%`[sort_order_sp1]),
    is.summary = rep(FALSE,15),
    zero = 1,
    xlab = "Accuracy",
    clip=c(0.7, 1),
    #hrzl_lines = gpar(col="#444444"),
    col=fpColors(box=c("royalblue","darkblue","red","darkred")),
    fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI, fpDrawDiamondCI, fpDrawPointCI),
    boxsize = .2,
    line.margin = .1,
  )
  print(fp) #workaround needed when embedded in a for loop
  dev.off()
  
  
  #### COHORT MISCLASSIFICATION ###
  
  #extract parameters of interest  
  log_OR = log(as.matrix(opt.JAGS)[, "OR.meta"])
  cohort_misclass = as.matrix(opt.JAGS)[, "cohort_misclass"]
  
  cor(log_OR, cohort_misclass)
  #plot(log_OR, cohort_misclass) #uncomment for scatterplot but may be slow
  
  
  ### ADDITIONAL PLOTS and DIAGNOSTICS ###
  
  summary(opt.JAGS)
  
  write.csv(round(MCMCsummary(opt.JAGS),3), file=paste(output_dir,"Table_results_summary_",scenario_title[scenario],".csv",sep=""), row.names=T, na="")
  
  pdf(paste(output_dir,"Supplement_traceplots_",scenario_title[scenario],".pdf",sep="")) 
  plot(NULL, xlim = c(0, 10), ylim = c(0, 10), xlab = "", ylab = "", axes=F, main=paste("Traces and Prior-to-Posterior Updating: ", scenario_title[scenario], sep=""), cex.main = 0.8)
  MCMCtrace(opt.JAGS, params=c("OR.1", "OR.2", "OR.3"), pdf=F)
  MCMCtrace(opt.JAGS, params=c("OR.4", "OR.5", "OR.6"), pdf=F)
  MCMCtrace(opt.JAGS, params=c("OR.7", "OR.8", "OR.9"), pdf=F)
  MCMCtrace(opt.JAGS, params=c("OR.10", "OR.11", "OR.12"), pdf=F)
  MCMCtrace(opt.JAGS, params=c("OR.13", "OR.14", "OR.15"), pdf=F)
  
  priors_n = 15000
  MCMCtrace(opt.JAGS, params=c("SN0.1", "SP0.1"), priors=data.frame(rbeta(priors_n, SN0_cohort$shape1, SN0_cohort$shape2), rbeta(priors_n, SP0_cohort$shape1, SP0_cohort$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN1.1", "SP1.1"), priors=data.frame(rbeta(priors_n, SN1_cohort$shape1, SN1_cohort$shape2), rbeta(priors_n, SP1_cohort$shape1, SP1_cohort$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN0.2", "SP0.2"), priors=data.frame(rbeta(priors_n, SN0_cohort$shape1, SN0_cohort$shape2), rbeta(priors_n, SP0_cohort$shape1, SP0_cohort$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN1.2", "SP1.2"), priors=data.frame(rbeta(priors_n, SN1_cohort$shape1, SN1_cohort$shape2), rbeta(priors_n, SP1_cohort$shape1, SP1_cohort$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN0.3", "SP0.3"), priors=data.frame(rbeta(priors_n, SN0_cohort$shape1, SN0_cohort$shape2), rbeta(priors_n, SP0_cohort$shape1, SP0_cohort$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN1.3", "SP1.3"), priors=data.frame(rbeta(priors_n, SN1_cohort$shape1, SN1_cohort$shape2), rbeta(priors_n, SP1_cohort$shape1, SP1_cohort$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN0.4", "SP0.4"), priors=data.frame(rbeta(priors_n, SN0_cohort$shape1, SN0_cohort$shape2), rbeta(priors_n, SP0_cohort$shape1, SP0_cohort$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN1.4", "SP1.4"), priors=data.frame(rbeta(priors_n, SN1_cohort$shape1, SN1_cohort$shape2), rbeta(priors_n, SP1_cohort$shape1, SP1_cohort$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN0.5", "SP0.5"), priors=data.frame(rbeta(priors_n, SN0_cc$shape1, SN0_cc$shape2), rbeta(priors_n, SP0_cc$shape1, SP0_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN1.5", "SP1.5"), priors=data.frame(rbeta(priors_n, SN1_cc$shape1, SN1_cc$shape2), rbeta(priors_n, SP1_cc$shape1, SP1_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN0.6", "SP0.6"), priors=data.frame(rbeta(priors_n, SN0_cc$shape1, SN0_cc$shape2), rbeta(priors_n, SP0_cc$shape1, SP0_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN1.6", "SP1.6"), priors=data.frame(rbeta(priors_n, SN1_cc$shape1, SN1_cc$shape2), rbeta(priors_n, SP1_cc$shape1, SP1_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN0.7", "SP0.7"), priors=data.frame(rbeta(priors_n, SN0_cc$shape1, SN0_cc$shape2), rbeta(priors_n, SP0_cc$shape1, SP0_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN1.7", "SP1.7"), priors=data.frame(rbeta(priors_n, SN1_cc$shape1, SN1_cc$shape2), rbeta(priors_n, SP1_cc$shape1, SP1_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN0.8", "SP0.8"), priors=data.frame(rbeta(priors_n, SN0_cc$shape1, SN0_cc$shape2), rbeta(priors_n, SP0_cc$shape1, SP0_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN1.8", "SP1.8"), priors=data.frame(rbeta(priors_n, SN1_cc$shape1, SN1_cc$shape2), rbeta(priors_n, SP1_cc$shape1, SP1_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN0.9", "SP0.9"), priors=data.frame(rbeta(priors_n, SN0_cc$shape1, SN0_cc$shape2), rbeta(priors_n, SP0_cc$shape1, SP0_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN1.9", "SP1.9"), priors=data.frame(rbeta(priors_n, SN1_cc$shape1, SN1_cc$shape2), rbeta(priors_n, SP1_cc$shape1, SP1_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN0.10", "SP0.10"), priors=data.frame(rbeta(priors_n, SN0_cc$shape1, SN0_cc$shape2), rbeta(priors_n, SP0_cc$shape1, SP0_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN1.10", "SP1.10"), priors=data.frame(rbeta(priors_n, SN1_cc$shape1, SN1_cc$shape2), rbeta(priors_n, SP1_cc$shape1, SP1_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN0.11", "SP0.11"), priors=data.frame(rbeta(priors_n, SN0_cc$shape1, SN0_cc$shape2), rbeta(priors_n, SP0_cc$shape1, SP0_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN1.11", "SP1.11"), priors=data.frame(rbeta(priors_n, SN1_cc$shape1, SN1_cc$shape2), rbeta(priors_n, SP1_cc$shape1, SP1_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN0.12", "SP0.12"), priors=data.frame(rbeta(priors_n, SN0_cc$shape1, SN0_cc$shape2), rbeta(priors_n, SP0_cc$shape1, SP0_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN1.12", "SP1.12"), priors=data.frame(rbeta(priors_n, SN1_cc$shape1, SN1_cc$shape2), rbeta(priors_n, SP1_cc$shape1, SP1_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN0.13", "SP0.13"), priors=data.frame(rbeta(priors_n, SN0_cc$shape1, SN0_cc$shape2), rbeta(priors_n, SP0_cc$shape1, SP0_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN1.13", "SP1.13"), priors=data.frame(rbeta(priors_n, SN1_cc$shape1, SN1_cc$shape2), rbeta(priors_n, SP1_cc$shape1, SP1_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN0.14", "SP0.14"), priors=data.frame(rbeta(priors_n, SN0_cc$shape1, SN0_cc$shape2), rbeta(priors_n, SP0_cc$shape1, SP0_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN1.14", "SP1.14"), priors=data.frame(rbeta(priors_n, SN1_cc$shape1, SN1_cc$shape2), rbeta(priors_n, SP1_cc$shape1, SP1_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN0.15", "SP0.15"), priors=data.frame(rbeta(priors_n, SN0_cc$shape1, SN0_cc$shape2), rbeta(priors_n, SP0_cc$shape1, SP0_cc$shape2)), pdf=F)
  MCMCtrace(opt.JAGS, params=c("SN1.15", "SP1.15"), priors=data.frame(rbeta(priors_n, SN1_cc$shape1, SN1_cc$shape2), rbeta(priors_n, SP1_cc$shape1, SP1_cc$shape2)), pdf=F)
  
  MCMCtrace(opt.JAGS, params=c("OR.meta"), pdf=F)
  dev.off()
  
}

