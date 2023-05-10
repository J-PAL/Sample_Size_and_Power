# This R script:
#
# Computes sample size and effect size for a given power and treatment to control size ratio
# Includes variations with controls, clusters and take-up rate
# Uses the Balsakhi dataset to illustrate how to calculate power
# About the data: The Balsakhi program was a remedial education program that was conducted in  Indian schools to increase literacy and numeracy skills.
# You can learn more about the Balsakhi dataset from the documentation and data here at https://doi.org/10.7910/DVN/UV7ERB


# Variables:
# - Outcome of interest is in the "normalised total score." This is represented by: 
#     - "pre_totnorm" at baseline
#     - "post_totnorm" at the endline
# - Treatment: bal
#     - 0=control
#     - 1=treatment
# - Clustering variable (by school): divid

# Note: key inputs for calculating power like the mean and the standard deviation at baseline, icc etc. are calculated using the 
# specified dataset but they can also be specified manually
#

# Binary outcome variable:
# Use power.prop.test https://stat.ethz.ch/R-manual/R-devel/library/stats/html/power.prop.test.html to calculate the sample 
# size or Minimum effect size for a binary outcome variable. The command is not suitable for clustered data. 
# The section on covariates is not applicable to binary outcome variables due to the different
# model specification for binary variables. See McConnell and Vera-Hernandez (2015) 
# for a discussion of how the power calculations change with covariates when the outcome variable is binary


# About the file:

# - This file contains sample code for the following:
#   0. Housekeeping and load data
#   1. No covariates
#     1a. Sample size for a given effect size
#     1b. MDE for a given sample size
#   2. With covariates (not applicable to binary outcome variable)
#     2a. Sample size for a specified effect - with covariates 
#     2b. MDE for a given sample size - with covariates
#   3. Sample size with Partial Take-up
#   4. Clustered designs
#     4a. Compute number of clusters for a given effect size and size of cluster 
#     4b. Compute cluster size given the number of clusters and effect size 
#     4c. Compute effect size for a given cluster size and number of clusters
# 
# - To use this file for your own power calculations, change the dataset/variables/values marked as "SPECIFY" 
# 
# Created by: Sabhya Gupta with input from Jack Cavanaugh, Maya Duru, Mike Gibson, Sarah Kopper
# All errors are the author's alone


# Contact: sagupta@povertyactionlab.org
# Last edited: 07/20/2021


### 0. Housekeeping and load data ##################################

#It may take a long time to run the install command. You can try to install them separately

# install.packages(c("haven", "ICC", "randomizr", "multiwayvcov", "lmtest","devtools", "Hmisc", "gsubfn"), 
#      dependencies=TRUE, INSTALL_opts = c('--no-lock'))

library(devtools)

#remotes::install_github("tidyverse/magrittr")
#devtools::install_github('vikjam/pwrcalc')

# The pwrcalc package was developed by vikjam: https://github.com/vikjam/pwrcalc
# Documentation: https://pwrcalc.readthedocs.io/en/latest/intro.html 
# If devtools fails, download the package file from the github page and install 
# by following the instructions here: http://outmodedbonsai.sourceforge.net/InstallingLocalRPackages.html


# Once you installed the packages, load them to get started:
library(haven)
library(ICC)
library(randomizr)
library(multiwayvcov)
library(lmtest)
library(magrittr)
library(pwrcalc)
library(Hmisc)
library(gsubfn)


# Load the dataset and specify the outcome and treatment variables
# The balsakhi data is already included in the pwrcalc package

rm(list=ls())
data(balsakhi)                                                                  #SPECIFY - load the dataset                                                                  

dataset = balsakhi                                                              #SPECIFY - the dataset
outcome = dataset$pre_totnorm                                                   #SPECIFY - the outcome variable
treatment = dataset$bal                                                         #SPECIFY - treatment variable


### 1. No covariates ######################################################


# The following code assumes the unit of randomization is the same as the unit of observation 


### 1a. Sample size for a given minimum detectable effect size ############

power = 0.8                                                                     #SPECIFY - desired power
nratio = 1                                                                      #SPECIFY - the ratio of the size of the treatment group to control group
alpha =0.05                                                                     #SPECIFY - significance level
p = nratio/(1+nratio)

N_base <- function (dataset, outcome, treatment){
  
  baseline_mean <- mean(outcome, na.rm = TRUE)                                  #Record the mean of the outcome variable at baseline   
  baseline_sd <- sd(outcome, na.rm = TRUE)                                      #Record the standard deviation of the outcome variable at baseline   
  
  expected_effect = 0.3*baseline_sd                                             #The expected effect should be specified based on the intervention and the cost. 
                                                                                #Here it is 0.3 times the sd
  treated_mean <- expected_effect + baseline_mean
  
  base_model = twomeans(m1 = baseline_mean, m2 = treated_mean, sd = baseline_sd, 
                        nratio=nratio, power=power, sig.level = alpha)
  
  # Here we assume that the standard deviation does not change with the treatment but 
  # you can also specify different standard deviations for the control and treatment groups

  print(base_model)
  
  cat("n1 and n2 are the control and treatment sample sizes respectively. 
m1 and m2 are the control and treatment means \n\n")
  
  cat("We need a minimum treatment size of",base_model$n2,"and control size of", 
      base_model$n1, "to detect an effect of", 
      expected_effect, "with a probability of", 
      power,  "if the effect is true and the ratio of the treatment and control is",nratio)
  
  return(base_model)
}



base_model <- N_base(dataset, outcome, treatment)  

# Change the parameters of the function to see how the sample size changes


### 1b. MDE for a given sample size ####################################

power = 0.8                                                                     #SPECIFY - desired power
nratio = 1                                                                      #SPECIFY - the ratio of the size of the treatment group to control group
alpha =0.05                                                                     #SPECIFY - significance level
p = nratio/(1+nratio)

mde_base <- function (dataset, outcome, treatment, N){
  t_power = qt(power, df=N-2)
  t_alpha = qt(1-alpha/2, df=N-2)
  
  baseline_sd <- sd(outcome, na.rm = TRUE)                                      #Record the standard deviation of the outcome variable at baseline   
  
  mde <- (t_power + t_alpha) * sqrt(1 /(p*(1-p))) * sqrt(1 / N) * baseline_sd
  mde = round(mde, digits=2)
  print(mde)
  
  cat("Given our sample size of",N,
      "and ratio of treatment and control group as,",
      nratio, ",the effect needs to be higher than", 
      mde, "for us to detect it with a probability of",power)
  
  return(mde)
  
}

N = nrow(dataset)                                                               #SPECIFY N - this is taken from the specified dataset but can be specified by the researcher

mde <- mde_base(dataset, outcome, treatment, N)              

# Change the parameters of the function to see how the MDE changes

### 2. Adding covariates #################################################


# To see how potential controls affect power,  we would ideally have access to a sample data set 
# (e.g. historical or pilot data). With these data, we would want to:

#   1. Regress Y_i (the outcome) on X_i (the controls) 
#   2. Use the residual standard deviation of the outcome variable from this regression to evaluate 
#      how much variance is explained by the set of covariates we plan to include

# - In practice, this residual SD becomes the new SD we include in our parametric power calculations
# 
# With access to historical data, for example, this would involve regressing last year's test scores 
# on test scores from the year before. Using balsakhi data, this would be as follows.

# Note that this section is not applicable for power calculations with a binary outcome variable. 
# See McConnell and Vera-Hernandez 2015 for a discussion of covariates for binary outcomes

### 2a. Sample size for a given effect size - with covariates #############

power = 0.8                                                                     #SPECIFY - desired power
nratio = 1                                                                      #SPECIFY - the ratio of the size of the treatment group to control group
alpha =0.05                                                                     #SPECIFY - significance level
p = nratio/(1+nratio)


N_cov <- function(dataset,covariates, outcome, treatment){
  
  cov_list <- paste(cov, collapse = " + ")
  formula <- as.formula(paste("outcome ~ ",cov_list,sep = ""))
  fit <- lm(formula, data = dataset)
  summary(fit)
  
  res_baseline_sd <- sd(summary(fit)$residuals, na.rm=TRUE)
  res_baseline_sd                                                               #the new SD for power calculations
  
  
  expected_effect = 0.3*res_baseline_sd                                         #The expected effect should be specified based on the intervention and the cost. 
                                                                                #Here it is 0.3 times the sd
  
  baseline_mean <- mean(outcome, na.rm = TRUE)                                  #Record the mean of the outcome variable at baseline   
  
  treated_mean <- expected_effect + baseline_mean
  
  cov_model <- twomeans(m1 = baseline_mean, m2 = treated_mean, sd = res_baseline_sd, 
                              power=power, nratio= nratio, sig.level = alpha)
  
  # Here we assume that the standard deviation does not change with the treatment but 
  # you can also specify different standard deviations for the control and treatment groups
  
  print(cov_model)
  
  cat("We need a minimum treatment size of",cov_model$n2,
      "and control size of", cov_model$n1, "to detect an effect of", 
      expected_effect, "with a probability of", 
      power,  "if the effect is true and the ratio of the treatment and control is",nratio)
  
  return(cov_model)

}

cov= c("gender", "std", "sessiond")                                                  #SPECIFY- a vector of covariate names- use baseline values of covariates

cov_model <- N_cov(dataset, cov, outcome, treatment) 
# Change the parameters of the function to see how the sample size changes


### 2b. MDE for a given sample size - with covariates ###################

power = 0.8                                                                     #SPECIFY - desired power
nratio = 1                                                                      #SPECIFY - the ratio of the size of the treatment group to control group
alpha =0.05                                                                     #SPECIFY - significance level
p = nratio/(1+nratio)

cov= c("gender", "std", "sessiond")                                             #SPECIFY- a vector of covariate names - use baseline values of covariates

mde_cov <- function (dataset, outcome, covariates, treatment, N){
  t_power = qt(power, df=N-2)
  t_alpha = qt(1-alpha/2, df=N-2)
  
  cov_list <- paste(covariates, collapse = " + ")
  formula <- as.formula(paste("outcome ~ ",cov_list,sep = ""))
  fit <- lm(formula, data = dataset)
  summary(fit)
  
  res_baseline_sd <- sd(summary(fit)$residuals, na.rm=TRUE)
  res_baseline_sd                                                                #the new SD for power calculations
  
  mde_res <- (t_power + t_alpha) * sqrt(1 /(p*(1-p))) * sqrt(1 / N) * res_baseline_sd
  mde_res = round(mde_res, digits=2)
  
  print(mde_res)
  
  cat("Given our sample size of",N,
      "and ratio of treatment and control group as,",
      nratio, ",the effect needs to be higher than", 
      mde_res, "for us to detect it with a probability of",power)
  
  return(mde_res)
  
}


N = nrow(dataset)                                                               #SPECIFY N - this is taken from the specified dataset but can be specified by the researcher

mde_res <- mde_cov(dataset, outcome, cov, treatment, N)   

# Change the parameters of the function to see how the MDE changes


### 3.Sample size with Partial Take-up ################################

## when we have imperfect compliance in the treatment or the control group, 
## the expected effect is reduced by a factor of the effective take-up 
## effective take-up = take-up in treatment - take-up in control

power = 0.8                                                                     #SPECIFY - desired power
nratio = 1                                                                      #SPECIFY - the ratio of the size of the treatment group to control group
alpha =0.05                                                                     #SPECIFY - significance level
p = nratio/(1+nratio)


N_partial <- function(dataset, outcome, treatment, takeup_treat, takeup_control){
  
  baseline_mean <- mean(outcome, na.rm = TRUE)                                  #Record the mean of the outcome variable at baseline   
  baseline_sd <- sd(outcome, na.rm = TRUE)                                      #Record the standard deviation of the outcome variable at baseline   
  
  
  expected_effect = 0.3*baseline_sd                                             #The expected effect with perfect take-up should be specified based on the intervention and the cost. 
                                                                                #Here it is 0.3 times the sd
  
  tu = takeup_treat - takeup_control													                  #effective take-up
  effect_tu = expected_effect*tu											                          #effect size after adjusting for take-up. 
                                                                                #This will be the effect size you expect to measure with a true effect size of 
                                                                                #`effect' and a take-up rate of `tu'. effect_tu < effect for imperfect take-up rates. 
  treat_tu = baseline_mean + effect_tu                                          #treatment mean after adjusting for take-up
  
  partial_model <- twomeans(m1 = baseline_mean, m2 = treat_tu, nratio=nratio, sd = baseline_sd, 
                            power=, sig.level = alpha)
  
  # Here we assume that the standard deviation does not change with the treatment but 
  # you can also specify different standard deviations for the control and treatment groups
  
  print(partial_model)
  
  cat("we need a higher sample size to have the same power because 
  the expected effect has decreased by the take-up rate \n\n")
  
  
  cat("We need a minimum treatment size of",partial_model$n2,
    "and control size of", partial_model$n1, "to detect an effect of", 
    expected_effect, "with a probability of", 
    power, "if the ratio of the treatment and control is",
    nratio,"and take-up rate is", tu, "and the effect is true")
  
  return(partial_model)

}



takeup_treat = 0.75                                                             #SPECIFY - take up in the treatment
takeup_control =  0.2                                                           #SPECIFY - take up in the control

partial_model <- N_partial(dataset, outcome, treatment,takeup_treat, takeup_control )

# Change the parameters of the function to see how the sample size changes


### 4.Clustered RCTs ####################################################

# The code presented so far has been for when the unit of randomization is the same
# as the unit of observation. The following code is for clustered designs, when there are
# multiple units of observation contained in a single unit of randomization 
# (e.g., randomization is at the school level but outcomes measured at the student level)


# Calculate ICC 

# Here the ICC is calculated from the dataset but it can also be manually defined

baseline_subset <- subset(dataset, !is.na(outcome))                             #remove the NA values 

cluster_var_subset <- as.factor(baseline_subset$divid)                           #SPECIFY- change "divid" to the cluster variable
outcome_subset <- baseline_subset$pre_totnorm                                    #SPECIFY- change "pre_totnorm" to the outcome variable

icc <- ICCest(cluster_var_subset, outcome_subset, data = baseline_subset)
rho <- icc$ICC                                                                

rho

# rho = 0.1                                                                     #SPECIFY - Manually define rho if required

### 4a. The number of clusters given cluster size and effect size #############

power = 0.8                                                                     #SPECIFY - desired power
nratio = 1                                                                      #SPECIFY - the ratio of the total number of units in the treatment and the control group
alpha =0.05                                                                     #SPECIFY - significance level

# maybe an iterative process is more appropriate since the k is a part of the df calc
number_of_clusters <- function(dataset, outcome, treatment, m, icc){
  
  baseline_mean <- mean(outcome, na.rm = TRUE)                                  #Record the mean of the outcome variable at baseline   
  baseline_sd <- sd(outcome, na.rm = TRUE)                                      #Record the standard deviation of the outcome variable at baseline   
  
  expected_effect = 0.3*baseline_sd                                              #The expected effect should be specified based on the intervention and the cost. 
                                                                                #Here it is 0.3 times the sd
  treated_mean = baseline_mean+expected_effect
  
  
  cluster_number <- twomeans(m1 = baseline_mean, m2 = treated_mean, sd = baseline_sd, nratio=nratio, 
                             sig.level = alpha, power=power)%>%
    clustered(obsclus = cluster_size, rho = icc)
  
  # Here we assume that the standard deviation does not change with the treatment but 
  # you can also specify different standard deviations for the control and treatment groups
  
  print(cluster_number)
  
  cat("Adjusted n1 and n2 indicate the sample size in the control and the treatment group respectively. 
    sample size is the total number of units across the clusters \n\n")
  
  
  cat("Given the size of each cluster as", 
     cluster_size,"and ratio of the number of units in the treatment to control as", nratio,
     "we need a minimum of", cluster_number$`Minimum number of clusters`, 
     "clusters to detect an effect of", 
     expected_effect,"with a probability of", power, "if the effect is true")
  
  return(cluster_number)
  
}

cluster_size = 50												                                        #SPECIFY - number of people in each cluster

cluster_number <- number_of_clusters(dataset, outcome, treatment, 
                                    cluster_size, rho)


### 4b. The cluster size given number of clusters and effect size ##############

# This assumes that the number of clusters in the control and the treatment arm is the same

power = 0.8                                                                     #SPECIFY - desired power
nratio = 2                                                                      #SPECIFY - the ratio of the total number of units in the treatment and the control group
alpha =0.05                                                                     #SPECIFY - significance level

cluster_size <- function(dataset, outcome, treatment, total_clusters, icc){
  
  baseline_mean <- mean(outcome, na.rm = TRUE)                                  #Record the mean of the outcome variable at baseline   
  baseline_sd <- sd(outcome, na.rm = TRUE)                                      #Record the standard deviation of the outcome variable at baseline   
  
  expected_effect = 0.3*baseline_sd                                             #The expected effect should be specified based on the intervention and the cost. 
                                                                                #Here it is 0.3 times the sd
  
  treated_mean = baseline_mean+expected_effect
  

  cluster_size_model <- twomeans(m1 = baseline_mean, m2 = treated_mean, sd = baseline_sd, 
                       power=power, nratio=nratio, sig.level= alpha)%>%
   clustered(numclus = total_clusters, rho = icc)
  
  print(cluster_size_model)
  
  cat("Given", total_clusters,"clusters, and the ratio of units in the treatment and the control as", nratio,
  "the minimum size of each cluster
    should be", cluster_size_model$`Average per cluster`, "for us to detect an effect
    of",expected_effect, "with a probability of", power, "if the effect is true")
  
  return(cluster_size_model)

}

total_clusters = 50                                                             #SPECIFY - the total number of clusters

cluster_size_model <- cluster_size(dataset, outcome, 
                                   treatment, total_clusters, rho)


### 4c. The effect size in a clustered design ##########################

# This assumes that the number of clusters in the control and the treatment arm is the same

power = 0.8                                                                     #SPECIFY - desired power
nratio = 1                                                                      #SPECIFY - the ratio of the total number of units in the treatment and the control group
alpha =0.05                                                                     #SPECIFY - significance level
p = nratio/(1+nratio)


mde_cluster <- function (dataset, outcome, treatment, cluster_size, number_clusters, icc){
  N = cluster_size*number_clusters
  t_power = qt(power, df=2*(number_clusters-1))
  t_alpha = qt(1-alpha/2, df=2*(number_clusters-1))
  
  t_stat <- t_power + t_alpha

  baseline_sd <- sd(outcome, na.rm = TRUE)                                      #Record the standard deviation of the outcome variable at baseline   
  
  cluster_mde <- t_stat * sqrt(1 / (p * (1 - p) * total_clusters )) * 
    sqrt(icc + (1 - icc) / cluster_size) * baseline_sd
  
  print(cluster_mde)
  
  cat("Given", total_clusters,"clusters and the ratio of the treatment and control size as", 
      nratio,"and", 
      cluster_size, "people in each cluster",
      "the effect needs to be higher than", cluster_mde,
      "for us to detect it with a probability of", power)
  
  return(cluster_mde)
  
}

total_clusters = 100                                                            #SPECIFY - the number of total number of clusters
cluster_size = 50												                                        #SPECIFY - number of people in each cluster

cluster_mde <- mde_cluster(dataset, outcome, treatment, cluster_size,
                                   total_clusters, rho)
