# Power by Simulation for an experiment with equal clusters
##################################################################################

# This R-script:
# - Generates fake data for control and treatment clusters with pre-determined cluster-specific errors, effect size and design features
# - Generates fake data for individuals with individual and cluster-specific errors
# - Runs the regression and clusters the errors
# - Repeats this over and over and record how many times the effect of the treatment is statistically significant 
# - The percentage of simulations with a statistically significant effect is the power
# 
# 
# In a clustered design, our model is:
#   
#   Outcome = intercept + effect*treatment_dummy + cluster_error + individual_error
# 
# The model can be changed based on the research design
# 
# The total variance of the outcome in the model = cluster-specific variance + individual specific variance
# We can define the ICC as the ratio of the variance of the cluster specific error and the total variance
# 
# The cluster-specific and individual-specific error are both assumed to be continuously and normally distributed but can be 
# changed based on the sample population. 
# 
# Contents:
#   0. Set-up
#   1. SPECIFY design factors
#   2. Simulation on generated data
#   3. Load results of simulation and estimate power
# 
# 
# - To use this file for your own power calculations, change the dataset/variables/values marked as "SPECIFY" 
# - The seed and all the locals in section 1 can be changed based on the research design 
# 
# Created by: Sabhya Gupta with input from Jack Cavanaugh, Mike Gibson, Sarah Kopper
# All errors are the author's alone

# To report errors/make suggestions or ask questions contact: Sabhya Gupta - sagupta@povertyactionlab.org / sabhya154@gmail.com
# Last edited: 05/19/2021

################################################################################


### Set-up ###########################################################

#install.packages(c("Hmisc","dplyr","ICC", "plm","lmtest"))

library(Hmisc)
library(dplyr)
library(ICC)
library(plm)
library(lmtest)
library(multiwayvcov)

rm(list=ls())

#### 1. Specify design factors ############################################

set.seed(123456)	                                      												#seed for replication 		
effect=0.25                                                                     #hypothesized treatment effect
prop=0.5                                                                        #ratio of the size of the treatment and control
alpha = 0.05                                                                    #define the significance value. Here we are conducting a two-sided test.
cluster_size = 500
num_clusters = 200                                                              #total number of clusters
sample_size = num_clusters*cluster_size
control_intercept = 10

ind_err_var = 10															                                  #This is the variation among individuals in a clusters
icc = 0.4																                                        #correlation between two randomly chosen individuals within a cluster
cluster_err_var = (icc*ind_err_var)/(1-icc)	      						                  #Variation across clusters. 
side="two"                                                                      #the kind of test (two, left or right)

sim.size = 1000                                                                 #number of simulations

#### 2. Calculate the power by Simulation ####################################

## initialize vectors and data frames to store results of the simulation
simulated_clusters <- data.frame()[1:num_clusters,]
simulated_data <- data.frame()[1:sample_size,]

rownames(simulated_clusters) <- 1:nrow(simulated_clusters)
rownames(simulated_data) <- 1:nrow(simulated_data)


reject_t <- rep(0, sim.size)
t_value <- rep(0L, sim.size)
rho_value<- rep(0L, sim.size)

critical_u <- qt(1-alpha/2,2*(num_clusters-1))                                      #the upper critical value
critical_l <- qt(alpha/2,2*(num_clusters-1))                                        #the lower critical value


i<- 1                                                                           #initialize iteration

while (i <= sim.size){
  simulated_clusters$cluster<- rep((1:num_clusters))
  
  simulated_clusters$cluster_error <- rnorm(n=num_clusters, mean=0, 
                                            sd=sqrt(cluster_err_var))           #SPECIFY distribution of cluster specific errors. 
                                                                                #The distribution can be changed based on the sample 
  simulated_clusters$treat <- ifelse(runif
                                     (n=num_clusters, min=0, max=1)<prop,1,0)   #random assignment of clusters
  
  simulated_data$u <- rnorm(n=sample_size, mean=0, sd=1)
  simulated_data$cluster <- as.numeric(cut2(simulated_data$u,
                                            g=num_clusters, m=cluster_size))    #dividing individuals into clusters
  
  full_data <- left_join(simulated_data, simulated_clusters, 
                         by=c("cluster"="cluster"),suffix = c(".x", ".y") )
  
  full_data$individual_error <-  rnorm(n=sample_size, 
                                       mean=0, sd=sqrt(ind_err_var))            #SPECIFY distribution of individual specific error
  
  full_data$outcome <- control_intercept+full_data$cluster_error+
                                        full_data$individual_error              #SPECIFY - outcome for the control
  
  full_data$outcome[full_data$treat==1]<-                                       #potential outcome for the treatment is higher than the control 
    full_data$outcome[full_data$treat==1] + effect                              #by the effect size

  control_subset <- subset(full_data,treat==0)
  rownames(control_subset) <- 1:nrow(control_subset)
  control_subset$cluster <- as.factor(control_subset$cluster)
  
  icc<- ICCest(cluster, outcome, data = control_subset)                         #testing the ICC in the simulation
  rho_value[i]<- icc$ICC                                                        #this is the ICC estimated from the simulated data.
                                                                                #As the sample size increase, this value will approach the ICC specified above 
  
  fit.sim <- lm(outcome~treat, data = full_data)                                #simple regression - the true effect is as specified above
  robust_SE <- cluster.vcov(fit.sim, full_data$cluster,df_correction = TRUE)    #cluster the standard errors
  robust_coef<-coeftest(fit.sim, robust_SE)                                     
  t_value[i]<- robust_coef[2,3]                                                 #the t value for the t test
  
  if (side == "two"){
    
    critical_u <- qt(1-alpha/2,2*(num_clusters-1))                              #the upper critical value
    critical_l <- qt(alpha/2,sample_size)                                       #the lower critical value
    reject_t[i] = ifelse(abs(t_value[i])> abs(critical_u),1,0)                  #reject if the t-value lies in the critical region 
    
  }
  
  else if (side == "right"){
    critical_u <- qt(1-alpha,2*(num_clusters-1))                                  #the upper critical value
    reject_t[i] = ifelse(abs(t_value[i])> abs(critical_u),1,0)                  #reject if the t-value is more than the upper critical value 
  }
  
  else if (side == "left"){
    critical_l <- qt(alpha,2*(num_clusters-1))                                  #the lower critical value
    reject_t[i] = ifelse(abs(t_value[i])> abs(critical_l),1,0)                  #reject if the abs value of t is higher than the lower critical level. 
    
  }
  
  i <- i +1 
  
  
}

### 3. Load results of simulation and estimate power ######################

power = mean(reject_t)
power

rho_calculated = mean(rho_value)

cat("If the treatment effect and the assumption about the disribution of the error is true, the study with", num_clusters,"of size", cluster_size,
    "will detect the true treatment effect of", effect, 
    "with probability,", power, "in a", side,"-sided test. Hence the power of the study is", power,
    ".The calculated ICC is", rho_calculated)

##### To improve the power, you can:
#   increase sample size (sample_size)
#   increase number of clusters (num_clusters). This will be more effective than increasing cluster size if the individuals are highly correlated
#   increase cluster size (cluster_size)
#   adjust the ratio of treatment to control (prop)
#   increase significance level (alpha)
#   add covariates to the regression (should be done carefully when simulating data)
# 
# note that increasing the effect size will mechanically increase power, but to ensure a study is adequately 
# powered it is important to use a reasonable effect size (i.e., the effect size specified here should not be 
# increased just to make the calculations turn out in your favor). So, to increase power, you should only increase
# the effect size specified if you believe you can increase the effect size in the real world, such as by tweaking
# the intervention




  