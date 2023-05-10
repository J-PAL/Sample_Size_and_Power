# Power by Simulation for an experiment without clusters
##################################################################################
# 
# This R-script:
# - Generates fake data for control and treatment with pre-determined effect size and error distribution
# - Runs the regression
# - Repeats this over and over and record how many times the treatment effect is statistically significant 
# - The percentage of simulations with a statistically significant effect is the power
# 
# Our model is:
#   
#   outcome = intercept + effect*treatment_dummy + individual_error
# 
# We assume that the error is normally distributed but the researcher can specify their own distribution
# 
# 
# Contents:
#   
#   1. SPECIFY design factors
#   2. Simulation on generated data
#   3. Load results of simulation and estimate power
# 
# 
# - To use this file for your own power calculations, change the dataset/variables/values marked as "SPECIFY" 
# - The seed and all the locals in section 1 can be changed based on the research design 
# 
# Created by: Sabhya Gupta with input from Jack Cavanaugh, Mike Gibson, Sarah Kopper
# To report errors/make suggestions or ask questions contact: Sabhya Gupta - sagupta@povertyactionlab.org / sabhya154@gmail.com
# Last edited: 05/19/2021

###############################################################################


#### 1. SPECIFY design factors ####################################

## set seed, specify design factors and the number of simulations
set.seed(123456)	                                      												#seed for replication 		
sample_size = 500                                                               #sample size
effect = 0.2                                                                    #hypothesized treatment effect
prop = 0.5                                                                      #the proportion of the total sample in the treatment group
alpha = 0.05                                                                    #define the significance value. Here we are conducting a two-sided test.
sim.size = 2000                                                               #number of simulations
side = "two"                                                                    #the kind of test (two, left or right)


### 2. Calculate the power by Simulation ##########################

## initialize vectors and data frames to store results of the simulation
simulated_data <- data.frame()[1:sample_size,]


reject_t <- rep(0, sim.size)
t_value <- rep(0L, sim.size)

## The simulation loop

  i<- 1                                                                         #initialize iteration
  
  while (i <= sim.size){
    simulated_data$treat <- ifelse(runif(n=sample_size, min=0, max=1)<prop,1,0) #random assignment
    
    simulated_data$error <- rnorm(n=sample_size, mean=0, sd=1)                  #SPECIFY - simulated error is assumed to follow a normal distribution 
                                                                                #but the researcher can decide the appropriate distribution 
                                                                                #for the sample population. 
                                                                                #the sd is the expected standard deviation of the outcome
    
    simulated_data$outcome <- 10 + simulated_data$error                         #SPECIFY - assumed outcome for control based on the design specification
    
    simulated_data$outcome[simulated_data$treat==1]<-                           #potential outcome for the treatment is higher than the control 
                   simulated_data$outcome[simulated_data$treat==1] + effect     #by the effect size

    
     fit.sim <- lm(outcome~treat, data = simulated_data)                        #simple regression - the true effect is as specified above
     t_value[i]<- summary(fit.sim)$coef[2,1]/summary(fit.sim)$coef[2,2]                                    #the t value for the t test
     
     if (side == "two"){
       
       critical_u <- qt(1-alpha/2,(sample_size-2))                              #the upper critical value
       critical_l <- qt(alpha/2,(sample_size-1))                                #the lower critical value
       reject_t[i] = ifelse(abs(t_value[i])> abs(critical_u),1,0)               #reject if the t-value lies in the critical region 
  
     }
     
     else if (side == "right"){
       critical_u <- qt(1-alpha,(sample_size-2))                                #the upper critical value
       reject_t[i] = ifelse(t_value[i]> critical_u,1,0)                         #reject if the t-value is more than the upper critical value 
     }
     
     else if (side == "left"){
       critical_l <- qt(alpha,(sample_size-2))                                  #the lower critical value
       reject_t[i] = ifelse(t_value[i]< critical_l,1,0)                         #reject if the  t-value is lower than the lower critical level. 

     }
     
     i <- i +1 
  
  }

  
### 3. Load results of simulation and estimate power ######################
  
power = mean(reject_t)

cat("If the treatment effect and the assumption about the disribution of the error is true, the study with the sample size",
      sample_size, "will detect the true treatment effect of", effect, 
      "with probability", power, "in a", side,"-sided test. Hence the power of the study is", power)

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
# 
# the same is true of the ICC: a lower ICC will mechanically increase study power, but it's important to use a 
# reasonable ICC rather than one that is favorable in your calculations
