** Power by Simulation for an experiment with equal clusters

*Created by: Sabhya Gupta with input from Jack Cavanagh, Mike Gibson, Sarah Kopper
*All errors are the author's alone

*To report errors/make suggestions or ask questions contact: Sabhya Gupta - sagupta@povertyactionlab.org
*Last edited: 06/02/2021

********************************************************************************

/* 
In a clustered design, the model (which can be changed based on the research design) is:

outcome = intercept + effect*treatment_dummy + cluster_error + individual_error

Total outcome variance = cluster-specific variance + individual-specific variance
The intra-cluster correlation coefficient (ICC) is the ratio of the variance of the cluster-specific error and the total variance

The cluster-specific and individual-specific error are both assumed to be continuously and normally distributed but can be 
changed based on the sample population. 


This do-file:
*Generates fake data for control and treatment clusters with pre-determined cluster-specific errors, effect size and design features
*Generates fake data for individuals with individual and cluster-specific errors
*Runs the regression and clusters the errors
*Repeats this over and over and record how many times the effect of the treatment is statistically significant 
*The percentage of simulations with a statistically significant effect is the power

It contains sample code for the following:
	0. Housekeeping and specify temp files 
	1. SPECIFY design factors and the number of simulations
	2. Simulation on generated data
	3. Load results of simulation and estimate power


To use this file for your own power calculations, change the dataset/variables/values marked as "SPECIFY" 
The seed and all the locals in section 1 can be changed based on the research design 
*/



********************************************************************************************************************
********************************* 0. Housekeeping and specify temp files  ******************************************
********************************************************************************************************************

cd "" 																			//SPECIFY - working directory
capture log close "power_by_simulation_clusters"
log using "power_by_simulation_clusters", replace
 
*Create a temporary file that will store the results of the simulation
tempname sim_name
tempfile sim_results
postfile `sim_name' reject_t rho_calculated  using `sim_results'

clear
set seed 123456																	//SPECIFY - seed for replication 



********************************************************************************************************************
************************* 1. SPECIFY design factors and the number of simulations **********************************
********************************************************************************************************************

local effect = 0.5																//SPECIFY - hypothesized treatment effect
local prop = 0.5 																//SPECIFY - ratio of the size of the treatment and control
local alpha = 0.05																//SPECIFY - the significance value. 
local side "two"																//SPECIFY   "two", "left", or "right" for a two-sided, left-sided or right-sided test

local control_intercept = 10													//SPECIFY - intercept for the design specification
local num_clusters = 100														//SPECIFY - total number of clusters
local cluster_size = 500														//SPECIFY - size of each cluster							
local sample_size = `num_clusters'*	`cluster_size'								//SPECIFY - sample size i.e. the number of units of observation across all clusters

local within_err_var = 1														//SPECIFY - this is the within-cluster variation
local icc = 0.4																	//SPECIFY - correlation between two randomly chosen individuals within a cluster
local between_err_var = (`icc'*`within_err_var')/(1-`icc')						//variation between clusters 

/*Cluster-specific variance is 0 when the ICC is 0. That is, when individuals in a cluster are not correlated with each other, 
the variance in the outcome is only from the variation across individuals. Cluster-specific variance increases as the ICC increases. 
The ICC can also be specified as a function of the within-cluster variance and between-cluster variance by solving for it in the 
between_err_var equation above. 

The ICC can be calculated using data from a similar population or based on the literature. 
If using existing data, the ICC can be estimated as follows:

loneway outcome cluster_var if treatment==0										//The loneway command calculates the one-way ANOVA by a group variable. 
																				//It gives the within-group variation and the between group variation of a variable. 
local rho = `r(rho)'															//rho gives the ICC 
*/
	

local sims=1000																	//SPECIFY number of simulations



********************************************************************************************************************
****************************************** 2. Simulation on generated data *****************************************
********************************************************************************************************************

*Generate fake data with clusters, specified distribution and effect, regress outcome on treatment and record if significant

local it = 1																	//iteration number

while  `it'<=`sims'{
	clear
	
	*cluster errors
	quietly set obs `num_clusters'
	quietly gen cluster_group= _n
	quietly gen cluster_error=rnormal(0,sqrt(`between_err_var'))				//we assume that the cluster-specific error is distributed normally with mean 0 
																				//SPECIFY - should be decided by the researcher based on the sample population

	quietly gen treat= runiform()<`prop'										//random assignment at the cluster level
	
	sort cluster_group
	tempfile cluster_error_g
	quietly save `cluster_error_g', replace										//save the file with cluster-specific errors

	
	*Assign individuals to clusters
	clear
	quietly set obs `sample_size'
	quietly gen u=invnormal(uniform())											//create a random number to randomly divide into clusters
	quietly egen cluster_group = cut(u), group(`num_clusters')					//divide individuals into clusters
	quietly replace cluster_group = cluster_group+1								//replace the numbering from 0-9 to 1-10
	sort cluster_group															//sort by cluster
	
	*merge in cluster errors
	quietly merge m:1 cluster_group using  `cluster_error_g'					//merge with the dataset with cluster-specific errors
	quietly gen individual_error=rnormal(0,sqrt(`within_err_var'))				//SPECIFY - individual specific error that is assumed to be normally distributed with mean 0. 
																				//This should be decided based on the sample population
																				
	quietly gen outcome = `control_intercept' + cluster_error + individual_error //SPECIFY - assumed outcome for control
	
	quietly replace outcome = outcome + `effect' if treat==1					//potential outcome for the treatment is higher than the control by the effect size
	
	quietly loneway outcome cluster_group if treat==0							//The loneway command calculates the one-way ANOVA by a group variable. 
																				//It gives the 	within-group variation and the between group variation of a variable. 
																				//It also  produces the intra-cluster correlation coefficient (ICC). 
																				//This is computed at the baseline or for the control (as a proxy for the baseline)
																				
	local rho_calculated=`r(rho)'												//this is the ICC estimated from the simulated data. As the sample size increase, this value will approach the ICC specified above 
	
	quietly regress outcome treat, vce(cluster cluster_group)					//simple regression - the true effect is as specified above. The errors are clustered
	
	local t_value = _b[treat]/_se[treat]
	local df=2*((`num_clusters'/2)-1)											//degrees of freedom is a function of the number of clusters
	
	if "`side'" == "two" {
		local critical_l = invt(`df', `alpha'/2)								//the lower critical value
		local critical_u = invt(`df', 1-`alpha'/2)								//the upper critical value
		local reject_t=(`t_value'>`critical_u')|(`t_value'<`critical_l')		//reject if the t-value lies in the critical level. It takes value 1 if the null is rejected and 0 if not
	}
	
	else if "`side'" == "left" {
		local critical_l = invt(`df', `alpha')									//the lower critical value
		local reject_t=(`t_value'<`critical_l')									//reject if the t-value is less than the lower critical level. It takes value 1 if the null is rejected and 0 if not
	}
	
	else if "`side'" == "right" {
		local critical_u = invt(`df', 1-`alpha')								//the upper critical value
		local reject_t=(`t_value'>`critical_u')									//reject if the t-value is more than the upper critical level. It takes value 1 if the null is rejected and 0 if not
	}
	
	
	
	quietly post `sim_name' (`reject_t') (`rho_calculated') 					//write output from simulation to the temporary file
	
	tempfile simulated_data_clusters`it'										//save the data from the iterations
	save `simulated_data_clusters`it'', replace
	
	local it=`it'+1
	clear
}

use `simulated_data_clusters4', clear											//load any of the simulated data files by changing the number - this may be useful to make sure that the simulation is working as intended


*****************************************************************************************
******************** 3. Load results of simulation and estimate power *******************
*****************************************************************************************

postclose `sim_name'
use `sim_results',clear

sum reject_t
local power = `r(mean)'

sum rho_calculated
local rho_calculated = round(`r(mean)',0.01)

di as error "If the treatment effect and the assumption about the distribution of the cluster and individual specific error is true, the study with `num_clusters' clusters of size `cluster_size' will detect the true treatment effect of `effect' with probability `power' in a `side'-sided test. Hence the power of the study is `power'. The calculated ICC from the simulation is `rho_calculated'"



/* To improve the power, you can:
*increase sample size (sample_size)
*increase number of clusters (num_clusters). This will be more effective than increasing cluster size if the individuals are highly correlated
*increase cluster size (cluster_size)
*adjust the ratio of treatment to control (prop)
*increase significance level (alpha)
*add covariates to the regression (should be done carefully when simulating data)

*note that increasing the effect size will mechanically increase power, but to ensure a study is adequately 
powered it is important to use a reasonable effect size (i.e., the effect size specified here should not be 
increased just to make the calculations turn out in your favor). So, to increase power, you should only increase
the effect size specified if you believe you can increase the effect size in the real world, such as by tweaking
the intervention
 
*the same is true of the ICC: a lower ICC will mechanically increase study power, but it's important to use a 
reasonable ICC rather than one that is favorable in your calculations

*/

cap log close
