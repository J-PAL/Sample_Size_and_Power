*Power by Simulation for an experiment without clusters

*Created by: Sabhya Gupta with input from Jack Cavanagh, Mike Gibson, Sarah Kopper
*All errors are the author's alone

*To report errors/make suggestions or ask questions contact: Sabhya Gupta - sagupta@povertyactionlab.org
*Last edited: 06/02/2021

********************************************************************************

/*

Our model is:

outcome = intercept + effect*treatment_dummy + individual_error

The error is assumed to be normally distributed, though the researcher can specify their own distribution
 

This do-file:
*Generates fake data for control and treatment with pre-determined effect size and error distribution
*Runs the regression
*Repeats this over and over and record how many times the treatment effect is statistically significant 
*The percentage of simulations with a statistically significant effect is the study's power



Contents:

	0. Housekeeping and specify temp files 
	1. SPECIFY design factors and the number of simulations
	2. Simulation on generated data
	3. Load results of simulation and estimate power
	
	
To use this file for your own power calculations, change the dataset/variables/values marked as "SPECIFY" 
The seed and all the locals in section 1 can be changed based on the research design 


*/

****************************************************************************************
***************************** 0. Housekeeping and specify temp files *******************
****************************************************************************************


cd "" 																			//SPECIFY - working directory
capture log close "power_by_simulation_no_clusters"
log using "power_by_simulation_no_clusters", replace
 

*Create a temporary file that will store the results of the simulation
tempname sim_name
tempfile sim_results
postfile `sim_name' reject_t using `sim_results'

clear
set seed 123456																	//SPECIFY - seed for replication 


******************************************************************************************
************* 1. SPECIFY design factors and the number of simulations ********************
******************************************************************************************

local sample_size = 500															//SPECIFY - sample size
local effect = 0.2																//SPECIFY - hypothesized treatment effect
local prop = 0.5 																//SPECIFY - ratio of the size of the treatment and control
local alpha = 0.05																//SPECIFY - define the significance value. 
local side "two"																//SPECIFY - "two", "left", "right" for a two-sided, left-sided or right-sided test

local sims=1000																	//SPECIFY - number of simulations


******************************************************************************************
******************************* 2. Simulation on generated data **************************
******************************************************************************************

*Generate fake data with specified distribution and effect, regress outcome on treatment and record if significant
local it = 1																	//iteration number

while `it' <=`sims'{
		clear
		quietly set obs `sample_size'
		quietly gen treat= runiform()<`prop'									//random assignment

		drawnorm error, m(0) sd(1)												//SPECIFY - the error is assumed to follow a standard normal but the researcher can decide the appropriate distribution for the sample population
		quietly gen outcome = 10 + error										//SPECIFY - assumed outcome for control. This can change based on the design specification
		
		quietly replace outcome = outcome + `effect' if treat					//potential outcome for the treatment is higher than the control by the effect size
		
		quietly regress outcome treat											//simple regression - the true effect is as specified above

		local t_value = _b[treat]/_se[treat]									//the t-value for the t-test
		local df=`sample_size'-2												//degrees of freedom is a function of the sample size
			
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
	

		post `sim_name' (`reject_t')											//write output from simulation to the temporary file
		
		tempfile simulated_data_no_clusters`it'									//save the data from the iterations
		save `simulated_data_no_clusters`it'', replace
	
		
		local it = `it' +1 

}

use `simulated_data_no_clusters4', clear										//load any of the simulated data files by changing the number - this may be useful to make sure that the simulation is working as intended



*****************************************************************************************
******************** 3. Load results of simulation and estimate power *******************
*****************************************************************************************


postclose `sim_name'
use `sim_results',clear

sum reject_t

di as error "If the treatment effect and the assumption about the distribution of the error is true, the study with sample size `sample_size' will detect the true treatment effect of `effect' with probability `r(mean)' in a `side'-sided test. Hence the power of the study is `r(mean)'"

cap log close

/* To improve the power, you can:
*increase sample size (sample_size)
*adjust the ratio of treatment to control (prop)
*increase significance level (alpha)
*add covariates to the regression (should be done carefully when simulating data)

*note that increasing the effect size will mechanically increase power, but to ensure a study is adequately 
powered it is important to use a reasonable effect size (i.e., the effect size specified here should not be 
increased just to make the calculations turn out in your favor). So, to increase power, you should only increase
the effect size specified if you believe you can increase the effect size in the real world, such as by tweaking
the intervention

********************************************************************************\
