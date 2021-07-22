# Power_code
 Sample code for conducting power calculations using in-built commands and simulations

## Description

This repository contains some sample code on conducting power calculations in either STATA or R. All files are self contained and can be  run independently from the other scripts. Please read the code preamble for more details on each file.  


## Files

J-PAL_Power_built_in_commands: Uses in-built power commands in STATA and R to calculate sample size and minimum detectable effect size with or without covariates and with or without imperfect compliance in individual and clustered models. Both files can be run with any baseline dataset with a continuous outcome and binary treatment variable. See the code preamble for more instructions on how to adapt the code to your context. The sample code uses the Balsakhi dataset (baroda_0102_1obs.dta) for illustration purposes. 

J-PAL_Power_by_simulation_no_clusters: Calculates power using a dummy dataset simulated using an underlying sample distribution and a few design parameters. The underlying distribution and the design factors can be changed to suit the context of use. 

J-PAL_Power_by_simulation_clusters: Calculates power with simulated dataset as the previous file but with a clustered design. 

baroda_0102_1obs.dta: Data file required to run J-PAL_Power_built_in_command as written. The file can also be run with other similar datasets with a continuous outcome and treatment variable.You can learn more about the Balsakhi dataset from the documentation and data here at https://doi.org/10.7910/DVN/UV7ERB. 

## Support

Please use the [issue tracker] for all support requests


## License

See [license file]


