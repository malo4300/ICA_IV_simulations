
####################################################################
# implementation of instrument validity test developed in 
# Burauel, Patrick F. "Evaluating Instrument Validity using the Principle of Independent Mechanisms." 
# Journal of Machine Learning Research 24.176 (2023): 1-56.
####################################################################

rm(list = ls()) # To clear all


#### define the simulation #######################

### specify the simulation that needs to be run ###############
simulation = "different_confounding_levels"
###########################



# load required packages
packages <- c("boot",
              "corpcor",
              "MASS",
              "parallel",
              "pracma",
              "knockoff")

lapply(packages,require,character.only=TRUE)

# load necessary functions
source('data_generation/IV_code/fn_2SLS.R')
source('data_generation/IV_code//estimate_kappa_bootstrap_wrapper.R')
source('data_generation/IV_code//fn_complement.R')
source('data_generation/IV_code//fn_complement_knockoff.R')
source('data_generation/IV_code//fn_test_instrument_validity.R')
source('data_generation/IV_code//estimate_confounding_via_kernel_smoothing.R')
source('data_generation/IV_code//estimate_confounding_sigmas.R')
source('data_generation/IV_code//estimate_confounding_via_kernel_smoothing.R')
source('data_generation/IV_code//estimate_confounding_sigmas.R')



# this line imports the iv_test_generation.R file implementing the calculate_iv_values function for each setting
source(paste0("data_generation/" , simulation,"/iv_test_generation.R" ))
#b number of datasets
B =100
# stop the timer
start_time <- Sys.time()  
calculate_iv_values(B = B)
end_time <- Sys.time()  

time_taken <- as.numeric(difftime(end_time, start_time, units = "hours"))
print(paste("Time taken:", time_taken, "hours"))


