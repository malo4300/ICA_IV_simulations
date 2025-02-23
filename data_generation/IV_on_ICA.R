
####################################################################
# implementation of instrument validity test developed in 
# Burauel, Patrick F. "Evaluating Instrument Validity using the Principle of Independent Mechanisms." 
# Journal of Machine Learning Research 24.176 (2023): 1-56.
####################################################################

rm(list = ls()) # To clear all

options(repos = c(CRAN = "https://cloud.r-project.org"))

#### define the simulation #######################
lib_path <- file.path("..", "..", "R_libs")

# Create the directory if it doesn't exist
if (!dir.exists(lib_path)) {
  dir.create(lib_path, recursive = TRUE)
}

install.packages("knockoff", lib = lib_path)
install.packages("pracma", lib = lib_path)
install.packages("MASS", lib = lib_path)
install.packages("corpcor", lib = lib_path)
install.packages("boot", lib = lib_path)
install.packages("parallel", lib = lib_path)


# load required packages
packages <- c("boot",
              "corpcor",
              "MASS",
              "parallel",
              "pracma",
              "knockoff")

lapply(packages,require,character.only=TRUE)

# load necessary functions
source('IV_code/fn_2SLS.R')
source('IV_code//estimate_kappa_bootstrap_wrapper.R')
source('IV_code//fn_complement.R')
source('IV_code//fn_complement_knockoff.R')
source('IV_code//fn_test_instrument_validity.R')
source('IV_code//estimate_confounding_via_kernel_smoothing.R')
source('IV_code//estimate_confounding_sigmas.R')
source('IV_code//estimate_confounding_via_kernel_smoothing.R')
source('IV_code//estimate_confounding_sigmas.R')

get_yes_or_no <- function(simulation) {
  repeat {
    user_input <- tolower(trimws(readline(prompt = paste0("Do you want do calculate to run the IV test for ", simulation, ". Please enter 'yes' or 'no': "))))
    if (user_input %in% c("yes", "no")) {
      return(user_input)
    } else {
      cat("Invalid input. Please type 'yes' or 'no'.\n")
    }
  }
}

# this line imports the iv_test_generation.R file implementing the calculate_iv_values function for each setting


### specify the simulation that needs to be run ###############
simulation = "bounded_treatment_effect"
###########################

source(paste0("data_generation/" , simulation,"/iv_test_generation.R" ))

#B number of datasets, can be reduced to test the pipline
B =1


if (get_yes_or_no(simulation) == "yes") {
  start_time <- Sys.time()
  calculate_iv_values(B = B)
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "hours"))
  print(paste("Time taken:", time_taken, "hours"))
} else {
  cat("Operation canceled.\n")
}




