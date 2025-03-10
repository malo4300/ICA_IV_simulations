library(progress)
source("data_generation/IV_code/helper.R")
Bootstraps = 200 # number of bootstrap draws
ncpus = 10
n = 1000
synthetic_D_method = 'standard'
kappa_method = 'sigmas'


get_path = function(str, covs, i){
  return(paste(str ,covs, "_" ,i, ".csv", sep = ""))
}





p_values_ = function(method, covs, B,J){
  results = matrix(0, nrow = B, ncol = J+1)
  
  pb <- progress_bar$new(
    format = "[:bar] :percent | ETA: :eta",
    total = B, 
    clear = FALSE, 
    width = 60
  )
  
  for (i in 0:(B-1)) {
    dt_path = get_path("data_generation/different_numbers_of_covariates/data/data_obs_covs_", covs, i)
    sg_path = get_path(paste0("data_generation/different_numbers_of_covariates/data/estimated_signals_", method ,"_covs_"), covs, i)
    data =  read.csv(dt_path, header = 1)
    signals = read.csv(sg_path, header = 1) 
    ordered_data = order_data(data)
    results[i+1,1] = i
    results[i+1,2:ncol(results)] = p_values(ordered_data, signals)
    pb$tick()
  }
  return(results)
}

p_values_true_sources = function(covs, B,J){
  results = matrix(0, nrow = B, ncol = J+1)
  
  pb <- progress_bar$new(
    format = "[:bar] :percent | ETA: :eta",
    total = B, 
    clear = FALSE, 
    width = 60
  )
  
  for (i in 0:(B-1)) {
    sg_path = get_path("data_generation/different_numbers_of_covariates/data/true_signals_covs_", covs, i) 
    dt_path = get_path("data_generation/different_numbers_of_covariates/data/data_obs_covs_", covs, i)
    data =  read.csv(dt_path, header = 1)
    signals = read.csv(sg_path, header = 1)  +  matrix(rnorm(n*J,0,.1), n,J) # required for IV test to run
    ordered_data = order_data(data)
    results[i+1,1] = i
    results[i+1,2:ncol(results)] = p_values(ordered_data, signals)
    pb$tick()
  }
  return(results)
}





calculate_iv_values = function(B){
  
  covs = c(3,6,9)

  methods = cbind("VarEM", "CausalVarEM")
  
  for (i in 1:length(covs  )) {
    J = 3 + covs[i]
    for (j in 1:length(methods)) {
      results = p_values_(method = methods[j], covs = covs[i], B, J)
      write.csv(results, file = paste0("data_generation/different_numbers_of_covariates/", methods[j], "_covs_", covs[i], ".csv"))
    }
    results = p_values_true_sources(covs[i], B, J)   
    write.csv(results, file = paste0("data_generation/different_numbers_of_covariates/true_source_covs_", covs[i], ".csv"))
    
  }
  
}



