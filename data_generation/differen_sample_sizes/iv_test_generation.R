source("data_generation/IV_code/helper.R")
library(progress)
Bootstraps = 1 # number of bootstrap draws
ncpus = 30
J = 9

synthetic_D_method = 'standard'
kappa_method = 'sigmas'


get_path = function(str, n, i){
  return(paste(str ,n, "_" ,i, ".csv", sep = ""))
}




p_values_ = function(method, n, B,J){
  results = matrix(0, nrow = B, ncol = J+1)
  
  pb <- progress_bar$new(
    format = "[:bar] :percent | ETA: :eta",
    total = B, 
    clear = FALSE, 
    width = 60
  )
  
  for (i in 0:(B-1)) {
    dt_path = get_path("data_generation/differen_sample_sizes/data/data_obs_n:", n, i)
    sg_path = get_path(paste0("data_generation/differen_sample_sizes/data/estimated_mixing_CausalVarEM_", method ,"_n:"), n, i)
    data =  read.csv(dt_path, header = 1)
    signals = read.csv(sg_path, header = 1) 
    ordered_data = order_data(data)
    results[i+1,1] = i
    results[i+1,2:ncol(results)] = p_values(ordered_data, signals)
    pb$tick()
  }
  return(results)
}

p_values_true_sources = function(n, B,J){
  results = matrix(0, nrow = B, ncol = J+1)
  
  pb <- progress_bar$new(
    format = "[:bar] :percent | ETA: :eta",
    total = B, 
    clear = FALSE, 
    width = 60
  )
  
  for (i in 0:(B-1)) {
    sg_path = get_path("data_generation/differen_sample_sizes/data/true_signals_n:", n, i) 
    dt_path = get_path("data_generation/differen_sample_sizes/data/data_obs_n:", n, i)
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
  
  ns = c(100,1000,10000)

  methods = cbind("VarEM", "CausalVarEM")
  
  for (i in 1:length(ns)) {
    for (j in 1:length(methods)) {
      results = p_values_(method = methods[j], n = ns[i], B, J)
      write.csv(results, file = paste0("data_generation/differen_sample_sizes/", methods[j], "_n:", n[i], ".csv"))
    }
    results = p_values_true_sources(ns[i], B, J)   
    write.csv(results, file = paste0("data_generation/differen_sample_sizes/true_source_n:", n[i], ".csv"))
    
  }
  
}



