library(progress)
source("data_generation/IV_code/helper.R")
Bootstraps = 250 # number of bootstrap draws
ncpus = 12
J = 9
n = 1000
synthetic_D_method = 'standard'
kappa_method = 'sigmas'


get_path = function(str, i){
  return(paste(str ,i, ".csv", sep = ""))
}





p_values_ = function(method, B,J){
  results = matrix(0, nrow = B, ncol = J+1)
  
  pb <- progress_bar$new(
    format = "[:bar] :percent | ETA: :eta",
    total = B, 
    clear = FALSE, 
    width = 60
  )
  for (i in 0:(B-1)) {
    dt_path = get_path("data_generation/bounded_treatment_effect/data/data_obs_init_flipp_",i)
    sg_path = get_path(paste0("data_generation/bounded_treatment_effect/data/estimated_signals_", method ,"_init_flipp_"),  i)
    data =  read.csv(dt_path, header = 1)
    signals = read.csv(sg_path, header = 1) 
    ordered_data = order_data(data)
    results[i+1,1] = i
    results[i+1,2:ncol(results)] = p_values(ordered_data, signals)
    pb$tick()
  }
  return(results)
}

p_values_true_sources = function(cn, B,J){
  results = matrix(0, nrow = B, ncol = J+1)
  
  pb <- progress_bar$new(
    format = "[:bar] :percent | ETA: :eta",
    total = B, 
    clear = FALSE, 
    width = 60
  )
  
  for (i in 0:(B-1)) {
    sg_path = get_path("data_generation/bounded_treatment_effect/data/true_signals_init_flipp_", i) 
    dt_path = get_path("data_generation/bounded_treatment_effect/data/data_obs_init_flipp_",i)
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
  
  

  methods = cbind("CausalVarEM", "VarEM")
  
  
    for (j in 1:length(methods)) {
      results = p_values_(method = methods[j], B, J)
      write.csv(results, file = paste0("data_generation/bounded_treatment_effect/", methods[j],".csv"))
    }
    results = p_values_true_sources(cn[i], B, J)   
    write.csv(results, file = paste0("data_generation/bounded_treatment_effect/true_source.csv"))
    
  
  
}



