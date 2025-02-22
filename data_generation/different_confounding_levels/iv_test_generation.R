library(progress)
Bootstraps = 1 # number of bootstrap draws
ncpus = 12
J = 9
n = 1000
synthetic_D_method = 'standard'
kappa_method = 'sigmas'


get_path = function(str, cn, i){
  return(paste(str ,cn, "_" ,i, ".csv", sep = ""))
}


order_data = function(data){
  ordered_data = matrix(0, nrow = nrow(data), ncol = ncol(data))
  ordered_data[,1] = data[, ncol(data)] # Y has to be the first column, in data it is the last
  ordered_data[,2:(ncol(data)-1)] = as.matrix(data[,1:(ncol(data)-2)]) # controls in the middle
  ordered_data[,ncol(data)] = data[, ncol(data)-1] # Treatment as last
  return(ordered_data)
}

p_values = function(ordered_data, signals){
  p_val = rep(0,ncol(signals))
  for (i in 1:ncol(signals)) {
    p_val[i] = fn_test_instrument_validity(ordered_data, signals[[names(signals)[i]]], Bootstraps, 
                                           ncpus, 
                                           kappa_method,
                                           synthetic_D_method)$pseudo_p
  }
  return(p_val)
}





p_values_ = function(method, cn, B,J){
  results = matrix(0, nrow = B, ncol = J+1)
  
  pb <- progress_bar$new(
    format = "[:bar] :percent | ETA: :eta",
    total = B, 
    clear = FALSE, 
    width = 60
  )
  
  for (i in 0:(B-1)) {
    dt_path = get_path("data_generation/different_confounding_levels/data/data_obs_large_conf_", cn, i)
    sg_path = get_path(paste0("data_generation/different_confounding_levels/data/estimated_signals_", method ,"_large_conf_"), cn, i)
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
    sg_path = get_path("data_generation/different_confounding_levels/data/true_signals_large_conf_", cn, i) 
    dt_path = get_path("data_generation/different_confounding_levels/data/data_obs_large_conf_", cn, i)
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
  
  cn = c(1,3,6)

  methods = cbind("VarEM", "CausalVarEM")
  
  for (i in 1:length(cn)) {
    for (j in 1:length(methods)) {
      results = p_values_(method = methods[j], cn = cn[i], B, J)
      write.csv(results, file = paste0("data_generation/different_confounding_levels/", methods[j], "_conf_", cn[i], ".csv"))
    }
    results = p_values_true_sources(cn[i], B, J)   
    write.csv(results, file = paste0("data_generation/different_confounding_levels/true_source_conf_", cn[i], ".csv"))
    
  }
  
}



