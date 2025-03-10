# the point is the combine the p-values of the IV test and of the independence test
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
get_path = function(str, i){
  return(paste(str,i, ".csv", sep = "" ))
}




estimated_confounder_index = function(p_values_iv,p_values_indp){
  # IV test should reject all sources as potential instruments leaving only the treatment source 
  # in practice pick the largest value if it is unique
  mx = max(p_values_iv)
  ind_iv_max = which(p_values_iv == mx)
  if(length(ind_iv_max)>1){
    warning("IV test has selcted non-unique candidate for treatment source")
  }
  
  # for the independence test, the only not independent sources should be the treatment and confounder source, pick the smallest two p-values if unique
  
  ordered_p_values =  sort(as.numeric(p_values_indp))
  if(ordered_p_values[2] == ordered_p_values[3]){
    warning("Independence test return non-unique candidates for treatment and confounde source")
  }
  
  candidates = which(ordered_p_values[2]  >= as.numeric(p_values_indp))

  # remove the source selected by the IV test
  final_candidate = candidates[!candidates %in% ind_iv_max]
  
  if(length(final_candidate)>1){
    warning("Final candidate not unique")
  }
  if(length(final_candidate) == 0){
    warning("No candidate: return NA")
    return(NA)
  }
  return(final_candidate)
}


estimated_confounder_index_v2 = function(p_values_iv,p_values_indp){
  
  # for the independence test, the only not independent sources should be the treatment and confounder source, pick the smallest two p-values if unique
  
  ordered_p_values =  sort(as.numeric(p_values_indp))
  if(ordered_p_values[2] == ordered_p_values[3]){
    warning("Independence test return non-unique candidates for treatment and confounde source")
  }
  
  candidates = which(ordered_p_values[2]  >= as.numeric(p_values_indp))
  
  # find the max p-values for the candidates and remove this 
   
  mn = min(p_values_iv[candidates])
  
  final_candidates = candidates[p_values_iv[candidates] == mn] 
  
  if(length(final_candidates)>1){
    warning("Final candidate not unique")
  }
  if(length(final_candidates) == 0){
    warning("No candidate: return NA")
    return(NA)
  }
  return(final_candidates)
}

estimate_treatment_index = function(p_values_iv,p_values_indp){
  
  # for the independence test, the only not independent sources should be the treatment and confounder source, pick the smallest two p-values if unique
  
  ordered_p_values =  sort(as.numeric(p_values_indp))
  if(ordered_p_values[2] == ordered_p_values[3]){
    warning("Independence test return non-unique candidates for treatment and confounde source")
  }
  
  candidates = which(ordered_p_values[2]  >= as.numeric(p_values_indp))
  
  # find the max p-values for the candidates and remove this 
  
  mn = min(p_values_iv[candidates])
  
  final_candidates = candidates[p_values_iv[candidates] != mn] 
  
  if(length(final_candidates)>1){
    warning("Final candidate not unique")
  }
  if(length(final_candidates) == 0){
    warning("No candidate: return NA")
    return(NA)
  }
  return(final_candidates)
}


online_extraction = function(p_values_indp, data_obs, signals){
  ordered_p_values =  sort(as.numeric(p_values_indp))
  if(ordered_p_values[2] == ordered_p_values[3]){
    warning("Independence test return non-unique candidates for treatment and confounde source")
  }
  
  candidates = which(ordered_p_values[2]  >= as.numeric(p_values_indp))
  
  # run iv test for the candidates an 
  

  order_data = function(data){
    ordered_data = matrix(0, nrow = nrow(data), ncol = ncol(data))
    ordered_data[,1] = data[, ncol(data)] # Y has to be the first column, in data it is the last
    ordered_data[,2:(ncol(data)-1)] = as.matrix(data[,1:(ncol(data)-2)]) # controls in the middle
    ordered_data[,ncol(data)] = data[, ncol(data)-1] # Treatment as last
    return(ordered_data)
  }
  
  
 
  
  B = 250 # number of bootstrap draws
  ncpus = 10
  synthetic_D_method = 'standard'
  kappa_method = 'sigmas'
  
  p_values = function(ordered_data, signals){
    p_val = rep(0,ncol(signals))
    for (i in 1:ncol(signals)) {
      p_val[i] = fn_test_instrument_validity(ordered_data, signals[[names(signals)[i]]], B, 
                                             ncpus, 
                                             kappa_method,
                                             synthetic_D_method)$pseudo_p
    }
    return(p_val)
  }
  ordered_data = order_data(data_obs)
  p_vals = p_values(ordered_data = ordered_data , signals[candidates])
  mn = min(p_vals)
  if (sum(p_vals == mn)> 1){
    warning("Final candidate not unique")
    
  }
  return(candidates[p_vals == mn])
}


estimated_treatmet_and_outcome_ind = function(p_values_iv,p_values_indp_undcond){
  

  candidate_outcome=  which.max(as.numeric(p_values_indp_undcond))
  candidate_treatment = which.max(as.numeric(p_values_iv))
  
  if(candidate_outcome == candidate_treatment) {
    warning("Final candidates are the same")
    
  }
  
  return(c(candidate_outcome,candidate_treatment ))
}

classic_ols = function(data_obsorved){
  I = ncol(data_obsorved)
  y = data_obsorved[,I]
  Treatment = data_obsorved[,I-1]
  controls = data_obsorved[,1:(I-2)]
  df = data.frame(y, Treatment = Treatment, cn = controls)
  fit = lm(y~.-1, df)
  return(coef(fit)["Treatment"])
}

column_extraction = function(mixing_matrix){
  J = ncol(mixing_matrix)
  I = nrow(mixing_matrix)
  return(mixing_matrix[I,J-1])
}



get_level_of_confounding = function(true_mixing_matrix){
  # level of confounding ca be extressed by the product of coefficents of the confounder on T and Y
  J = ncol(true_mixing_matrix)
  I = nrow(true_mixing_matrix)
  edge_U_T = true_mixing_matrix[I-1,1]
  treatment_effect =  as.numeric(column_extraction(true_mixing_matrix))
  edge_T_Y =  true_mixing_matrix[I,1] - treatment_effect*edge_U_T
  return(abs(edge_U_T*edge_T_Y))
}
  
  
rmse = function(true_treatment_effect, estimated_treatment_efect){
  return(sqrt(mean((true_treatment_effect-estimated_treatment_efect)^2)))
}




non_sense_method = function(p_values_iv,p_values_indp_undcond){
  
  # for the independence test, the only not independent sources should be the treatment and confounded source, pick the smallest two p-values if unique
  
  
  candidate_outcome=  which.min(as.numeric(p_values_indp_undcond))
  candidate_treatment = which.max(as.numeric(p_values_iv))
  # find the max p-values for the candidates and remove this 
  
  if(candidate_outcome == candidate_treatment) {
    warning("Final candidates are the same")
    
  }
  
  return(c(candidate_outcome,candidate_treatment ))
}

iv_only = function(p_values_iv){
  mx = max(p_values_iv)
  
  if (sum(p_values_iv == mx)> 1){
    warning("Final candidate not unique")
    
  }
  return(which(p_values_iv == mx))
}




iv_only_estimation = function(cand_source_iv_only){
  ind = which(sapply(cand_source_iv_only, function(x) length(x) == 1))
  
  l = length(ind)
  print(l)
  
  
  
  true_treatment_effect_iv_only = rep(0,l)
  iv_only_effect = rep(NA,l)
  ols_biased = rep(NA,l)
  estimated_treatment_efect_column_extraction = rep(NA,l)
  level_of_confounding = rep(NA, l)
  
  
  for (i in 1:l){
    seed = ind[i]-1
    indx = ind[i]
    signals = get_signals(seed)
    data = get_data(seed)
    remove = unlist(cand_source_iv_only[indx])
    df = data.frame(y = data[,I], treatment = data[, treatment_col], sign = signals[,-remove])
    fit = lm(y~.-1,df)  
    iv_only_effect[i] = coef(fit)["treatment"] 
    true_mm = get_true_mixing_matrix(seed)
    true_treatment_effect_iv_only[i] = column_extraction(true_mm)
    ols_biased[i] = classic_ols(data)
    mm = get_estimated_mixing_matrix(seed)
    estimated_treatment_efect_column_extraction[i] = column_extraction(mm)
    level_of_confounding[i] = get_level_of_confounding(true_mm)
    
  }
  return(data.frame(ind,
                    true_treatment_effect_iv_only,
               iv_only_effect,
               ols_biased,
               estimated_treatment_efect_column_extraction,
               level_of_confounding))
}



remove_two_treatment_estimtaion = function(cand_source_non_sense){
  
  
  ind = which(sapply(cand_source_non_sense, function(x) x[1] != x[2]))
  
  l = length(ind)
  print(l)
  
  
  true_treatment_effect_nonsense = rep(0,l)
  effect_est = rep(NA,l)
  ols_biased = rep(NA,l)
  estimated_treatment_efect_column_extraction = rep(NA,l)
  level_of_confounding = rep(NA, l)
  
  
  for (i in 1:l){
    seed = ind[i]-1
    indx = ind[i]
    signals = get_signals(seed)
    data = get_data(seed)
    remove = unlist(cand_source_non_sense[indx])
    df = data.frame(y = data[,I], treatment = data[, treatment_col], sign = signals[,-remove])
    fit = lm(y~.-1,df)  
    effect_est[i] = coef(fit)["treatment"] 
    true_mm = get_true_mixing_matrix(seed)
    true_treatment_effect_nonsense[i] = column_extraction(true_mm)
    ols_biased[i] = classic_ols(data)
    mm = get_estimated_mixing_matrix(seed)
    estimated_treatment_efect_column_extraction[i] = column_extraction(mm)
    level_of_confounding[i] = get_level_of_confounding(true_mm)
    
  }  
  return(data.frame(
    ind,
    true_treatment_effect_nonsense ,
    effect_est ,
    ols_biased ,
    estimated_treatment_efect_column_extraction ,
    level_of_confounding 
  ))
}


remove_treatment_and_outcome_estimation = function(cand_source_idx){
  ind = which(sapply(cand_source_idx, function(x) x[1] != x[2]))
  
  l = length(ind)
  print(l)
  
  
  true_treatment_effect_confounder_idx = rep(0,l)
  estimated_treatment_efect_source_idx = rep(NA,l)
  ols_biased = rep(NA,l)
  estimated_treatment_efect_column_extraction = rep(NA,l)
  level_of_confounding = rep(NA, l)
  
  
  
  for (i in 1:l){
    seed = ind[i]-1
    indx = ind[i]
    signals = get_signals(seed)
    data = get_data(seed)
    remove = unlist(cand_source_idx[indx])
    df = data.frame(y = data[,I], treatment = data[, treatment_col], sign = signals[,-remove])
    fit = lm(y~.-1,df)  
    estimated_treatment_efect_source_idx[i] = coef(fit)["treatment"] 
    true_mm = get_true_mixing_matrix(seed)
    true_treatment_effect_confounder_idx[i] = column_extraction(true_mm)
    ols_biased[i] = classic_ols(data)
    mm = get_estimated_mixing_matrix(seed)
    estimated_treatment_efect_column_extraction[i] = column_extraction(mm)
    level_of_confounding[i] = get_level_of_confounding(true_mm)
    
  }  
  
  return(data.frame(
    ind, 
    
    true_treatment_effect_confounder_idx ,
    estimated_treatment_efect_source_idx ,
    ols_biased ,
    estimated_treatment_efect_column_extraction ,
    level_of_confounding
  ))
  
}


finding_confounder_source_estimation= function(cand_confounder_idx){
  ind = which(sapply(cand_confounder_idx, function(x) length(x) == 1 && !is.na(x)))
  
  l = length(ind)
  print(l)
  
  true_treatment_effect_confounder_idx = rep(0,l)
  estimated_treatment_efect_confounder_idx = rep(NA,l)
  ols_biased = rep(NA,l)
  estimated_treatment_efect_column_extraction = rep(NA,l)
  level_of_confounding = rep(NA, l)
  
  
  
  for (i in 1:l){
    seed = ind[i]-1
    indx = ind[i]
    signals = get_signals(seed)
    data = get_data(seed)
    confounder_source_col = as.numeric(cand_confounder_idx[indx])
    df = data.frame(y = data[,I], treatment = data[, treatment_col], data[,c(-treatment_col, -I)], conf = signals[,confounder_source_col])
    fit = lm(y~.-1,df)  
    estimated_treatment_efect_confounder_idx[i] = coef(fit)["treatment"] 
    true_mm = get_true_mixing_matrix(seed)
    true_treatment_effect_confounder_idx[i] = column_extraction(true_mm)
    ols_biased[i] = classic_ols(data)
    mm = get_estimated_mixing_matrix(seed)
    estimated_treatment_efect_column_extraction[i] = column_extraction(mm)
    level_of_confounding[i] = get_level_of_confounding(true_mm)
  }  
  
  return(data.frame(
    ind,
    true_treatment_effect_confounder_idx ,
    estimated_treatment_efect_confounder_idx ,
    ols_biased ,
    estimated_treatment_efect_column_extraction ,
    level_of_confounding
    
  ))
}
