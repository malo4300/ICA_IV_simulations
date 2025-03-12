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
  
  # for the independence test, the only not independent sources should be the treatment and confounder source, pick the smallest two p-values if unique
  
  ordered_p_values =  sort(as.numeric(p_values_indp))
  if(ordered_p_values[2] == ordered_p_values[3]){
    warning("Independence test return non-unique candidates for treatment and confounde source")
    return(NA)
  }
  
  candidates = which(ordered_p_values[2]  >= as.numeric(p_values_indp))
  
  # find the max p-values for the candidates and remove this 
   
  mn = min(p_values_iv[candidates])
  
  final_candidates = candidates[p_values_iv[candidates] == mn] 
  
  if(length(final_candidates)>1){
    warning("Final candidate not unique")
    return(NA)
  }
  if(length(final_candidates) == 0){
    warning("No candidate: return NA")
    return(NA)
  }
  return(final_candidates)
}

estimated_confounder_index_proper_test = function(p_values_iv,p_values_indp){
  
  candidates = which(p_values_indp <.05)
    
  if (length(candidates) != 2) {
    warning("Independence test did not return exactly 2 significant p-values")
    return(NA)
  }
  
  out_idx = candidates[p_values_iv[candidates] <.05]
  
  if (length(out_idx) != 1) {
    warning("IV test did not return exactly one siginifcant p-value for the proposed candidates")
    return(NA)
  }
  return(out_idx)
    
}




estimated_treatmet_and_outcome_ind = function(p_values_iv,p_values_indp_undcond){
  
  mx = max(as.numeric(p_values_indp_undcond))
  candidate_outcome=  which(as.numeric(p_values_indp_undcond) == mx)
  if (length(candidate_outcome)!=1) {
    warning("candidate_outcome not unique")
    return(NA)
  }
  
  mx  = max(as.numeric(p_values_iv))
  candidate_treatment =  which(as.numeric(p_values_iv) == mx)
  if (length(candidate_treatment)!=1) {
    warning("candidate_treatment not unique")
    return(NA)
  }
  
  if(candidate_outcome == candidate_treatment) {
    warning("Final candidates are the same")
    return(NA)
  }
  
  return(c(candidate_outcome,candidate_treatment ))
}



estimated_treatmet_and_outcome_ind_proper_test = function(p_values_iv,p_values_indp_undcond){
  
  outcome_idx = which(p_values_indp_undcond>.05)
  if (length(outcome_idx)!=1) {
    warning("Independence test did not return exactly 1 significant p-values")
    return(NA)
  }
  
  conf_index = which(p_values_iv>.05)
  if (length(conf_index)!=1) {
    warning("IV test did not return exactly one siginifcant p-value for the proposed candidates")
    return(NA)
  }
  if(outcome_idx == conf_index) {
    warning("Final candidates are the same")
    return(NA)
    
  }
  return(c(outcome_idx,conf_index))
  
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
  
  mn = min(as.numeric(p_values_indp_undcond))
  candidate_outcome=  which(as.numeric(p_values_indp_undcond) == mn)
  if (length(candidate_outcome)>1) {
    warning("candidate_outcome not unique")
    return(NA)
  }
  
  mx  = max(as.numeric(p_values_iv))
  candidate_treatment =  which(as.numeric(p_values_iv) == mx)
  if (length(candidate_treatment)>1) {
    warning("candidate_treatment not unique")
    return(NA)
  }
  
  if(candidate_outcome == candidate_treatment) {
    warning("Final candidates are the same")
    return(NA)
  }
  
  return(c(candidate_outcome,candidate_treatment ))
}

iv_only = function(p_values_iv){
  mx = max(p_values_iv)
  
  if (sum(p_values_iv == mx)> 1){
    warning("Final candidate not unique")
    return(NA)
    
  }
  return(which(p_values_iv == mx))
}




iv_only_estimation = function(cand_source_iv_only){
    ind = which(sapply(cand_source_iv_only, function(x) !is.na(x)))
  
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
  
  
  ind = which(sapply(cand_source_non_sense, function(x) all(!is.na(x))))
  
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
  ind =  which(sapply(cand_source_idx, function(x) all(!is.na(x))))
  
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
