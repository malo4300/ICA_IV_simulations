source("treatment_effect_estimation/treatment_effect_estimators/treatment_effect_estimators.R")

n = 5000
# change the functions for different mehtods


get_signals = function(seed){
    read.csv(get_path(paste0("data_generation/differen_sample_sizes/data/true_signals_n:", n, "_"), seed))
}
get_estimated_mixing_matrix = function(seed){
  read.csv(paste0("data_generation/differen_sample_sizes/data/true_mixing_n:", n, "_", seed, ".csv"))
}

get_data = function(seed){
  read.csv(get_path(paste0("data_generation/differen_sample_sizes/data/data_obs_n:", n, "_"), seed))
}

get_true_mixing_matrix = function(seed){
  read.csv(get_path(paste0("data_generation/differen_sample_sizes/data/true_mixing_n:", n, "_"), seed))
}


save_treatment_estimation = function(list_of_df, name){
  write.csv(data.frame(list_of_df), file = paste("treatment_effect_estimation/treatment_effect_estimations/", name, ".csv", sep = ""))
}

p_values_iv = read.csv(paste("data_generation/differen_sample_sizes/true_source_n:",n,".csv", sep = ""), row.names = 1)[-1]
p_values_indp_cond = read.csv(paste("data_generation/differen_sample_sizes/p_values_conditional_true_sources_n:",n,".csv", sep = ""), row.names = NULL, header= FALSE)[-1]
p_values_indp_undcond = read.csv(paste("data_generation/differen_sample_sizes/p_values_unconditional_true_sources_n:", n, ".csv", sep = "") ,header= FALSE)[-1]


cand_confounder_idx <- vector("list", 100)
cand_source_idx = vector("list", 100)


for (i in 1:100) {
  if (i == 24 & n == 5000) {
    cand_confounder_idx[[i]] =
    cand_source_idx[[i]] =  NA
    
    next
  }
  candidates <- estimated_confounder_index_v2(p_values_iv[i,], p_values_indp_cond[i,])
  cand_confounder_idx[[i]] <- candidates
  candidates <- non_sense_method(p_values_iv[i,], p_values_indp_undcond[i,])
  cand_source_idx[[i]] <- candidates
}



# confounder source -----

ind = which(sapply(cand_confounder_idx, function(x) length(x) == 1 && !is.na(x)))

l = length(ind)
print(l)

true_treatment_effect_confounder_idx = rep(0,l)
estimated_treatment_efect_confounder_idx = rep(NA,l)
ols_biased = rep(NA,l)
estimated_treatment_efect_column_extraction = rep(NA,l)
level_of_confounding = rep(NA, l)
J = 9
I = J-1
treatment_col = I-1 # by construction
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


plot(true_treatment_effect_confounder_idx,estimated_treatment_efect_confounder_idx ,xlab = "True treatment effect", ylab = "Estimated treatment effect")
points(true_treatment_effect_confounder_idx, ols_biased, col ="red")
points(true_treatment_effect_confounder_idx, estimated_treatment_efect_column_extraction, col ="blue")
abline(a = 0, b = 1)


rmse(true_treatment_effect_confounder_idx, estimated_treatment_efect_confounder_idx)
rmse(true_treatment_effect_confounder_idx, ols_biased)
rmse(true_treatment_effect_confounder_idx, estimated_treatment_efect_column_extraction)

save_treatment_estimation(list(seed = ind-1,
                               true_treatment_effect = true_treatment_effect_confounder_idx,
                               estimated_treatment_using_confounder_source = estimated_treatment_efect_confounder_idx,
                               ols_biased = ols_biased,
                               column_extraction = estimated_treatment_efect_column_extraction,
                               level_of_confounding = level_of_confounding), paste("true_source_confounder_source_n:", n, sep =""))

##############################
# estimation on 7 sources -----


ind = which(sapply(cand_source_idx, function(x) x[1] != x[2]))

l = length(ind)
print(l)

true_treatment_effect_confounder_idx = rep(0,l)
estimated_treatment_efect_source_idx = rep(NA,l)
ols_biased = rep(NA,l)
estimated_treatment_efect_column_extraction = rep(NA,l)
level_of_confounding = rep(NA, l)


J = 9
I = J-1
treatment_col = I-1 # by construction
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


plot(true_treatment_effect_confounder_idx,estimated_treatment_efect_source_idx ,xlab = "True treatment effect", ylab = "Estimated treatment effect")
points(true_treatment_effect_confounder_idx, ols_biased, col ="red")
points(true_treatment_effect_confounder_idx, estimated_treatment_efect_column_extraction, col ="blue")
abline(a = 0, b = 1)


rmse(true_treatment_effect_confounder_idx, estimated_treatment_efect_source_idx)
rmse(true_treatment_effect_confounder_idx, ols_biased)
rmse(true_treatment_effect_confounder_idx, estimated_treatment_efect_column_extraction)



save_treatment_estimation(list(seed = ind-1,
                               true_treatment_effect = true_treatment_effect_confounder_idx,
                               estimated_treatment_on_sources = estimated_treatment_efect_source_idx,
                              ols_biased = ols_biased,
                               column_extraction = estimated_treatment_efect_column_extraction,
                               level_of_confounding = level_of_confounding), paste("true_sources_7_sources_n:", n, sep =""))




