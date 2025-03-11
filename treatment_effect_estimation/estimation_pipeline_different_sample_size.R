source("treatment_effect_estimation/treatment_effect_estimators/treatment_effect_estimators.R")

n = 1000
J = 9
I = J-1
treatment_col = I-1 
# change the functions for different mehtods


get_signals = function(seed){
    read.csv(get_path(paste0("data_generation/differen_sample_sizes/data/true_signals_n:", n, "_"), seed))
}
get_estimated_mixing_matrix = function(seed){
  read.csv(paste0("data_generation/differen_sample_sizes/data/estimated_mixing_CausalVarEM_n:", n, "_", seed, ".csv"))
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
cand_source_non_sense = vector("list", 100)
cand_source_iv_only =vector("list", 100)




for (i in 1:100) {
  candidates <- estimated_confounder_index_proper_test(p_values_iv[i,], p_values_indp_cond[i,])
  cand_confounder_idx[[i]] <- candidates
  candidates <- estimated_treatmet_and_outcome_ind(p_values_iv[i,], p_values_indp_undcond[i,])
  cand_source_idx[[i]] <- candidates
  candidates <- non_sense_method(p_values_iv[i,], p_values_indp_undcond[i,])
  cand_source_non_sense[[i]] <- candidates
  candidates <- iv_only(p_values_iv[i,])
  cand_source_iv_only[[i]] <- candidates
  
}


# confounder source -----


est = finding_confounder_source_estimation(cand_confounder_idx)

ind = est$ind

true_treatment_effect_confounder_idx = est$true_treatment_effect_confounder_idx
estimated_treatment_efect_confounder_idx =est$estimated_treatment_efect_confounder_idx
ols_biased = est$ols_biased
estimated_treatment_efect_column_extraction = est$estimated_treatment_efect_column_extraction
level_of_confounding = est$level_of_confounding


plot(true_treatment_effect_confounder_idx,estimated_treatment_efect_confounder_idx ,xlab = "True treatment effect", ylab = "Estimated treatment effect")
points(true_treatment_effect_confounder_idx, ols_biased, col ="red")
points(true_treatment_effect_confounder_idx, estimated_treatment_efect_column_extraction, col ="blue")
abline(a = 0, b = 1)


rmse(true_treatment_effect_confounder_idx, estimated_treatment_efect_confounder_idx)
rmse(true_treatment_effect_confounder_idx, ols_biased)
rmse(true_treatment_effect_confounder_idx, estimated_treatment_efect_column_extraction)

#save_treatment_estimation(list(seed = ind-1,
                               true_treatment_effect = true_treatment_effect_confounder_idx,
                               estimated_treatment_using_confounder_source = estimated_treatment_efect_confounder_idx,
                               ols_biased = ols_biased,
                               column_extraction = estimated_treatment_efect_column_extraction,
                               level_of_confounding = level_of_confounding), paste("true_sources_confounder_source_n:", n, sep =""))

##############################
# estimation on 7 sources -----


est = remove_treatment_and_outcome_estimation(cand_source_idx)
ind = est$ind
true_treatment_effect_confounder_idx = est$true_treatment_effect_confounder_idx
estimated_treatment_efect_source_idx = est$estimated_treatment_efect_source_idx
ols_biased = est$ols_biased
estimated_treatment_efect_column_extraction = est$estimated_treatment_efect_column_extraction
level_of_confounding = est$level_of_confounding




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




### non-sense 



####remove the most information about T 

est = remove_two_treatment_estimtaion(cand_source_non_sense)
ind = est$ind
true_treatment_effect_nonsense = est$true_treatment_effect_nonsense
two_treatment_effect = est$effect_est
ols_biased =est$ols_biased
estimated_treatment_efect_column_extraction = est$estimated_treatment_efect_column_extraction
level_of_confounding = est$level_of_confounding

plot(true_treatment_effect_nonsense,two_treatment_effect ,xlab = "True treatment effect", ylab = "Estimated treatment effect", ylim = c(-3,3))
points(true_treatment_effect_nonsense, ols_biased, col ="red")
points(true_treatment_effect_nonsense, estimated_treatment_efect_column_extraction, col ="blue")
abline(a = 0, b = 1)


rmse(true_treatment_effect_nonsense, two_treatment_effect)
rmse(true_treatment_effect_nonsense, ols_biased)
rmse(true_treatment_effect_nonsense, estimated_treatment_efect_column_extraction)


mean(abs(true_treatment_effect_nonsense-two_treatment_effect))
mean(abs(true_treatment_effect_nonsense-ols_biased))



save_treatment_estimation(list(seed = ind-1,
                               true_treatment_effect = true_treatment_effect_nonsense,
                               mistaken_scheme = two_treatment_effect,
                               ols_biased = ols_biased,
                               column_extraction = estimated_treatment_efect_column_extraction,
                               level_of_confounding = level_of_confounding), paste("true_sources_7_sources_nonsense_n:", n, sep =""))




### IV only extraction
est = iv_only_estimation(cand_source_iv_only)
ind = est$ind
true_treatment_effect_iv_only = est$true_treatment_effect_iv_only
iv_only_effect = est$iv_only_effect
ols_biased = est$ols_biased
estimated_treatment_efect_column_extraction = est$estimated_treatment_efect_column_extraction
level_of_confounding = est$level_of_confounding



plot(true_treatment_effect_iv_only,iv_only_effect ,xlab = "True treatment effect", ylab = "Estimated treatment effect", ylim = c(-3,3))
points(true_treatment_effect_iv_only, ols_biased, col ="red")
points(true_treatment_effect_iv_only, estimated_treatment_efect_column_extraction, col ="blue")
abline(a = 0, b = 1)


rmse(true_treatment_effect_iv_only, iv_only_effect)
rmse(true_treatment_effect_iv_only, ols_biased)
rmse(true_treatment_effect_iv_only, estimated_treatment_efect_column_extraction)


mean(abs(true_treatment_effect_iv_only-iv_only_effect))
mean(abs(true_treatment_effect_iv_only-ols_biased))



save_treatment_estimation(list(seed = ind-1,
                               true_treatment_effect = true_treatment_effect_iv_only,
                               iv_only_effect = iv_only_effect,
                               ols_biased = ols_biased,
                               column_extraction = estimated_treatment_efect_column_extraction,
                               level_of_confounding = level_of_confounding), paste("true_sources_iv_only_n:", n, sep =""))



