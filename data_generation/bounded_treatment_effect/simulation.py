from tqdm import tqdm
import pandas as pd
import numpy as np
import fcit

B = 100
n = 1000
J = 9
I = J-1
iter = 100

name = "init_flipp"
methods = ["CausalVarEM"]

def create_ICA_output(estimator, dgp):
    data = dgp(noise_dict= {"loc" : 0, "scale" : 0}, prior= {"loc" : 0, "scale" : 1/np.sqrt(2)}, level_of_confounding = 1)
    for i in tqdm(range(B)):
        data.generate_data(n,I, J, random_state=i, bounded_treatment = True)
        
        pd.DataFrame(data.mixing_matrix_observed).to_csv(f"bounded_treatment_effect/data/true_mixing_{name}_{i}.csv", index = False)
        
        pd.DataFrame(data.signals).to_csv(f"bounded_treatment_effect/data/true_signals_{name}_{i}.csv", index = False)
        pd.DataFrame(data.data_observed).to_csv(f"bounded_treatment_effect/data/data_obs_{name}_{i}.csv", index = False)
        est = estimator(update_sigma=False, true_A=None, max_iter = iter,
                           random_seed= i,  mode = "lower_triangular",
                           init_range = [-3,3]) 
        est.fit(data.data_observed,J, noise_params= {"mean" : 0, "std" : 1}, progress_bar=False)

        pd.DataFrame(est.A).to_csv(f"bounded_treatment_effect/data/estimated_mixing_{"CausalVarEM"}_{name}_{i}.csv", index = False)
        pd.DataFrame(est.Signals).to_csv(f"bounded_treatment_effect/data/estimated_signals_{"CausalVarEM"}_{name}_{i}.csv", index = False)

        est = estimator(update_sigma=False, true_A=None, max_iter = iter,
                           random_seed= i,  mode = "VarEM",
                           init_range = [-3,3]) 
        est.fit(data.data_observed,J, noise_params= {"mean" : 0, "std" : 1}, progress_bar=False)

        pd.DataFrame(est.A).to_csv(f"bounded_treatment_effect/data/estimated_mixing_{"VarEM"}_{name}_{i}.csv", index = False)
        pd.DataFrame(est.Signals).to_csv(f"bounded_treatment_effect/data/estimated_signals_{"VarEM"}_{name}_{i}.csv", index = False)


def calculate_p_values(conditional_function, unconditional_function):
    print("calculate_p_values")
    for method in methods:
        p_values_unconditional = np.ones((B,J))
        p_values_conditional = np.ones((B,J))

        for i in tqdm(range(B)):
            data = pd.read_csv(f"bounded_treatment_effect/data/data_obs_init_flipp_{i}.csv", header=0).values
            signals = pd.read_csv(f"bounded_treatment_effect/data/estimated_signals_{method}_init_flipp_{i}.csv", header=0).values
            true_signals = pd.read_csv(f"bounded_treatment_effect/data/true_signals_init_flipp_{i}.csv", header=0).values
            p_values_unconditional[i,:] = unconditional_function(data, signals, n, J)
            p_values_conditional[i,:] = conditional_function(data, signals, n, J)

                
        pd.DataFrame(p_values_unconditional).reset_index().to_csv(f"bounded_treatment_effect/p_values_uncnditional_{method}.csv", header=False, index = False)
        
        pd.DataFrame(p_values_conditional).reset_index().to_csv(f"bounded_treatment_effect/p_values_conditional_{method}.csv", header=False, index = False)
    
    p_values_conditional_true_signals = np.ones((B,J))
    p_values_unconditional_true_signals = np.ones((B,J))

    for i in tqdm(range(B)):
        data = pd.read_csv(f"bounded_treatment_effect/data/data_obs_init_flipp_{i}.csv", header=0).values
        true_signals = pd.read_csv(f"bounded_treatment_effect/data/true_signals_init_flipp_{i}.csv", header=0).values
        p_values_conditional_true_signals[i,:] = conditional_function(data, true_signals, n, J)
        p_values_unconditional_true_signals[i,:] = unconditional_function(data, true_signals, n, J)
        
    pd.DataFrame(p_values_conditional_true_signals).reset_index().to_csv(f"different_confounding_levels/p_values_conditional_true_signals.csv", header=False, index = False)
    
    pd.DataFrame(p_values_unconditional_true_signals).reset_index().to_csv(f"different_confounding_levels/p_values_unconditional_true_signals.csv", header=False, index = False)
    