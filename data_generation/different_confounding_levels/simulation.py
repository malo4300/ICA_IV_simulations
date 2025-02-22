from tqdm import tqdm
import pandas as pd
import numpy as np
import fcit

B = 100
n = 1000
J = 9
I = J-1
iter = 100
level_of_confounding =[1,3,6]
def create_ICA_output(estimator, dgp):
    print("construct dataset and fit ICA")
    
    for cn in level_of_confounding:
        data_generator = dgp(noise_dict= {"loc" : 0, "scale" : 0}, prior= {"loc" : 0, "scale" : 1/np.sqrt(2)}, level_of_confounding = cn)
        for i in tqdm(range(B)):
            data_generator.generate_data(n, I, J, random_state=i, init_range = [-3,3])
            pd.DataFrame(data_generator.mixing_matrix_observed).to_csv(f"different_confounding_levels/data/true_mixing_matrix_cn_{cn}_{i}.csv", index = False)
            pd.DataFrame(data_generator.signals).to_csv(f"different_confounding_levels/data/true_signals_cn_{cn}_{i}.csv", index = False)
            pd.DataFrame(data_generator.data_observed).to_csv(f"different_confounding_levels/data/data_obs_cn_{cn}_{i}.csv", index = False)


            est = estimator(update_sigma=False, true_A=None, max_iter = iter,
                            random_seed= i,  mode = "lower_triangular",
                            init_range = [-3,3])
            est.fit(data_generator.data_observed,J, noise_params= {"mean" : 0, "std" : 1}, progress_bar=False)
            pd.DataFrame(est.A).to_csv(f"different_confounding_levels/data/estimated_mixing_CausalVarEM_cn_{cn}_{i}.csv", index = False)
            pd.DataFrame(est.Signals).to_csv(f"different_confounding_levels/data/estimated_signals_CausalVarEM_cn_{cn}_{i}.csv", index = False)

            est = estimator(update_sigma=False, true_A=None, max_iter = iter,
                            random_seed= i,  mode = "VarEM",
                            init_range = [-3,3])
            est.fit(data_generator.data_observed,J, noise_params= {"mean" : 0, "std" : 1}, progress_bar=False)
            pd.DataFrame(est.A).to_csv(f"different_confounding_levels/data/estimated_mixing_VarEM_cn_{cn}_{i}.csv", index = False)
            pd.DataFrame(est.Signals).to_csv(f"different_confounding_levels/data/estimated_signals_VarEM_cn_{cn}_{i}.csv", index = False)


def calculate_p_values(conditional_function, unconditional_function):
    print("calculate_p_values")
    methods = ["VarEM","CausalVarEM"]

    for method in methods:
        for cn in level_of_confounding:
            p_values_unconditional = np.ones((B,J))
            p_values_conditional = np.ones((B,J))
            for i in tqdm(range(B)):
                data = pd.read_csv(f"different_confounding_levels/data/data_obs_cn_{cn}_{i}.csv", header=0).values
                signals = pd.read_csv(f"different_confounding_levels/data/estimated_signals_{method}_cn_{cn}_{i}.csv", header=0).values
                p_values_unconditional[i,:] = unconditional_function(data, signals, n, J)
                p_values_conditional[i,:] = conditional_function(data, signals, n, J)
            pd.DataFrame(p_values_unconditional).reset_index().to_csv(f"different_confounding_levels/p_values_unconditional_{method}_{cn}", header=False, index = False)
            pd.DataFrame(p_values_conditional).reset_index().to_csv(f"different_confounding_levels/p_values_conditional_{method}_{cn}", header=False, index = False)


    