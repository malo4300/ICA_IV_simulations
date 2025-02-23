import numpy as np
import pandas as pd
def calc_error(df):
    """
    Calculate the error of the estimated treatment effect for each simulation"""
    true_treatment = df["true_treatment_effect"]
    seed = df["seed"]
    temp = df.drop(["level_of_confounding", "seed", "true_treatment_effect"], axis=1)
    err = temp.apply(lambda x: x - true_treatment, axis=0) 
    return err,seed


def calculate_conditioning_number_(path  = "xxx/mixing_matrix", seed = 0):
    """
    Calculate the conditioning number of the mixing matrix A"""
    A = pd.read_csv(f"{path}_{seed}.csv").values
    return np.linalg.cond(A)

def calculate_conditioning_number(path = "xxx/mixing_matrix_"):
    return [calculate_conditioning_number_(path, seed) for seed in range(100)]


def level_of_confounding(mm):
    I, J = mm.shape
    mm = mm.values
    treatment_effec = mm[I-1, J-2]
    confounder_effect_on_treatment = mm[I-2,0]
    total_confounder_effect = mm[I-1,0]
    confounder_edge = total_confounder_effect - treatment_effec*confounder_effect_on_treatment
    return np.abs(confounder_edge*confounder_effect_on_treatment)