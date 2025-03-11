import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

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



def load_and_calculate_errors(file_path, estimated_col, true_col):
    """Loads data and calculates errors."""
    data = pd.read_csv(file_path, header=0, index_col=0)
    return calc_erros(data[estimated_col], data[true_col])

def gather_statistics(errors, rmse, mae, std):
    """Extracts RMSE, MAE, and STD from error dictionary and appends to lists."""
    rmse.append(errors["RMSE"])
    mae.append(errors["MAE"])
    std.append(errors["STD"])


def plot_errors(err, rmse, mae, std, ax, level, name_of_simulation):
    """Creates boxplots and adds legend with RMSE, MAE, and STD."""
    df_long = err.melt(var_name="variable", value_name="value").dropna()
    sns.boxplot(data=df_long, y="value", hue="variable", ax=ax)
    
    ax.set_title(f"{name_of_simulation}:  {level}", fontsize=20)
    ax.set_ylabel("Absolute deviation from true treatment effect", fontsize=20)
    ax.set_xlabel("Methods", fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_ylim((0,8))
    palette = sns.color_palette(n_colors=err.shape[1])
    #handles = [
    #    mpatches.Patch(color=palette[i], label=f"RMSE = {rmse[i]:.2f}, MAE = {mae[i]:.2f}, Std = {std[i]:.2f}  : {err.columns[i]}")
    #    for i in range(len(err.columns))
    #]
    handles = [
        mpatches.Patch(color=palette[i], label=f"{err.columns[i]}")
        for i in range(len(err.columns))
    ]
    ax.legend(handles=handles, loc="upper center", fontsize=15)

from sklearn.linear_model import LinearRegression
def calculate_ols_estimate(path = "data_generation/different_confounding_levels/data/", name = "data_obs_large_conf_3_"):
    ols_est = []
    for i in range(100):
        data = pd.read_csv(f"../{path}{name}{i}.csv")
        est = LinearRegression()
        # regress last column on all other columns
        est.fit(data.iloc[:, :-1], data.iloc[:, -1])
        ols_est.append(est.coef_[-1])

    return ols_est

def calc_erros(est, true):
    err = est-true
    rmse = np.sqrt(np.mean(err**2))
    mae = np.mean(np.abs(err))
    std = np.std(err)
    return {"RMSE": rmse, "MAE": mae, "STD": std, "error": np.abs(err)}


def column_extraction(path = "data_generation/different_confounding_levels/data/", name  = "estimated_mixing_CausalVarEM_large_conf_1_"):
    est = []
    for i in range(100):
        data = pd.read_csv(f"../{path}{name}_{i}.csv")
        est.append(data.iloc[-1, -2])
    return est


def get_true_treatment_effect(path, name):
    true = []
    for i in range(100):
        data = pd.read_csv(f"../{path}{name}{i}.csv")
        true.append(data.iloc[-1, -2])
    return true

def print_method_summary(method_summary):
    for level in method_summary.keys():
        print(f"Level {level}")
        col = method_summary[level]["cases"].keys()
        cases = method_summary[level]["cases"].values()

        print(pd.DataFrame([method_summary[level]["rmse"], method_summary[level]["mae"], method_summary[level]["std"], cases], columns=col, index=["RMSE", "MAE", "STD", "Cases"]))



def plot_cases_vs_performance(method_summary, metric, simulation):
    fig, ax = plt.subplots(1, 3, figsize=(15, 5), sharey=True)


  


    for i, level in enumerate(method_summary.keys()):
        methods = method_summary[level]["cases"].keys()
        cases = method_summary[level]["cases"].values()
        rmse = method_summary[level][metric]
        palette = sns.color_palette("Set1", n_colors=len(methods))

        ax[i].scatter(cases, rmse, c=palette, label=methods)
    

        ax[i].set_title(f"{simulation}: {level}")


    legend_labels = list(methods)
    legend_colors = palette[:len(legend_labels)]
    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10) for color in legend_colors]


    fig.legend(handles, legend_labels, loc='upper center', ncol=len(legend_labels), bbox_to_anchor=(0.5, 1.1))
    plt.tight_layout()

    plt.show()