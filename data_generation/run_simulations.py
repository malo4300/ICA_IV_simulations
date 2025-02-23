
from dgp.dgp import *
from ICA_estimator.estimators import *
from independence_tests.independence_tests import *

import time
# this import determines which simulation is run
from bounded_treatment_effect.simulation import create_ICA_output as create_ICA_output_bounded, calculate_p_values as calculate_p_values_bounded
from differen_sample_sizes.simulation import create_ICA_output as create_ICA_output_differen_sample_sizes, calculate_p_values as calculate_p_values_differen_sample_sizes



if __name__ == "__main__":
    # Run the simulation stop the time
    start = time.time()
    create_ICA_output_differen_sample_sizes(CausalVarEM, dgp_extended)
    calculate_p_values_differen_sample_sizes(p_vals_conditional, p_vals_unconditional)
    end = time.time()
    print("Time taken: ", (end - start)/3600, "hours")