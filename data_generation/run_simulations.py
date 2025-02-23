
from dgp.dgp import *
from ICA_estimator.estimators import *
from independence_tests.independence_tests import *

import time
# this import determines which simulation is run
from bounded_treatment_effect.simulation import *




if __name__ == "__main__":
    # Run the simulation stop the time
    start = time.time()
    create_ICA_output(CausalVarEM, dgp_extended)
    calculate_p_values(p_vals_conditional, p_vals_unconditional)
    end = time.time()
    print("Time taken: ", (end - start)/3600, "hours")