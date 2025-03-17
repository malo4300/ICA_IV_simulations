# ICA_IV_simulations


In 'data_generation' each simulation has its own folder with one python file 'simulation.py' and one R file 'iv_test_generation.R'. Each of them implement one function that can be run from the file 'IV_on_ICA.R' or 'run_simulation.py'. In combination they run everything for the specific simulation from generating the data, running ICA, doing the independence and IV tests. 

These datasets are the foundation of the treatment estimators implemented in 'treatment_effect_estimation'.
The IV test is not running without the code in the folder data_generation/IV_code. 