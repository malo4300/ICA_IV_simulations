# Combining ICA with an Instrument Validity Test  

This repository contains the code for my semester project, where I combine Independent Component Analysis (ICA) with an Instrumental Variable (IV) validity test to estimate the treatment effect in the presence of a hidden confounder.  

## Repository Structure  

### `data_generation/`  
Each simulation has its own folder, containing:  
- A Python script (`simulation.py`)  
- An R script (`iv_test_generation.R`)  

Each script implements a function that can be executed from either `IV_on_ICA.R` or `run_simulation.py`. Together, these scripts handle the entire pipeline for a given simulation, including:  
1. Data generation  
2. Running ICA  
3. Performing independence and IV tests  

### `treatment_effect_estimation/`  
This folder contains the treatment effect estimators, which rely on the datasets generated in `data_generation/`.  

### `data_generation/IV_code/`  
The IV test requires the code in this folder to run properly.  For my simulations it was provided by Patrick Burauel. 
