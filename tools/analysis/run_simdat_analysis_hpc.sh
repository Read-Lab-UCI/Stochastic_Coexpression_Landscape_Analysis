#!/bin/bash

module load anaconda/3.6-5.0.1

TRIALS_PATH="../../Simulation_Data-py/"
TRIAL_FOLDERS=("Trial_0028-py" "Trial_0029-py") 

for folder in "${TRIAL_FOLDERS[@]}"
do
	python make_simdat_analysis_hdf5.py -i "${TRIALS_PATH}${folder}"
done
