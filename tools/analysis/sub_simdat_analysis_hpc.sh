#!/bin/bash
#$ -S /bin/bash           # run with this shell
#$ -N analysis         # this name shows in qstat
#$ -q rxn                # run in this Q
#$ -j y                   # specify where the error messages get written to
#$ -cwd                   # run the job out of the current directory
#$ -m beas
#$ -M cgalliva@uci.edu

module load gcc/6.4.0 
module load Cluster_Defaults
module load anaconda/3.6-5.0.1

TRIALS_PATH="/pub/cgalliva/Simulation_Data-py/"
# TRIAL_FOLDERS=("Trial_0023-py" "Trial_0024-py" "Trial_0025-py" "Trial_0026-py")
TRIAL_FOLDERS=("Trial_0030-py") 

for folder in "${TRIAL_FOLDERS[@]}"
do
	python make_simdat_analysis_hdf5.py -i "${TRIALS_PATH}${folder}"
done
