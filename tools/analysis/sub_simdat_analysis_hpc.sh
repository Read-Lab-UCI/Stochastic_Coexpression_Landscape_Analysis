#!/bin/bash
#$ -S /bin/bash           # run with this shell
#$ -N analysis         # this name shows in qstat
#$ -q rxn                # run in this Q
#$ -j y                   # specify where the error messages get written to
#$ -cwd                   # run the job out of the current directory;

#module load gcc/6.4.0 
#module load Cluster_Defaults
#module load openmpi-3.1.2/gcc-6.4.0
#module load anaconda/3.6-5.0.1

TRIALS_PATH="/Users/camerongallivan/Research_Data/Simulation_Data/"
TRIAL_FOLDERS=("Trial_0021-py" "Trial_0022-py")

for folder in "${TRIAL_FOLDERS[@]}"
do
	python make_simdat_analysis_hdf5.py -i "${TRIALS_PATH}${folder}"
done
