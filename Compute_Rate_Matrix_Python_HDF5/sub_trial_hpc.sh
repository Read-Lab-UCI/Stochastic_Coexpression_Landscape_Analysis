#!/bin/bash
#$ -S /bin/bash           # run with this shell
#$ -N MISAChromatin         # this name shows in qstat
#$ -q rxn                # run in this Q
#$ -pe openmp 32
#$ -j y                   # specify where the error messages get written to
#$ -cwd                   # run the job out of the current directory;

module load gcc/6.4.0 
module load Cluster_Defaults
module load openmpi-3.1.2/gcc-6.4.0
module load anaconda/3.6-5.0.1

PARAM_CSV=parameter_files/paramValues_MISAChromatin_long_N30.csv
RUNPARALLEL=true # set this to true or false to attempt parallel ratematrix calculations
MODEL_FILE=Compute_RateMatrix_MISAChromatin # Python filename that calculates ratematrix, called from models/ folder

RESULTSDIR="/pub/cgalliva/Simulation_Data-py"
FILENAME=Trial

FILES=$(find "$RESULTSDIR/" -maxdepth 1 -name "$FILENAME*" | sort | wc -l)
FILES="$(echo "$FILES" | sed -e 's/^[ \t]*//')"

if [[ "$FILES" != "0" ]] ; then
    echo Creating new trial folder
	
	LATEST=$(basename "$(find "$RESULTSDIR/" -maxdepth 1 -name "$FILENAME*" | sort | tail -1)") #Gets latest trial folder
	HEAD=${LATEST#*_} #Deletes everything before _
	NUM=${HEAD%-*} #Deletes everything after -
	NUM=$(echo $NUM | sed 's/^0*//') #Removes padded zeros
	
	((NUM=NUM + 1))
	NEWFOLDER="${RESULTSDIR}/"${FILENAME}_$(printf "%04d" "$NUM")-py
	mkdir "${NEWFOLDER}"

else
	echo Creating initial trial folder
	NEWFOLDER="${RESULTSDIR}/"${FILENAME}_0001-py
	mkdir "${NEWFOLDER}"
fi

# Making output file folders
#mkdir "$NEWFOLDER"/RateMatrix
#mkdir "$NEWFOLDER"/ProbVec
#mkdir "$NEWFOLDER"/Prob2D
#mkdir "$NEWFOLDER"/EigenValues
#mkdir "$NEWFOLDER"/TimeScales
mkdir "$NEWFOLDER"/Analysis

cp ${PARAM_CSV} "${NEWFOLDER}"/paramValues.csv

if [[ "$RUNPARALLEL" = true ]]; then
    python SimulationWrapper.py -o "${NEWFOLDER}" -p "${PARAM_CSV}" -pe -m ${MODEL_FILE}
else
    python SimulationWrapper.py -o "${NEWFOLDER}" -p "${PARAM_CSV}" -m ${MODEL_FILE}
fi
