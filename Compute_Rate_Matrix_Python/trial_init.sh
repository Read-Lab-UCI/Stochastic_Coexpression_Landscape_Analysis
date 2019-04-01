#!/bin/bash
PARAM_CSV=parameter_files/paramValues_MISAFlex_v2.csv
RUNPARALLEL=true # set this to true or false to attempt parallel ratematrix calculations
MODEL_FILE=Compute_RateMatrix_MISAFlex # Python filename that calculates ratematrix, called from models/ folder

RESULTSDIR="../Simulation_Results/"
FILENAME=Trial

FILES=$(find "$RESULTSDIR/" -maxdepth 1 -name "$FILENAME*" | sort | wc -l)
FILES="$(echo "$FILES" | sed -e 's/^[ \t]*//')"

if [ "$FILES" != "0" ] ; then
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
mkdir "$NEWFOLDER"/RateMatrix
mkdir "$NEWFOLDER"/ProbVec
mkdir "$NEWFOLDER"/Prob2D
mkdir "$NEWFOLDER"/EigenValues
mkdir "$NEWFOLDER"/TimeScales
mkdir "$NEWFOLDER"/Analysis

cp ${PARAM_CSV} "${NEWFOLDER}"/paramValues.csv

if [ "$RUNPARALLEL" = true ]; then
    python SimulationWrapper.py -o "${NEWFOLDER}" -p "${PARAM_CSV}" -pe -m $MODEL_FILE
else
    python SimulationWrapper.py -o "${NEWFOLDER}" -p "${PARAM_CSV}" -m $MODEL_FILE
fi
