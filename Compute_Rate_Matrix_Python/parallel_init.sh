#!/bin/bash
PARAM_CSV=parameter_files/paramValues.csv
MODEL_FILE=Compute_RateMatrix_MISAEx
# NOT IMPLEMENTED: MODEL_FILE is the name of the rate matrix calculation script to be called within the ../models/ folder. 

#PARAM_CSV_REALPATH=$(realpath $PARAM_CSV)

#. ../../.directory_save.txt
RESULTSDIR="../outputs_tmp"
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
#mkdir "$NEWFOLDER"/Analysis

cp ${PARAM_CSV} "${NEWFOLDER}"/paramValues.csv

time python ParallelWrapper.py -o "${NEWFOLDER}" -p "${PARAM_CSV}"
