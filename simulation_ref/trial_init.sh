#!/bin/bash
PARAM_CSV=parameter_files/paramValues_test.csv
MODEL_FILE=Compute_RateMatrix_MISAEx
# MODEL_FILE is the name of the rate matrix calculation script to be called within the ../models/ folder. 
  # Do not include the .m file extension

RESULTSDIR="../outputs_tmp/"
FILENAME=Trial

FILES=$(find "$RESULTSDIR" -maxdepth 1 -name "Trial*" | sort | wc -l)
FILES="$(echo "$FILES" | sed -e 's/^[ \t]*//')"

if [ "$FILES" != "0" ] ; then
    echo Creating new trial folder
	
	LATEST=$(basename "$(find "$RESULTSDIR" -maxdepth 1 -name "Trial*" | sort | tail -1)") #Gets latest trial folder
	HEAD=${LATEST#*_} #Deletes everything before _
	NUM=${HEAD%-*} #Deletes everything after -
	NUM=$(echo $NUM | sed 's/^0*//') #Removes padded zeros
	
	((NUM=NUM + 1))
	NEWFOLDER="${RESULTSDIR}"${FILENAME}_$(printf "%04d" "$NUM")
	mkdir "$NEWFOLDER"

else
	echo Creating initial trial folder
	NEWFOLDER="${RESULTSDIR}"${FILENAME}_0001
	mkdir "$NEWFOLDER"
fi

# Making output file folders
mkdir "$NEWFOLDER"/RateMatrix
mkdir "$NEWFOLDER"/ProbVec
mkdir "$NEWFOLDER"/Prob2D
mkdir "$NEWFOLDER"/Analysis
mkdir "$NEWFOLDER"/TimeScales
mkdir "$NEWFOLDER"/EigenValues

# Runs python script to generate python parameters for the trial and appends to .trialdir_save.txt
cp ${PARAM_CSV} "$NEWFOLDER"/paramValues.csv
echo $NEWFOLDER
#matlab -nodisplay -nodesktop -nosplash -r "SimulationWrapper(\"${NEWFOLDER}\", \"${MODEL_FILE}\")"; 
