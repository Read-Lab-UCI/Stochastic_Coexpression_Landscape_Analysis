#!/bin/bash
TRIAL_PATH='../../Simulation_Results/Trial_0002-py/'
RUNPARALLEL=false # set this to true or false to attempt parallel ratematrix calculations

if [[ "$RUNPARALLEL" = true ]]; then
    python hdf5_wrapper_plot_probability_and_quasipotential.py -p "${TRIAL_PATH}" -pe
else
    python hdf5_wrapper_plot_probability_and_quasipotential.py -p "${TRIAL_PATH}"
fi
