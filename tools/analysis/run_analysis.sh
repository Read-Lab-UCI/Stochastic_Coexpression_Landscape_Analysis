#!/bin/bash
TRIAL_PATH='/Users/camerongallivan/Research_Data/Simulation_Data/Trial_0002-py/'
RUNPARALLEL=true # set this to true or false to attempt parallel ratematrix calculations

if [[ "$RUNPARALLEL" = true ]]; then
    python hdf5_wrapper_plot_probability_and_quasipotential.py -p "${TRIAL_PATH}" -pe
else
    python hdf5_wrapper_plot_probability_and_quasipotential.py -p "${TRIAL_PATH}"
fi
