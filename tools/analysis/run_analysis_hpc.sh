#!/bin/bash
TRIAL_PATH='/dfs3/pub/cgalliva/Simulation_Data-py/Trial_0008-py'
RUNPARALLEL=true # set this to true or false to attempt parallel ratematrix calculations

echo Starting $TRIAL_PATH
if [[ "$RUNPARALLEL" = true ]]; then
    python hdf5_wrapper_plot_probability_and_quasipotential.py -p "${TRIAL_PATH}" -pe
else
    python hdf5_wrapper_plot_probability_and_quasipotential.py -p "${TRIAL_PATH}"
fi
