#!/bin/bash
echo job_20
parallel --progress python script_plot_probability_and_quasipotential.py -c user_config_files/hpc_Trial_0020.cfg -s {} ::: {1..1296}
echo job_21
parallel --progress python script_plot_probability_and_quasipotential.py -c user_config_files/hpc_Trial_0021.cfg -s {} ::: {1..1296}
echo job_22
parallel --progress python script_plot_probability_and_quasipotential.py -c user_config_files/hpc_Trial_0022.cfg -s {} ::: {1..1296}
