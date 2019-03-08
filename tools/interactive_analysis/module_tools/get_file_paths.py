import configparser
import numpy as np
import pandas as pd


def load_config(cfg_file):
    config = configparser.ConfigParser()
    config.read(cfg_file)
    return config['PATHS']['simulation_folder_path']


def generate_set_paths(trial_path, paramSetNum):
    filename = 'set_{:05}.mat'.format(paramSetNum)
    paths = { 'simulation_folder': trial_path + '/',
              'analysis_folder': trial_path + '/Analysis/',
              'metrics_folder': trial_path + '/Analysis/Metrics/',
              'eigenvalues': trial_path + '/EigenValues/' + filename,
              'probvec': trial_path + '/ProbVec/' + filename,
              'ratematrix': trial_path + '/RateMatrix/' + filename,
              'stateslist': trial_path + '/StatesList/' + filename,
              'timescales': trial_path + '/TimeScales/' + filename,
              'tmat_properties': trial_path + '/TransitionMatrixProperties/' + filename,
              'simulation_parameters': trial_path + '/paramValues.csv',
              'phenotype_counts_csv': trial_path + '/Analysis/Metrics/landscape_phenotype_counts.csv',
              'phenotype_weights_csv': trial_path + '/Analysis/Metrics/landscape_phenotype_weights.csv',
              'system_entropies' : trial_path + '/Analysis/Metrics/system_entropies.txt'}
    return paths


def generate_trial_paths(trial_path):
    simulation_parameters = trial_path + '/paramValues.csv'
    df=pd.read_csv(simulation_parameters, index_col=0)
    filenames = ['set_{:05}.mat'.format(i) for i in range(1,df.index[-1]+1)]
    paths = { 'simulation_folder': trial_path + '/',
              'analysis_folder': trial_path + '/Analysis/',
              'metrics_folder': trial_path + '/Analysis/Metrics/',
              'eigenvalues': trial_path + '/EigenValues/',
              'probvec': trial_path + '/ProbVec/',
              'ratematrix': trial_path + '/RateMatrix/',
              'stateslist': trial_path + '/StatesList/',
              'timescales': trial_path + '/TimeScales/',
              'tmat_properties': trial_path + '/TransitionMatrixProperties/',
              'simulation_parameters': simulation_parameters,
              'phenotype_counts_csv': trial_path + '/Analysis/Metrics/landscape_phenotype_counts.csv',
              'phenotype_weights_csv': trial_path + '/Analysis/Metrics/landscape_phenotype_weights.csv',
              'system_entropies' : trial_path + '/Analysis/Metrics/system_entropies.txt'}
    return paths, df, filenames
