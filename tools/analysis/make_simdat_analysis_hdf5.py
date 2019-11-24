import os
import h5py
import numpy as np
import pandas as pd
import argparse

# Suppressing numpy error outputs
old_settings = np.seterr(all='ignore')


def get_logic_parameter_columns(model_name):
    if model_name == 'MISAFlex':
        logic_columns = ['g0','g1','g2','g3']
    elif model_name == 'MISAFlex_asym':
        logic_columns = ['g0','g1','g2','g3','g0_b','g1_b','g2_b','g3_b']
    elif model_name == 'MISAChromatin':
        logic_columns = ['g0','g1','g2','g3','g0_b','g1_b','g2_b','g3_b']
    elif model_name == 'TwoGeneFlex':
        logic_columns = ['g0','g1','g2','g3','g0_b','g1_b','g2_b','g3_b']
    else:
        logic_columns = None
        print("Couldn't get model logic columns vals")
    return logic_columns


def calc_sim_prob2d_stats(prob2d):
    # Computing Mutual Information from prob2D
    marg_x = np.sum(prob2d, axis=0)
    marg_y = np.sum(prob2d, axis=1, keepdims=True)
    m_x = np.meshgrid(marg_x, marg_x)[0]
    m_y = np.meshgrid(marg_y, marg_y)[1]
    mutual_info_array = np.multiply(prob2d, (np.log2(prob2d)-np.log2(m_x)-np.log2(m_y)))
    mutual_info=np.sum(mutual_info_array)
    
    # Computing correlation coefficient from prob2D
    x_vals = np.arange(0, prob2d.shape[1], 1)
    y_vals = np.arange(0, prob2d.shape[0], 1)
    prob_x = np.sum(prob2d, axis=0)
    prob_y = np.sum(prob2d, axis=1)
    mean_x = np.sum(np.multiply(prob_x, x_vals))
    mean_y = np.sum(np.multiply(prob_y, y_vals))
    moment_2x = np.sum(np.multiply(prob_x, x_vals**2))
    moment_2y = np.sum(np.multiply(prob_y, y_vals**2))
    variance_x = moment_2x-mean_x**2
    variance_y = moment_2y-mean_y**2
    sigma_x = np.sqrt(variance_x)
    sigma_y = np.sqrt(variance_y)
    mesh_x = np.meshgrid(x_vals, x_vals)[0]
    mesh_y = np.meshgrid(y_vals, y_vals)[1]
    covariance_array = np.multiply(np.multiply((mesh_x-mean_x), (mesh_y-mean_y)), prob2d)
    covariance = np.sum(covariance_array)
    correl_coeff = covariance/sigma_x/sigma_y
    
    # Calculating coexpression
    dropped_origin = prob2d.copy()
    dropped_origin[0, 0]=  0.
    co_exp = prob2d[1:, 1:].sum() / dropped_origin.sum()
    
    return mutual_info, correl_coeff, co_exp


def create_analysis_h5(simdat_hdf5_path, analysis_hdf5_output_path):
    # Openng HDF5 simdat file
    h5_file = h5py.File(simdat_hdf5_path, "r")
    
    # Initializing containers and counters 
    dimensions = np.array(h5_file['RateMatrix']['set_00001'].attrs['dimensions'])
    model_name = h5_file.attrs['model_name']
    total_sets = len(h5_file['RateMatrix'].keys())
    phenotype_count = dimensions[0:2].prod()
    microstate_count = dimensions.prod()
    probvec_array = np.zeros((total_sets, microstate_count))
    prob_2d_vector_array = np.empty((total_sets, phenotype_count))
    mutual_info_array = np.zeros(total_sets)
    correl_coeff_array = np.zeros(total_sets)
    co_exp_array = np.zeros(total_sets)

    # Looping over parameter sets
    for i in range(total_sets):
        set_formatted = 'set_{:05}'.format(i+1)    
        # Loading probvec and prob2d
        probvec_array[i] = h5_file['ProbVec/'+set_formatted].value
        prob_2d = h5_file['Prob2D/'+set_formatted].value
        prob_2d_vector_array[i] = prob_2d.reshape(phenotype_count)
        
        # Computing Mutual Information from prob2D
        mutual_info_array[i], correl_coeff_array[i], co_exp_array[i] = calc_sim_prob2d_stats(prob_2d)

    # Closing simdat file
    h5_file.close()

    # Opening parametes dataframe
    parameters_df = pd.read_csv(trial_path+'/paramValues.csv', index_col=0)

    # Cleaning probability data and calculating entropy's
    probvec_array_cleaned = np.abs(probvec_array + 1e-30)
    system_probvec_entropies = np.real(-np.sum(np.multiply(probvec_array_cleaned, np.log(probvec_array_cleaned)), axis=1))
    prob_2d_vector_array_cleaned = np.abs(prob_2d_vector_array)
    system_prob2d_entropies = np.real(-np.sum(np.multiply(prob_2d_vector_array_cleaned, np.log(prob_2d_vector_array_cleaned)), axis=1))
    relative_binding_strength = np.log10((parameters_df.fa*parameters_df.hr)/(parameters_df.fr*parameters_df.ha))
    X_a = np.log10((parameters_df.fa)/(parameters_df.ha))
    X_r = np.log10((parameters_df.fr)/(parameters_df.hr))

    # Appending columns to parameter dataframe with new parameter set metrics
    parameters_df['model_name'] = model_name
    parameters_df['probvec_entropy'] = np.nan_to_num(system_probvec_entropies)
    parameters_df['prob2d_entropy'] = np.nan_to_num(system_prob2d_entropies)
    parameters_df['mutual_information'] = np.nan_to_num(mutual_info_array)
    parameters_df['correl_coefficient'] = np.nan_to_num(correl_coeff_array)
    parameters_df['coexpression'] = np.nan_to_num(co_exp_array)
    parameters_df['relative_binding_strength'] = np.nan_to_num(relative_binding_strength)
    parameters_df['X_a'] = np.nan_to_num(X_a)
    parameters_df['X_r'] = np.nan_to_num(X_r)
    parameters_df['prob2d_vector'] = [row for row in prob_2d_vector_array]
    parameters_df['probvec_error']= np.array([(x < -1).any() for x in probvec_array])
    parameters_df['prob2d_error'] = np.array([(x < -1).any() for x in prob_2d_vector_array])
    
    parameters_df['errored'] = parameters_df[['probvec_error', 'prob2d_error', 'prob2d_error']].apply(lambda x: True if x.values.any() else False, axis=1)
    
    # Getting logic columns and assigning logic_id
    logic_columns = get_logic_parameter_columns(model_name)
    if logic_columns is not None:
        logic_vals = parameters_df[logic_columns]
        unique_logic_vals = np.unique(logic_vals, axis=0)
        unique_logic_vals = [tuple(i) for i in unique_logic_vals]

        logic_dict = {}
        for i, val in enumerate(unique_logic_vals):
            logic_dict[val] =  i+1
        parameters_df['logic_ID'] = parameters_df[logic_columns].apply(lambda x: logic_dict[tuple(x.values)] if tuple(x.values) in logic_dict.keys() else 0, axis=1)
    else:
        parameters_df['logic_ID'] = 0

    # Saving parameters_df to analysis h5 file
    parameters_df.to_hdf(analysis_hdf5_output_path, '/param_set_outputs_df', mode='w')
    return parameters_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="trial_folder",
                        help="Path to trial folder where simulation_data.h5 file is stored.", required=True)

    # Load input
    args = parser.parse_args()
    trial_path = os.path.join(args.trial_folder)

    # Simulation data and target trial folders
    simdat_hdf5_path = os.path.join(trial_path, 'simulation_data.h5')
    analysis_hdf5_output_path = os.path.join(trial_path, 'Analysis', 'analysis_output.h5')

    # Calling analysis creation function
    df = create_analysis_h5(simdat_hdf5_path, analysis_hdf5_output_path)

