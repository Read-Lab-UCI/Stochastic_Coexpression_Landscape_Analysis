import argparse
import numpy as np
import scipy.io as sio
import pathlib
from module_tools import generate_2d_heatmap
from module_tools import get_file_paths

# TODO: Turn this into documentation string


def calc_prob2d_and_quasipotential(probvec, dimensions, dimensions_to_reduce=(2,3)):
    # Collapse ProbVec to 2 desired dimensions and calculate quasipotential for visualization
    # dimensions_to_reduce are axis corresponding to number of gene states for MISA Ex model, TODO: fix hardcode

    prob_full_d = probvec.reshape(dimensions, order='F')  # Collapses probability along all system dimensions
    prob_2d = np.sum(prob_full_d, axis=dimensions_to_reduce)  # Reduces probability to two dimensions
    q_potential = -np.log(np.abs(prob_2d))  # Calculates quasi-potential
    
    return prob_2d, q_potential


if __name__ == "__main__":
    # Declaring argument parsers and flags
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", dest="trialPath", help="Path to folder containing the simulation trial data", required=True)
    parser.add_argument("-s", "--set", dest="paramSet", help="Parameter set number within trial to be visualized.", required=True, type=int)

    # Load input file
    args = parser.parse_args()
    paramSetNum = args.paramSet
    trialPath = args.trialPath
    paramSetPaths = get_file_paths.generate_set_paths(trialPath, paramSetNum)

    # Assigning Paths
    outputPath = paramSetPaths['analysis_folder']
    probVecPath = paramSetPaths['probvec']
    rateMatrixPath = paramSetPaths['ratematrix']

    # Creating folders
    foldersNeeded = ['Landscapes', 'Quasipotentials']
    for folder in foldersNeeded:
        pathlib.Path(outputPath+folder).mkdir(parents=True, exist_ok=True)

    # Loading mat files
    probVec = sio.loadmat(probVecPath)['ProbVec']
    dimensions = sio.loadmat(rateMatrixPath)['Dimensions'][0]
    
    # Reducing dimensions and calculating potential
    prob2D, quasipotential = calc_prob2d_and_quasipotential(probVec, dimensions)

    # Plotting and saving
    axisLabels = ['D1 States', 'D2 States']
    probabilitySavePath = outputPath + 'Landscapes/set_{:05}.png'.format(paramSetNum)
    quasipotentialSavePath = outputPath + 'Quasipotentials/set_{:05}.png'.format(paramSetNum)

    generate_2d_heatmap.plot_2d_heatmap(prob2D, title='Probability', axis_labels=axisLabels, save_fig=True, fig_name=probabilitySavePath)
    generate_2d_heatmap.plot_2d_heatmap(quasipotential, title='Quasipotential', axis_labels=axisLabels, save_fig=True, fig_name=quasipotentialSavePath, invert_color = True)
        
        
