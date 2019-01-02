import os
import argparse
import pandas as pd
import scipy.io as sio
from scipy.sparse import lil_matrix

from models import Compute_RateMatrix_MISAEx
from RateMatrix_Calcs import calc_probvec_prob2d


if __name__ == "__main__":
    # Declaring argument parsers and flags
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--parameter", dest="parameterFile", help="Parameter csv file that contains the parameters for the simulations.", required=True)
    parser.add_argument("-o", "--output", dest="outputPath", help="Folder path to output simulation results to.", required=True)

    # Load inputs
    args = parser.parse_args()
    parametersDF = pd.read_csv(args.parameterFile, index_col=0)
    outputPath = os.path.abspath(args.outputPath)
    
    rateMatrixPath = os.path.join(outputPath, 'RateMatrix_py.mat')
    ProbsPath = os.path.join(outputPath, 'Probs_py.mat')
    
    for row in parametersDF.itertuples():
        print(row)
    
        # Calculating and saving rate matrix
        RateMatrix, Dimensions = Compute_RateMatrix_MISAEx.main(row)
        RateMatrix = lil_matrix(RateMatrix)
        sio.savemat(rateMatrixPath, {'RateMatrix': RateMatrix, 'Dimensions': Dimensions})
    
        # Calculating and saving ProbVec, Prob2D and eigenvalues
        eigenValues, probVec, prob2D = calc_probvec_prob2d(RateMatrix, Dimensions)
        sio.savemat(ProbsPath, {'EigenValues': eigenValues, 'ProbVec': probVec, 'Prob2d': prob2D})

