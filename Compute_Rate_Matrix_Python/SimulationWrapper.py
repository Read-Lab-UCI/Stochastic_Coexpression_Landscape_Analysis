import os
import argparse
import pandas as pd
import scipy.io as sio

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
    
    rateMatrixPath = os.path.join(outputPath, 'RateMatrix')
    eigenValuesPath =  os.path.join(outputPath, 'EigenValues')
    probVecPath =  os.path.join(outputPath, 'ProbVec')
    prob2DPath = os.path.join(outputPath, 'Prob2D')
    timeScalesPath = os.path.join(outputPath, 'TimeScales')
    
    for row in parametersDF.itertuples():
        print(row)
        saveFileName = 'set_{:05}.mat'.format(row.Index)

        # Calculating and saving rate matrix
        RateMatrix, Dimensions, StatesDict = Compute_RateMatrix_MISAEx.main(row)
        sio.savemat(os.path.join(rateMatrixPath, saveFileName), {'RateMatrix': RateMatrix, 'Dimensions': Dimensions})
    
        # Calculating and saving ProbVec, Prob2D and eigenvalues
        eigenValues, probVec, prob2D, timeScales = calc_probvec_prob2d(RateMatrix, Dimensions)
        sio.savemat(os.path.join(eigenValuesPath, saveFileName), {'EigenValues': eigenValues})
        sio.savemat(os.path.join(probVecPath, saveFileName), {'ProbVec': probVec})
        sio.savemat(os.path.join(prob2DPath, saveFileName), {'Prob2d': prob2D})
        sio.savemat(os.path.join(timeScalesPath, saveFileName), {'TimeScales': timeScales})

    np.save(os.path.join(outputPath, 'StatesDict.npy'), StatesDict)
