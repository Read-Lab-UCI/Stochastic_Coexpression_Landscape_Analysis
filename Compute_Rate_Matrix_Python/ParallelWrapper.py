import os
import argparse
import multiprocessing as mp
import pandas as pd
import scipy.io as sio

from models import Compute_RateMatrix_MISAEx_parallel
from RateMatrix_Calcs import calc_probvec_prob2d


def _workers_count():
    cpu_count = 0
    try:
        cpu_count = len(os.sched_getaffinity(0))
    except AttributeError:
        cpu_count = os.cpu_count()
    return cpu_count


def parallel_wrapper(row):
    # Building Paths
    rateMatrixPath = os.path.join(row[-1], 'RateMatrix')
    eigenValuesPath =  os.path.join(row[-1], 'EigenValues')
    probVecPath =  os.path.join(row[-1], 'ProbVec')
    prob2DPath = os.path.join(row[-1], 'Prob2D')
    timeScalesPath = os.path.join(row[-1], 'TimeScales')
    saveFileName = 'set_{:05}.mat'.format(row[0])
    
    # Calculating
    RateMatrix, Dimensions, StatesDict = Compute_RateMatrix_MISAEx_parallel.main(row)
    eigenValues, probVec, prob2D, timeScales = calc_probvec_prob2d(RateMatrix, Dimensions)
    # Saving
    sio.savemat(os.path.join(rateMatrixPath, saveFileName), {'RateMatrix': RateMatrix, 'Dimensions': Dimensions})
    sio.savemat(os.path.join(eigenValuesPath, saveFileName), {'EigenValues': eigenValues})
    sio.savemat(os.path.join(probVecPath, saveFileName), {'ProbVec': probVec})
    sio.savemat(os.path.join(prob2DPath, saveFileName), {'Prob2d': prob2D})
    sio.savemat(os.path.join(timeScalesPath, saveFileName), {'TimeScales': timeScales})


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--parameter", dest="parameterFile", help="Parameter csv file that contains the parameters for the simulations.", required=True)
    parser.add_argument("-o", "--output", dest="outputPath", help="Folder path to output simulation results to.", required=True)

    # Load inputs
    args = parser.parse_args()
    parametersDF = pd.read_csv(args.parameterFile, index_col=0)
    outputPath = os.path.abspath(args.outputPath)
    parametersDF['OutputPath'] = outputPath
 
    tups = parametersDF.itertuples(name=None)
    
    # Workers
    num_workers = _workers_count()
    print(num_workers)
    # Initiating pool
    print('Starting pool')
    pool = mp.Pool(processes=num_workers)
    results = pool.map_async(parallel_wrapper, tups).get()
