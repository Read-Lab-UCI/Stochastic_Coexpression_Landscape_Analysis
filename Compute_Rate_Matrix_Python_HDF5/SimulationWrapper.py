import os
import argparse
import importlib
import h5py
import multiprocessing as mp
import numpy as np
import pandas as pd
import scipy.io as sio

from RateMatrix_Calcs import basic_calcs


def _workers_count():
    cpu_count = 1
    try:
        cpu_count = len(os.sched_getaffinity(0))
    except AttributeError:
        cpu_count = os.cpu_count()
    return cpu_count


def wrapper(row):
    statesDictPath = os.path.join(outputPath, 'StatesDict.npy')
    set_formatted = 'set_{:05}'.format(row[0])
    
    # Calculating
    RateMatrix, Dimensions, StatesDict = modelFunc.main(row)
    eigenValues, probVec, prob2D, timeScales = basic_calcs(RateMatrix, Dimensions)

    sim_results = (outputPath, set_formatted, RateMatrix, Dimensions, eigenValues, probVec, prob2D, timeScales)

    # Saving
    #sio.savemat(os.path.join(rateMatrixPath, saveFileName), {'RateMatrix': RateMatrix, 'Dimensions': Dimensions}, do_compression = True)
    #sio.savemat(os.path.join(eigenValuesPath, saveFileName), {'EigenValues': eigenValues})
    #sio.savemat(os.path.join(probVecPath, saveFileName), {'ProbVec': probVec})
    #sio.savemat(os.path.join(prob2DPath, saveFileName), {'Prob2d': prob2D})
    #sio.savemat(os.path.join(timeScalesPath, saveFileName), {'TimeScales': timeScales})

    if not os.path.isfile(statesDictPath):
        np.save(os.path.join(outputPath, 'StatesDict.npy'), StatesDict)
    return sim_results


def output_handler(sim_results):
    h5_save_path = os.path.join(sim_results[0][0], 'simulation_data.h5')
    f = h5py.File(h5_save_path, "a")
    f.create_dataset('ProbVec/'+sim_results[0][1], data=sim_results[0][5])
    f.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--parameter", dest="parameterFile", help="Parameter csv file that contains the parameters for the simulations.", required=True)
    parser.add_argument("-o", "--output", dest="outputPath", help="Folder path to output simulation results to.", required=True)
    parser.add_argument("-m", "--model", dest="modelFile", help="Model file used to determine and calculate the RateMatrix, placed inside models/ folder.", required=True)
    parser.add_argument("-pe", "--parallel", dest="runParallel", action='store_true', help="When -pe called, simulations will attempt to be run in parallel")

    # Load inputs
    args = parser.parse_args()
    modelFile = "models." + args.modelFile
    parametersDF = pd.read_csv(args.parameterFile, index_col=0)
    outputPath = os.path.abspath(args.outputPath)
    parametersDF['OutputPath'] = outputPath

    global modelFunc
    modelFunc = importlib.import_module(modelFile)

    if args.runParallel:
        # Creating blank hdf5 file
        savePathHDF5 = os.path.join(outputPath, 'simulation_data.h5')
        f = h5py.File(savePathHDF5, "w")
        f.close()

        # Getting worker count and initiating pool
        num_workers = _workers_count()
        print('Starting pool with', num_workers, 'workers')
        pool = mp.Pool(processes=num_workers)

        # Distributing tasks and executing
        for row in parametersDF.itertuples(name=None):
            pool.map_async(wrapper, (row, ), callback=output_handler)
        pool.close()
        pool.join()

    else:
        print('Running simulations sequentially')
        for row in parametersDF.itertuples(name=None):
            wrapper(row)
