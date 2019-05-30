import os
import logging
import argparse
import importlib
import h5py
import multiprocessing as mp
import numpy as np
import pandas as pd

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

    if not os.path.isfile(statesDictPath):
        np.save(os.path.join(outputPath, 'StatesDict.npy'), StatesDict)
    return sim_results


def output_handler(sim_results):
    h5_save_path = os.path.join(sim_results[0][0], 'simulation_data.h5')
    h5_file = h5py.File(h5_save_path, "a")

    for sim_result in sim_results:
        # Breaking down Rate Matrix to allow for HDF5 saving
        # See https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html for reconstruction info
#        import scipy.io as sio
#        sio.savemat('ratematrix_py.mat', {'RateMatrix':sim_result[2]})
        g = h5_file.create_group('RateMatrix/' + sim_result[1])
        g.create_dataset('data', data=sim_result[2].data, compression='gzip')
        g.create_dataset('indptr', data=sim_result[2].indptr, compression='gzip')
        g.create_dataset('indices', data=sim_result[2].indices, compression='gzip')
        g.attrs['shape'] = sim_result[2].shape
        g.attrs['dimensions'] = sim_result[3]

        # Saving other outputs
        h5_file.create_dataset('EigenValues/' + sim_result[1], data=sim_result[4])
        h5_file.create_dataset('ProbVec/' + sim_result[1], data=sim_result[5])
        h5_file.create_dataset('Prob2D/' + sim_result[1], data=sim_result[6])
        h5_file.create_dataset('TimeScales/' + sim_result[1], data=sim_result[7])
    h5_file.close()
    logging.info('Finished '+sim_results[-1][1])


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--parameter", dest="parameterFile",
                        help="Parameter csv file that contains the parameters for the simulations.", required=True)
    parser.add_argument("-o", "--output", dest="outputPath",
                        help="Folder path to output simulation results to.", required=True)
    parser.add_argument("-m", "--model", dest="modelFile",
                        help="Model file used to calculate the RateMatrix, found in models/ .", required=True)
    parser.add_argument("-pe", "--parallel", dest="runParallel", action='store_true',
                        help="When -pe called, simulations will attempt to be run in parallel")

    # Load inputs
    args = parser.parse_args()
    modelFile = "models." + args.modelFile
    parametersDF = pd.read_csv(args.parameterFile, index_col=0)
    outputPath = os.path.abspath(args.outputPath)
    parametersDF['OutputPath'] = outputPath

    # Setting up logger
    outputLogger = logging.getLogger()
    outputLogger.setLevel(logging.DEBUG)
    logOutputFile = os.path.join(outputPath, 'simulations.log')
    formatter = logging.Formatter('%(message)s')
    fh = logging.FileHandler(logOutputFile, 'w')
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    outputLogger.addHandler(fh)

    # Number of simulations each worker will run
    chunk_size = 4

    global modelFunc
    modelFunc = importlib.import_module(modelFile)

    # Creating blank hdf5 file
    savePathHDF5 = os.path.join(outputPath, 'simulation_data.h5')
    f = h5py.File(savePathHDF5, "w")
    f.attrs['model_name'] = modelFunc.model_name()
    f.close()

    if args.runParallel:
        # Getting worker count and initiating pool
        num_workers = _workers_count()
        print('Starting pool with', num_workers, 'workers')
        pool = mp.Pool(processes=num_workers)

        # Distributing tasks and executing
        for i in range(0, len(parametersDF), chunk_size):
            slc = parametersDF.iloc[i: i+chunk_size]
            pool.map_async(wrapper, slc.itertuples(name=None), callback=output_handler)
        pool.close()
        pool.join()

    else:
        print('Running simulations sequentially')
        for i in range(0, len(parametersDF), chunk_size):
            slc = parametersDF.iloc[i: i+chunk_size]
            sequential_sim_results = [wrapper(row) for row in slc.itertuples(name=None)]
            output_handler(sequential_sim_results)
