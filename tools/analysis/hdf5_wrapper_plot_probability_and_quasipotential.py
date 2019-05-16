import os
import io
import h5py
import multiprocessing as mp
import argparse
import numpy as np
import pandas as pd
import pathlib
from PIL import Image

from module_tools import generate_2d_heatmap
old_settings = np.seterr(all='ignore')

def _workers_count():
    cpu_count = 1
    try:
        cpu_count = len(os.sched_getaffinity(0))
    except AttributeError:
        cpu_count = os.cpu_count()
    return cpu_count


def output_handler(analysis_results):
    analysis_results = analysis_results[0]
    outputFile = analysis_results['outputFile']
    setName = analysis_results['set']
    analysis_h5_file = h5py.File(outputFile, "a")

    for key in analysis_results['buffers'].keys():
        # Selecting buffer and setting to first position
        buf = analysis_results['buffers'][key]
        buf.seek(0)

        # Reading binary as an image object
        img = Image.open(buf)
        img = np.asarray((img), dtype="uint8")
        
        # Saving image and HDFView attributes
        dset = analysis_h5_file.create_dataset(key+'/'+setName, data=img, shape=img.shape, maxshape=img.shape, dtype='uint8', compression="gzip")
        dset.attrs['CLASS'] = np.string_('IMAGE')
        dset.attrs['IMAGE_VERSION'] = np.string_('1.2')
        dset.attrs['IMAGE_MINMAXRANGE'] = np.asarray([0, 255], dtype='uint8')
        dset.attrs['IMAGE_SUBCLASS'] = np.string_('IMAGE_TRUECOLOR')
        dset.attrs['INTERLACE_MODE'] = np.string_('INTERLACE_PIXEL')

    analysis_h5_file.close()


def wrapper(args):
    # Assigning and calculating landscapes
    prob2d = args['prob2d']
    dimensions = args['dimensions']
    setName = args['set']
    outputFile = args['outputFile']
    q_potential = -np.log(np.abs(prob2d))

    # Generating plots and returning buffers
    axisLabels = ['D1 States', 'D2 States']
    prob_buffer = generate_2d_heatmap.plot_2d_heatmap_buffer(prob2d, title='Probability', axis_labels=axisLabels, save_fig=True)
    q_potent_buffer = generate_2d_heatmap.plot_2d_heatmap_buffer(q_potential, title='Quasipotential', axis_labels=axisLabels, save_fig=True, invert_color = True)

    # Returning packed values
    return {'buffers': {'Probability': prob_buffer, 'Quasipotential': q_potent_buffer}, 'outputFile': outputFile, 'set': setName}
    

if __name__ == "__main__":
    # Declaring argument parsers and flags
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", dest="trialPath", help="Path to folder containing the simulation trial data", required=True)
    parser.add_argument("-pe", "--parallel", dest="runParallel", action='store_true',
                        help="When -pe called, simulations will attempt to be run in parallel")
    
    # Load input file
    args = parser.parse_args()
    trialPath = os.path.abspath(args.trialPath)
    trialDataPath = os.path.join(trialPath, 'simulation_data.h5')

    # Creating blank hdf5 file
    savePathHDF5 = os.path.join(trialPath, 'Analysis', 'landscapes.h5')
    f = h5py.File(savePathHDF5, "w")
    f.close()

    # Getting simulation set names
    trialDataFile = h5py.File(trialDataPath, "r")
    trialSetNames = list(trialDataFile['RateMatrix'].keys())

    if args.runParallel:
        # Getting worker count and initiating pool
        num_workers = _workers_count()
        print('Starting pool with', num_workers, 'workers')
        pool = mp.Pool(processes=num_workers)

        
        # Distributing tasks 
        for setName in trialSetNames:
            prob2d = np.array(trialDataFile['Prob2D/'+setName])
            dimensions = np.array(trialDataFile['RateMatrix'][setName].attrs['dimensions'])
            argDict = {'prob2d': prob2d, 'dimensions': dimensions, 'set': setName, 'outputFile': savePathHDF5}
            pool.map_async(wrapper, (argDict,), callback=output_handler)
        
        # Closing simulation data file
        trialDataFile.close()
        
        # Executing pool tasks
        pool.close()
        pool.join()

    else:
        print('Running simulations sequentially')
        for setName in trialSetNames:
            prob2d = np.array(trialDataFile['Prob2D/'+setName])
            dimensions = np.array(trialDataFile['RateMatrix'][setName].attrs['dimensions'])
            argDict = {'prob2d': prob2d, 'dimensions': dimensions, 'set': setName, 'outputFile': savePathHDF5}
            results = wrapper(argDict)
            output_handler(results)

        # Closing simulation data file
        trialDataFile.close()

