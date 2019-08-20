import math
import sys
import itertools
import numpy as np
import pandas as pd
import scipy.io as sio

# This script will determine the combinations of parameters

# Values
allParameters = ["N", "kd", "g0", "g1", "g2", "g3", "g0_b", "g1_b", "g2_b", "g3_b", "ha", "hr", "fa", "fr", "c_c", "c_cr", "c_o"]
parametersToChange = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16] # Indexes of the parameters to be varied in allParameters
constants = {"N":30,
             "kd":1}  # The values for parameters not being varied

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print('Missing output file path, this will be placed in this  directory.')
        csvFileName = sys.argv[0]+'.csv'
    else:
        csvFileName=sys.argv[1]
        print('csv file will placed in {}, be sure to change the PARAM_CSV variable in the trial_init.sh  script to refer to this name. By default trial_init looks for paramValues.csv'.format(csvFileName))

    
    # Creating logic variation portion
    productionMatrix = np.array([[1e-3, 4.0, 1e-3, 1e-3, 1e-3, 4.0, 1e-3, 1e-3],
                                 [1e-3, 7.0, 1e-3, 1e-3, 1e-3, 4.0, 1e-3, 1e-3],
                                 [1e-3, 4.0, 1e-3, 1e-3, 1e-3, 7.0, 1e-3, 1e-3],
                                 [1e-3, 7.0, 1e-3, 1e-3, 1e-3, 7.0, 1e-3, 1e-3]])

    # Creating binding variation portion
    bindingRateValues = [1e-2, 1e0, 1e2, 1e4, 1e6]
    bindingMatrix = list(itertools.product(bindingRateValues, repeat=4))
    
    # Creating chromatin variation portion
    chromatinParameters = {"c_c": [1e-5,], "c_cr":[1e-5, 1e-3, 1e-1], "c_o":[1e-5, 1e-3, 1.]}
    chromatinMatrix = list(itertools.product(*chromatinParameters.values()))

    paramMatrix_seperated = list(itertools.product(productionMatrix, bindingMatrix, chromatinMatrix))
    paramMatrix = [list(itertools.chain.from_iterable(row)) for row in paramMatrix_seperated]

    setCount = len(paramMatrix)
    variedParamCount = len(parametersToChange)
    
    paramDF = pd.DataFrame(0, index=range(1,setCount+1), columns=allParameters)
    paramDF.index.rename("paramSetNum", inplace=True)
    
    for key in constants.keys():
        paramDF[key] = constants[key]

    paramDF[[allParameters[i] for i in parametersToChange]] = paramMatrix
    paramDF.to_csv(csvFileName)
