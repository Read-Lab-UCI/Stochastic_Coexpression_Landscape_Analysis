import math
import sys
import itertools
import numpy as np
import pandas as pd
import scipy.io as sio

# TODO ADD THIS SAVE FUNCTIONALITY TO THE parameterCreationScripts FOLDER AND CALL DURING WRITE
# This script will determine the combinations of parameters

# Values
allParameters = ["N", "kd", "g0", "g1", "g2", "g3", "g0_b", "g1_b", "g2_b", "g3_b", "ha", "hr", "fa", "fr"]
parametersToChange = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13] # Indexes of the parameters to be varied in allParameters
constants = {"N":20,
             "kd":1}  # The values for parameters not being varied


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print('Missing output file path, this will be placed in this  directory.')
        csvFileName = sys.argv[0]+'.csv'
    else:
        csvFileName=sys.argv[1]
        print('csv file will placed in {}, be sure to change the PARAM_CSV variable in the trial_init.sh  script to refer to this name. By default trial_init looks for paramValues.csv'.format(csvFileName))

    binding_values = np.array([1.00000000e-05, 1.00000000e-03, 1.00000000e-01, 1.00000000e+01, 1.00000000e+03])
    
    # Creating logic variation portion
    variedLogicMatrix = sio.loadmat('gmatrix_asym.mat')['gmatrix']

    # Creating logic variation portion
    variedParamMatrix = list(itertools.product(binding_values, repeat=4))
    variedLogicMatrix = [tuple(row) for row in variedLogicMatrix]
    
    setCount = len(variedLogicMatrix) * len(variedParamMatrix)
    variedParamCount = len(parametersToChange)
    paramMatrix = np.zeros((setCount, variedParamCount))

    i = 0
    for logic_set in variedLogicMatrix:
        for param_set in variedParamMatrix:
            paramMatrix[i] = logic_set+param_set
            i = i+1
    
    paramDF = pd.DataFrame(0, index=range(1,setCount+1), columns=allParameters)
    paramDF.index.rename("paramSetNum", inplace=True)
    
    for key in constants.keys():
        paramDF[key] = constants[key]
    
    paramDF[[allParameters[i] for i in parametersToChange]] = paramMatrix
    
    paramDF.to_csv(csvFileName)
