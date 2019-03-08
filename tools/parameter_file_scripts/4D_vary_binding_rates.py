import math
import sys
import itertools
import numpy as np
import pandas as pd

# TODO ADD THIS SAVE FUNCTIONALITY TO THE parameterCreationScripts FOLDER AND CALL DURING WRITE
# This script will determine the combinations of parameters

# Values
allParameters = ["N", "g0", "g1", "kd", "ha", "hr", "fa", "fr"]
parametersToChange = [4, 5, 6, 7] # Indexes of the parameters to be varied in allParameters
constants = {"N":20,
             "g0":2,
             "g1":12,
             "kd":1}  # The values for parameters not being varied
binCount = 5  # Number of bins that the parameter space should be divided over
param_min = 1e-5
param_max = 1e5

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print('Missing output file path, this will be placed in this  directory.')
        csvFileName = sys.argv[0]+'.csv'
    else:
        csvFileName=sys.argv[1]
        print('csv file will placed in {}, be sure to change the PARAM_CSV variable in the trial_init.sh  script to refer to this name. By default trial_init looks for paramValues.csv'.format(csvFileName))

    #values = np.power(10, (np.linspace(math.log(param_min, 10), math.log(param_max, 10), binCount)))
    values = np.array([1.00000000e-05, 1.00000000e-03, 1.00000000e+00, 1.00000000e+03, 1.00000000e+05])
    
    # Building variation grid
    #vGrid = []
    #for p in allParameters:
    #    if p in parametersToChange:
    #        vGrid.append(True)
    #    else:
    #        vGrid.append(False)

    setCount = len(list(itertools.product(values, repeat=len(parametersToChange))))
    allParamCount = len(allParameters)
    variedParamCount = len(parametersToChange)
    
    variedParamMatrix = np.zeros((setCount, variedParamCount))
    variedParamMatrix[:] = [pSet for pSet in list(itertools.product(values, repeat=len(parametersToChange)))]
    
    paramDF = pd.DataFrame(0, index=range(1,setCount+1), columns=allParameters)
    paramDF.index.rename("paramSetNum", inplace=True)
    
    for key in constants.keys():
        paramDF[key] = constants[key]
    
    paramDF[[allParameters[i] for i in parametersToChange]] = variedParamMatrix
    
    paramDF.to_csv('../../simulation/parameter_files/{}'.format(csvFileName))
