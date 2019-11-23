import math
import sys
import itertools
import numpy as np
import pandas as pd
import scipy.io as sio

# This script will determine the combinations of parameters

# Values
allParameters = ["N", "kd", "g0", "g1", "g2", "g3", "g0_b", "g1_b", "g2_b", "g3_b", "ha", "hr", "fa", "fr"]
parametersToChange = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13] # Indexes of the parameters to be varied in allParameters
constants = {"N":25,
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
                                 [1e-3, 4.0, 1e-3, 1e-3, 1e-3, 7.0, 1e-3, 1e-3],
                                 [1e-3, 7.0, 1e-3, 1e-3, 1e-3, 4.0, 1e-3, 1e-3],
                                 [1e-3, 7.0, 1e-3, 1e-3, 1e-3, 7.0, 1e-3, 1e-3]])

    # Creating binding variation portion
    expon_h = np.linspace(0, 2, 9)
    expon_f = np.linspace(0, 5, 9)
    print('expon_h:', expon_h)
    print('expon_f:', expon_f)

    h_vals = np.power(10, expon_h)
    f_vals = np.power(10, expon_f)
    print('h_vals:', h_vals)
    print('f_vals:', f_vals)

    h_vals_combinations = list(itertools.product(h_vals, repeat=2))
    f_vals_combinations = list(itertools.product(f_vals, repeat=2))

    variedBindingMatrix = list(itertools.product(h_vals_combinations, f_vals_combinations))
    variedBindingMatrix = [i[0] + i[1] for i in variedBindingMatrix]


    paramMatrix_seperated = list(itertools.product(productionMatrix, variedBindingMatrix))
    paramMatrix = [list(itertools.chain.from_iterable(row)) for row in paramMatrix_seperated]

    setCount = len(paramMatrix)
    variedParamCount = len(parametersToChange)
    
    paramDF = pd.DataFrame(0, index=range(1,setCount+1), columns=allParameters)
    paramDF.index.rename("paramSetNum", inplace=True)
    
    for key in constants.keys():
        paramDF[key] = constants[key]

    paramDF[[allParameters[i] for i in parametersToChange]] = paramMatrix
    paramDF.to_csv(csvFileName)
