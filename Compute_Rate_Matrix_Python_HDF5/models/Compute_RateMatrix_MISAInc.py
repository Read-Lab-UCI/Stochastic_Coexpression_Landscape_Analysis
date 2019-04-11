import numpy as np
from scipy.sparse import coo_matrix
from .common_calcs import Update_RateMatrix


def MISA_Inc_Rxn(parameters):
    # Parameters and model name
    model = 'MISAInc'
    N = parameters[1]
    g0 = parameters[2]
    g1 = parameters[3]
    kd = parameters[4]
    ha = parameters[5]
    hr = parameters[6]
    fa = parameters[7]
    fr = parameters[8]
    model_name = model + '_N{}'.format(N)

    # Copy number lists
    A = list(range(N+1))
    B = list(range(N+1))

    # Defining Gene States
    GeneA_00, GeneA_01, GeneA_10, GeneA_11 = [0,1], [0,1], [0,1], [0,1]
    GeneB_00, GeneB_01, GeneB_10, GeneB_11 = [0,1], [0,1], [0,1], [0,1]
    GeneA_States=[[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]]
    GeneB_States=[[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]]

    # Number of Microstate = CopyA * CopyB * GeneStatesA * GeneStatesB
    NumStates = (N+1) * (N+1) * 4 * 4

    # Initializing Rxn dict to hold Reactions, Species and Parameters
    Rxn = {}
    NumRxn = 26
    NumSpec = 10

    Rxn['Parameters'] = np.array([g0,g1,g0,g0,ha,hr,hr,ha,fa,fr,fa,fr,kd,g0,g1,g0,g0,ha,hr,hr,ha,fa,fr,fa,fr,kd])
    Rxn['Law'] = np.zeros((NumRxn, NumSpec), dtype=int)
    Rxn['Stoich'] = np.zeros((NumRxn, NumSpec), dtype=int)

    # Reaction Rate Laws, number of each species involved in the reaction
    Rxn['Law'][0,2]=1
    Rxn['Law'][1,3]=1
    Rxn['Law'][2,4]=1
    Rxn['Law'][3,5]=1
    Rxn['Law'][4,0]=2 
    Rxn['Law'][4,2]=1
    Rxn['Law'][5,1]=2 
    Rxn['Law'][5,2]=1
    Rxn['Law'][6,1]=2 
    Rxn['Law'][6,3]=1
    Rxn['Law'][7,0]=2 
    Rxn['Law'][7,4]=1
    Rxn['Law'][8,3]=1
    Rxn['Law'][9,4]=1
    Rxn['Law'][10,5]=1
    Rxn['Law'][11,5]=1
    Rxn['Law'][12,0]=1
    
    Rxn['Law'][13,6]=1
    Rxn['Law'][14,7]=1
    Rxn['Law'][15,8]=1
    Rxn['Law'][16,9]=1
    Rxn['Law'][17,1]=2 
    Rxn['Law'][17,6]=1
    Rxn['Law'][18,0]=2 
    Rxn['Law'][18,6]=1
    Rxn['Law'][19,0]=2 
    Rxn['Law'][19,7]=1
    Rxn['Law'][20,1]=2 
    Rxn['Law'][20,8]=1
    Rxn['Law'][21,7]=1
    Rxn['Law'][22,8]=1
    Rxn['Law'][23,9]=1
    Rxn['Law'][24,9]=1
    Rxn['Law'][25,1]=1

    # Reaction Stoichiometry, change in species resulting from reaction
    Rxn['Stoich'][0,0] = 1
    Rxn['Stoich'][1,0] = 1
    Rxn['Stoich'][2,0] = 1
    Rxn['Stoich'][3,0] = 1
    Rxn['Stoich'][4,2] = -1 
    Rxn['Stoich'][4,0] = -2 
    Rxn['Stoich'][4,3] = 1
    Rxn['Stoich'][5,2] = -1 
    Rxn['Stoich'][5,1] = -2 
    Rxn['Stoich'][5,4] = 1
    Rxn['Stoich'][6,3] = -1 
    Rxn['Stoich'][6,1] = -2 
    Rxn['Stoich'][6,5] = 1
    Rxn['Stoich'][7,4] = -1 
    Rxn['Stoich'][7,0] = -2 
    Rxn['Stoich'][7,5] = 1
    Rxn['Stoich'][8,3] = -1 
    Rxn['Stoich'][8,0] = 2 
    Rxn['Stoich'][8,2] = 1
    Rxn['Stoich'][9,4] = -1 
    Rxn['Stoich'][9,1] = 2 
    Rxn['Stoich'][9,2] = 1
    Rxn['Stoich'][10,5] = -1 
    Rxn['Stoich'][10,0] = 2 
    Rxn['Stoich'][10,4] = 1
    Rxn['Stoich'][11,5] = -1 
    Rxn['Stoich'][11,1] = 2 
    Rxn['Stoich'][11,3] = 1
    Rxn['Stoich'][12,0] = -1
    
    
    Rxn['Stoich'][13,1] = 1
    Rxn['Stoich'][14,1] = 1
    Rxn['Stoich'][15,1] = 1
    Rxn['Stoich'][16,1] = 1
    Rxn['Stoich'][17,6] = -1 
    Rxn['Stoich'][17,1] = -2 
    Rxn['Stoich'][17,7] = 1
    Rxn['Stoich'][18,6] = -1 
    Rxn['Stoich'][18,0] = -2 
    Rxn['Stoich'][18,8] = 1
    Rxn['Stoich'][19,7] = -1 
    Rxn['Stoich'][19,0] = -2 
    Rxn['Stoich'][19,9] = 1
    Rxn['Stoich'][20,8] = -1
    Rxn['Stoich'][20,1] = -2
    Rxn['Stoich'][20,9] = 1
    Rxn['Stoich'][21,7] = -1
    Rxn['Stoich'][21,1] = 2
    Rxn['Stoich'][21,6] = 1
    Rxn['Stoich'][22,8] = -1
    Rxn['Stoich'][22,0] = 2
    Rxn['Stoich'][22,6] = 1
    Rxn['Stoich'][23,9] = -1
    Rxn['Stoich'][23,1] = 2
    Rxn['Stoich'][23,8] = 1
    Rxn['Stoich'][24,9] = -1
    Rxn['Stoich'][24,0] = 2
    Rxn['Stoich'][24,7] = 1
    Rxn['Stoich'][25,1] = -1

    return Rxn, A, B, NumStates, NumSpec, NumRxn, GeneA_States, GeneB_States


def Determine_StatesDict(Dimensions, GeneA_States, GeneB_States):
    # original stateslist creation, returns ordered numpy array
    States_dict = {}

    for i in range(Dimensions[0]):
        for j in range(Dimensions[1]):
            for k in range(Dimensions[2]):
                for l in range(Dimensions[3]):
                    Cur = (i, j, *GeneA_States[k], *GeneB_States[l])
                    CurInd = np.ravel_multi_index((i, j, k, l), Dimensions, order='F')
                    States_dict[Cur] = CurInd
    return States_dict


def Calc_RateMatrix( Rxn, StatesDict, NumStates, NumRxn):
    StatesKeys = set(StatesDict.keys())
    MaxNumInteractions = NumStates * NumRxn
    CurInd_vals = np.zeros((MaxNumInteractions))
    DestInd_vals = np.zeros((MaxNumInteractions))
    RateMatrix_vals = np.zeros((MaxNumInteractions))

    n=0
    for state in StatesKeys:
        Cur = np.array(state, dtype=int)
        CurInd = StatesDict[state]
        for m in range(NumRxn):
            TestDest = tuple(Cur + Rxn['Stoich'][m,:])
            if TestDest in StatesKeys:
                n=n+1
                DestInd = StatesDict[TestDest]
                CurInd_vals, DestInd_vals, RateMatrix_vals = Update_RateMatrix(CurInd_vals, DestInd_vals, RateMatrix_vals, Cur, CurInd, DestInd, Rxn['Parameters'], Rxn['Law'], m, n)

    RateMatrix = coo_matrix((RateMatrix_vals, (DestInd_vals, CurInd_vals)), shape=(NumStates, NumStates)).tolil()
    RateMatrix.setdiag((RateMatrix.diagonal() - RateMatrix.sum(axis=0)).A[0])
    RateMatrix = RateMatrix.tocsc()

    return RateMatrix

def main(inputs):
    Rxn, A, B, NumStates, NumSpec, NumRxn, GeneA_States, GeneB_States = MISA_Inc_Rxn(inputs)
    Dimensions = [len(A), len(B), len(GeneA_States), len(GeneB_States)]
    StatesDict = Determine_StatesDict(Dimensions, GeneA_States, GeneB_States)
    RateMatrix = Calc_RateMatrix(Rxn, StatesDict, NumStates, NumRxn)

    return RateMatrix, Dimensions, StatesDict
