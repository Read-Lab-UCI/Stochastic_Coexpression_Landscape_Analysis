import numpy as np
from .common_calcs import rate_calc

def MISA_Ex_Rxn(parameters):
    # Parameters and model name
    model = 'MISAEx'
    paramSetNum = parameters.Index
    N = parameters.N
    g0 = parameters.g0
    g1 = parameters.g1
    kd = parameters.kd
    ha = parameters.ha
    hr = parameters.hr
    fa = parameters.fa
    fr = parameters.fr
    model_name = model + '_N{}'.format(N)

    # Copy number lists
    A = list(range(N+1))
    B = list(range(N+1))

    # Defining Gene States, Microstates, and their limits
    GeneA_00, GeneA_01, GeneA_10 = [0,1], [0,1], [0,1]
    GeneB_00, GeneB_01, GeneB_10 = [0,1], [0,1], [0,1]
    NumStates = len(A) * len(B) * 3 * 3 # New name for NS
    Smalls = np.array([A[0], B[0], GeneA_00[0], GeneA_01[0], GeneA_10[0], GeneB_00[0], GeneB_01[0], GeneB_10[0]], dtype=int)
    Bigs = np.array([A[-1], B[-1], GeneA_00[-1], GeneA_01[-1], GeneA_10[-1], GeneB_00[-1], GeneB_01[-1], GeneB_10[-1]], dtype=int)

    # Initializing Rxn dict to hold Reactions, Species and Parameters
    Rxn = {}
    NumRxn = 16
    NumSpec = 8

    Rxn['Parameters'] = [g0,g1,g0,ha,hr,fa,fr,kd,g0,g1,g0,ha,hr,fa,fr,kd]
    Rxn['Law'] = np.zeros((NumRxn, NumSpec), dtype=int)
    Rxn['Stoich'] = np.zeros((NumRxn, NumSpec), dtype=int)

    # Reaction Rate Laws, number of each species involved in the reaction
    Rxn['Law'][0,2]=1
    Rxn['Law'][1,3]=1
    Rxn['Law'][2,4]=1
    Rxn['Law'][3,0]=2
    Rxn['Law'][3,2]=1
    Rxn['Law'][4,1]=2
    Rxn['Law'][4,2]=1
    Rxn['Law'][5,3]=1
    Rxn['Law'][6,4]=1
    Rxn['Law'][7,0]=1

    Rxn['Law'][8,5]=1
    Rxn['Law'][9,6]=1
    Rxn['Law'][10,7]=1
    Rxn['Law'][11,1]=2
    Rxn['Law'][11,5]=1
    Rxn['Law'][12,0]=2
    Rxn['Law'][12,5]=1
    Rxn['Law'][13,6]=1
    Rxn['Law'][14,7]=1
    Rxn['Law'][15,1]=1

    # Reaction Stoichiometry, change in species resulting from reaction
    Rxn['Stoich'][0,0] = 1
    Rxn['Stoich'][1,0] = 1
    Rxn['Stoich'][2,0] = 1
    Rxn['Stoich'][3,0] = -2
    Rxn['Stoich'][4,1] = -2
    Rxn['Stoich'][3,2] = -1
    Rxn['Stoich'][3,3] = 1
    Rxn['Stoich'][4,2] = -1
    Rxn['Stoich'][4,4] = 1
    Rxn['Stoich'][5,0] = 2
    Rxn['Stoich'][6,1] = 2
    Rxn['Stoich'][5,3] = -1
    Rxn['Stoich'][5,2] = 1
    Rxn['Stoich'][6,4] = -1
    Rxn['Stoich'][6,2] = 1
    Rxn['Stoich'][7,0] = -1

    Rxn['Stoich'][8,1] = 1
    Rxn['Stoich'][9,1] = 1
    Rxn['Stoich'][10,1] = 1
    Rxn['Stoich'][11,1] = -2
    Rxn['Stoich'][12,0] = -2
    Rxn['Stoich'][11,5] = -1
    Rxn['Stoich'][11,6] = 1
    Rxn['Stoich'][12,5] = -1
    Rxn['Stoich'][12,7] = 1
    Rxn['Stoich'][13,1] = 2
    Rxn['Stoich'][14,0] = 2
    Rxn['Stoich'][13,6] = -1
    Rxn['Stoich'][13,5] = 1
    Rxn['Stoich'][14,7] = -1
    Rxn['Stoich'][14,5] = 1
    Rxn['Stoich'][15,1] = -1

    GeneA_States=[[1,0,0],
                  [0,1,0],
                  [0,0,1]];
    GeneB_States=[[1,0,0],
                  [0,1,0],
                  [0,0,1]];

    return Rxn, A, B, Smalls, Bigs, NumStates, NumSpec, NumRxn, GeneA_States, GeneB_States


def Calc_RateMatrix( Rxn, A, B, Smalls, Bigs, NumStates, NumSpec, NumRxn, GeneA_States, GeneB_States):
    Dimensions = [len(A), len(B), 3, 3]
    StatesList = np.zeros((NumStates, NumSpec))
    RateMatrix = np.zeros((NumStates, NumStates))
    
    TestinRange = [0 for i in range(NumSpec*2)]
    Cur = np.zeros((NumSpec,), dtype=int)
    for i in range(len(A)):
        for j in range(len(B)):
            for k in range(len(GeneA_States)):
                for l in range(len(GeneB_States)):
                    Cur[:] = A[i], B[j], *GeneA_States[k], *GeneB_States[l]
                    CurInd = np.ravel_multi_index((i, j, k, l), Dimensions, order='F')
                    StatesList[CurInd,:] = Cur
                    rate_total = 0.0
                    for m in range(NumRxn):
                        TestDest = Cur + Rxn['Stoich'][m,:]
                        TestinRange[0:NumSpec] = np.greater_equal(TestDest, Smalls)
                        TestinRange[NumSpec:] = np.less_equal(TestDest, Bigs)
                        if all(TestinRange):
                            GeneA_Dest = TestDest[2:5]
                            GeneB_Dest = TestDest[5:8]
                            GAind = GeneA_Dest.nonzero()[0][0]
                            GBind = GeneB_Dest.nonzero()[0][0]
                            Aind = TestDest[0]
                            Bind = TestDest[1]
                            DestInd = np.ravel_multi_index((Aind, Bind, GAind, GBind), Dimensions, order='F')
                            par = Rxn['Parameters'][m]
                            law = Rxn['Law'][m,:]
                            rate_total = rate_total + rate_calc(Cur, law, par)
                    RateMatrix[DestInd, CurInd] = rate_total
    RateMatrix = RateMatrix - np.diagflat(RateMatrix.sum(axis=0))
    return RateMatrix, Dimensions


def main(inputs):
    Rxn, A, B, Smalls, Bigs, NumStates, NumSpec, NumRxn, GeneA_States, GeneB_States = MISA_Ex_Rxn(inputs)
    RateMatrix, Dimensions = Calc_RateMatrix( Rxn, A, B, Smalls, Bigs, NumStates, NumSpec, NumRxn, GeneA_States, GeneB_States)
    return RateMatrix, Dimensions
