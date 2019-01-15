from numba.decorators import jit

@jit('Tuple((f8[:], f8[:], f8[:]))(f8[:], f8[:], f8[:], i8[:], i8, i8, f8[:], i8[:,:], i8, i8)')
def Update_RateMatrix(CurInd_vals, DestInd_vals, RateMatrix_vals, Cur, CurInd, DestInd, parameters, laws, m, n):
    par = parameters[m]
    law = laws[m,:]
    for i in range(Cur.shape[0]):
        tmp1 = Cur[i]**law[i]
        tmp2 = 1
        for j in range(law[i]):
            tmp2 = tmp2 * (j+1)
        par = par * (tmp1/tmp2)

    CurInd_vals[n] = CurInd
    DestInd_vals[n] = DestInd
    RateMatrix_vals[n] = par

    return CurInd_vals, DestInd_vals, RateMatrix_vals
