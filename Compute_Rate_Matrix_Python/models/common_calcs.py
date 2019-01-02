from numba.decorators import jit

@jit('f8(i8[:], i8[:], f8)')
def rate_calc(cur, law, rate):
    tmp1 = 0
    tmp2 = 1
    for i in range(cur.shape[0]):
        tmp1 = cur[i]**law[i]
        tmp2 = 1
        for j in range(law[i]):
            tmp2 = tmp2 * (j+1)
        rate = rate * (tmp1/tmp2)
    return rate 
