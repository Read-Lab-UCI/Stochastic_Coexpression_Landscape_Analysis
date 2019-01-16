import numpy as np
from scipy.sparse.linalg import eigs

def basic_calcs(matrix, dimensions, dims_to_reduce=(2,3), nsig=15, sigma=1e-12):

    eigen_values, right_eigen_vectors = eigs(matrix, k=nsig, sigma=sigma)
    prob_vec = right_eigen_vectors[:,1]/np.sum(right_eigen_vectors[:,1])

    prob_full_d = prob_vec.reshape(dimensions, order='F')
    prob_2d = np.sum(prob_full_d, axis=dims_to_reduce)

    time_scales = -np.divide(1, eigen_values[1:].real)

    return eigen_values, prob_vec, prob_2d, time_scales
