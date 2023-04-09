
import numpy as np


from itertools import product
from scipy.sparse import coo_array
# https://docs.scipy.org/doc/scipy/reference/sparse.html

def build_assembly_coo_array(matrices):
    '''
    Parameters
    ----------
    matrices : iterable
        Iterable containing 2D numpy array
        ex: [K1, K2, ..., KN]
    
    Returns
    -------
    assembly : coo_array
        An array containing the matrices assembly
        ex: [K1, 0, 0
             0,  K2, 0
             0,  0,  KN]
    '''

    rows = 0
    cols = 0
    offset = 0
    idx = []

    for matrix in matrices:
        r, c = matrix.shape
        rows += r
        cols += c
        ri = np.arange(r) + offset
        ci = np.arange(c) + offset

        idx.append(
            np.array(list(product(ri, ci)))
        )
        offset += r

    data = []
    row = []
    col = []

    for i, matrix in enumerate(matrices):
        m_values = list(
            matrix.reshape(matrix.size))
        
        data.extend(m_values)
        row.extend(list(
            idx[i][::, 0])
        )
        col.extend(list(
            idx[i][::, 1])
        )

    return coo_array((data, (row, col)))


if __name__ == '__main__':

    a = np.array([[0.1, 0.2],
                  [0.3, 0.4]])

    b = np.array([[0.5, 0.6],
                [0.7, 0.8]])

    c = np.array([[0.9, 0.11, 0.12],
                [0.71, 0.81, 0.55],
                [0.71, 0.81, 0.55]])     
    
    print(
        build_assembly_coo_array([a, b, c]).toarray())

