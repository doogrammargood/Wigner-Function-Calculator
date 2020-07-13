import numpy as np
def is_unitary(m):
    return np.allclose(np.eye(m.shape[0]), m.H * m)
def is_hermitian(m):
    return np.allclose(np.zeros(m.shape[0]), m.H - m)
