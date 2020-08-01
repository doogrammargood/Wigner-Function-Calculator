import numpy as np
def is_unitary(m):
    return np.allclose(np.eye(m.shape[0]), m.H * m)
def is_hermitian(m):
    return np.allclose(np.zeros(m.shape[0]), m.H - m)
def eigen_vector_of_hermitian(m):
    w, v = np.linalg.eigh(m)
    print(w)
    return v[:,-1]
