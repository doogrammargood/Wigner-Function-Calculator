import numpy as np
def is_unitary(m):
    return np.allclose(np.eye(m.shape[0]), m.H * m)
def is_hermitian(m):
    return np.allclose(np.zeros(m.shape[0]), m.H - m)
def eigen_vector_of_hermitian(m):
    w, v = np.linalg.eigh(m)
    print(w)
    return v[:,-1]

def pure_state_from_density_matrix(mat):
    assert np.allclose(mat @ mat, mat)
    w,v = np.linalg.eigh(mat)
    i = len(v)-1
    vect= v[:,i].H
    assert np.allclose(np.matrix(vect).H @ np.matrix(vect), mat)
    return vect

def collapse_to_eigenstate(density_matrix, unitary, eigenvalue):
    #Returns a density matrix which is the normalized projection onto the eigenspace of unitary with the given eigenvalue.
    #As if starting in state 'density matrix', and applying the eigenvalue-finding method for 'unitary', then post-selecting for 'eigenvalue'
    #For now, assume density matrix is rank 1
    assert np.allclose(density_matrix, density_matrix @ density_matrix)

    pure_state = pure_state_from_density_matrix(density_matrix)
    eig = np.linalg.eig(unitary)
    indices = [i for i, val in enumerate(eig[0]) if np.allclose(val, eigenvalue)]
    print(np.shape(pure_state), np.shape(eig[1][:,indices[0]]))
    eig_tot = sum([(np.matrix(eig[1][:,i]) @ np.matrix(pure_state)) * np.matrix(eig[1][:,i]) for i in indices])
    eig_tot = eig_tot/np.linalg.norm(eig_tot)
    return np.matmul(eig_tot,eig_tot.H)

def antiunitary_normal_form(A):
    #From https://aip.scitation.org/doi/pdf/10.1063/1.1703672
    #Our variable A = AK from the paper.
    def special_sqrt(omega):
        x = omega**(.5)
        if np.imag(x) * np.imag(omega) < 0:
            x = -x
        return x

    eig = np.linalg.eig(A @ np.conj(A))
    ones_index = [i for i, val in enumerate(eig[0]) if np.allclose(val, 1)]
    dict = {1: ones_index} #from eigenvalues to a list of indices
    for i, value in enumerate(eig[0]): #prepares dict
        key = [key for key in dict.keys() if np.allclose(value, key) and not np.allclose(value,1)]
        if len(key)>0:
            dict[key[0]].append(i)
        elif not np.allclose(value, 1):
            dict[value]=[i]
    return_dict = {}
    for key in dict.keys():
        if any([np.allclose(key, k) for k in return_dict.keys()]):
            continue #key is already in return_dict
        if key == 1:
            return_dict[1]=[]
            indices = dict[key]
            vects = [eig[1][:,i] for i in indices]
            for v in vects:
                x = A@np.conj(np.vstack(v))+np.vstack(v)
                #print("xis ",x)
                if not np.allclose(np.linalg.norm(x),0):
                    x = x/ np.linalg.norm(x)
                    return_dict[key].append(x)
                else:
                    return_dict[key].append(1j*np.vstack(v))
        elif not np.allclose(key,-1):
            if np.imag(key) > 0:
                return_dict[key] = [np.vstack(eig[1][:,index]) for index in dict[key]]
                new_key = np.conj(key)
                return_dict[new_key] = [special_sqrt(key)*A@np.conj(np.vstack(v)) for v in return_dict[key]]
        else:
            print("untested") #In the test cases, this never occurs, so I haven't been able to test it. It probably works.
            assert(False)
            return_dict[key] = []
            current_values = [[],[]]
            indices = dict[key]
            for index in indices:
                vector = np.vstack(eig[1][:,index])
                s = sum(current_values[0])+sum(current_values[1])
                vector = vector - np.vdot(s,vector)*s
                vector = vector / np.linalg.norm(vector)
                current_values[0].append(vector)
                conj_vector = 1j*A @ np.conj(vector)
                current_values[1].append(conj_vector)
            return_dict[key] = current_values[0] + current_values[1]
    return return_dict
