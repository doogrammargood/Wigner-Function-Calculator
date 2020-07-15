import itertools
import numpy as np
import math
import traceback
#from sl_class import *
from covariance_functions import *
from wigner_function import *
def is_unitary(m):
    return np.allclose(np.eye(m.shape[0]), np.matmul(m.H, m) )

def is_hermitian(m):
    return np.allclose(np.zeros(m.shape[0]), m.H - m)

def matrix_from_list(vector):
    mat = np.matrix(vector)
    mat = mat / np.linalg.norm(mat)
    return np.matmul(mat.H,mat)

def random_pure_state(p,n):
    # for line in traceback.format_stack():
    #     print(line.strip())
    vector1 = np.array(np.random.rand(p**n)) - np.array([0.5 for i in range(p**n)])
    vector2 = 1j*( np.array(np.random.rand(p**n)) - np.array([0.5 for i in range(p**n)]) )
    vector = vector1 + vector2
    return matrix_from_list(vector)
def zero_state(p,n):
    vector = np.array([1.0 if i ==0 else 0 for i in range(p**n)])
    #vector = np.array(weil_elementary(0,1,p)).dot(vector)
    return matrix_from_list(vector)

def maximally_mixed_state(p,n):
    identity = np.identity(p**n)
    return 1/(p**n) * np.matrix(identity)

def super_position_state(p,n):
    vector = np.array([1.0 for i in range(p**n)])
    #vector = np.array(weil_elementary(1,0,p)).dot(vector)
    #vector = np.matmul(weil_elementary( -1, 0,p),vector.T).T
    return matrix_from_list(vector)
def super_position_zero_negative(p,n):
    vector = [-1 if i ==0 else 1 for i in range(p**n)]
    return matrix_from_list(vector)

def super_position_state_negatives(p,n):
    #returns a superposition of all positions
    vectors = [ np.matrix([1.0 if i ==j else 0 for i in range(p**n)])  for j in range(p**n)]
    total = np.zeros(p**n)
    for i,v in enumerate(vectors):
        if i%2==0:
            total = total + v
        else:
            total = total - v
    mat = np.matrix(total) / np.linalg.norm(total)
    return np.matmul(mat.H,mat)

#def superpositon_state_zero_negative

def cat_state(p,n):
    def all_equal(l):
        return all([l[i]==l[i+1] for i in range(len(l)-1)])
    vector = [1.0 if all_equal(finite_field_element.from_int(i,p,n).coordinates) else 0 for i in range(p**n)]
    return matrix_from_list(vector)

def state_from_sl(sl, index, p, n):
    unitary = unitary_from_sl(sl, p, n)
    assert(is_unitary(unitary))
    eig = np.linalg.eig(unitary)
    assert np.allclose(unitary @ eig[1][:,index], eig[0][index] * eig[1][:,index])
    eig = eig[1][:,index]
    #print(eig)
    eig = np.matrix(eig)
    assert np.allclose((unitary @ eig).H, eig.H @ unitary.H)
    assert np.allclose(eig @ eig.H, unitary @ eig @ eig.H @ unitary.H)
    return np.matmul(eig,eig.H)
    #return matrix_from_list(eig)

def matrix_from_gross():
    p=3
    vector = np.array([0., 1., -1. ])
    vector = 2**(-.5)*vector
    vectors = [vector]
    #vectors.append(np.matmul(weil_elementary( -1, 0,p),vector.T).T)
    # print(np.array(weil_elementary(-1,0,p)))
    # print(np.array(vector))
    vectors.append(np.array(weil_elementary(1,0,p)).dot(np.array(vector)))
    #vectors.append(np.matmul(weil_elementary(-1, -1,p),vector.T).T)
    vectors.append(np.array(weil_elementary(1,-1,p)).dot(np.array(vector)))
    vectors = [np.matrix(v) for v in vectors]
    c = matrix_from_list(vector)
    a = weil_elementary(-1,0,p) @ c @ weil_elementary(1,0,p)
    b = weil_elementary(-1,-1,p) @ c @ weil_elementary(1,1,p)
    #dense_matrix = 1/3. * np.matmul(vectors[0].H, vectors[0]) + 1/3. * np.matmul(vectors[1].H, vectors[1]) +1/3. * np.matmul(vectors[2].H, vectors[2])
    #return matrix_from_list(vector)
    dense_matrix = 1/3*(a + b + c)
    return dense_matrix

def unitary_from_sl(sl, p, n):
    sp = sl.sl_matrix_to_prime_matrix(p,n)
    unitary = sp.weil_representation()
    return unitary

# def unitary_from_sl(mat, p): This is deprecated. It only handled prime matrices.
#     #Assumes p is prime.
#     alpha = int(mat.matrix[0][0])
#     beta = int(mat.matrix[0][1])
#     gamma = int(mat.matrix[1][0])
#     epsilon = int(mat.matrix[1][1])
#     omega = np.exp((2*np.pi*1j)/float(p))
#     tau = omega ** modinv(2,p)
#     if p==2:
#         tau = 1j
#     array = []
#     if beta == 0:
#         array = [[tau**(alpha*gamma*r**2) if (alpha*r)%p==c else 0 for c in range(p)] for r in range(p)]
#         return np.matrix(array)
#     else:
#         for j in range(p):
#             row = []
#             for k in range(p):
#                 #print (alpha*k**2 - 2*j*k + epsilon*j**2)%p *(beta**-1)
#                 row.append(p**-0.5 * tau** ( (alpha*k**2 - 2*j*k + epsilon*j**2) *(beta**-1) ))
#             array.append(row)
#         return np.matrix(array)

def test_states():
    epsilon = 10**-7
    test_primes = [5]
    test_n = [2]
    #test_functions = [random_pure_state, zero_state, maximally_mixed_state,super_position_state,cat_state,super_position_state_negatives]
    test_functions = [super_position_state_negatives]
    for p,n in itertools.product(test_primes, test_n):
        for f in test_functions:
            density_matrix = f(p,n)
            print(f.__name__)
            assert is_hermitian(density_matrix)
            assert abs(np.trace(density_matrix)-1)< epsilon
            assert np.all(np.linalg.eigvals(density_matrix) > -epsilon)
#test_states()
