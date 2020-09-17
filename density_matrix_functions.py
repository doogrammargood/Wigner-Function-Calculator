import itertools
import numpy as np
import math
import traceback
#from sl_class import *
from covariance_functions import *
from wigner_function import *
from finite_sp_matrix_class import *
def is_unitary(m):
    return np.allclose(np.eye(m.shape[0]), np.matmul(m.H, m) )

def is_hermitian(m):
    return np.allclose(np.zeros(m.shape[0]), m.H - m)

def matrix_from_list(vector):
    mat = np.matrix(vector)
    mat = mat / np.linalg.norm(mat)
    return np.matmul(mat.H,mat)

def random_pure_state(p,n):
    vector1 = np.array(np.random.normal(size = p ** n))
    vector2 = 1j*( np.array(np.random.normal( size = p ** n)) )
    vector = vector1 + vector2
    return matrix_from_list(vector)
def random_pure_state_strange(p,n):
    #This is not the usual distrubition of random states.
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
    #vector = np.array(weil_elementary(1,0,p)00).dot(vector)
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

def stabilizer_state_from_line(l):
    #TODO: use covariance property to define these faster.
    p=l.coefficients[0].p
    n=l.coefficients[0].n
    total = sum([phase_pt_general(pt) for pt in l.gen_points()])
    return total

def superpositon_of_stabilizers(p,n):
    zero = finite_field_element.zero(p,n)
    one = Z = finite_field_element.one(p,n)
    l = line_of_plane((zero,one,one))
    sl = finite_sp_matrix.get_element_of_sl_2_from_field_extension(p,n)
    total = sum([stabilizer_state_from_line(lp) for lp in sl.orbit(l)])
    return matrix_from_list(total.T)


def cat_state(p,n):
    def all_equal(l):
        return all([l[i]==l[i+1] for i in range(len(l)-1)])
    vector = [1.0 if all_equal(finite_field_element.from_int(i,p,n).coordinates) else 0 for i in range(p**n)]
    return matrix_from_list(vector)

def labyrinth_state(p,n):
    sl = finite_sp_matrix.get_element_of_sl_2_from_field_extension(p,n)
    return state_from_sl(sl, 0, p, n)

def mirror_state(p,n):
    m = finite_field_element([1 if i == 1 else 0 for i in range(2*n)], p, 2*n)
    mat = finite_matrix.from_finite_field_element(m, new_n = n)
    esl = mat ** ((p**n-1)/ 2)
    one = finite_field_element.one(p,n)
    zero = finite_field_element.zero(p,n)
    J = finite_matrix([[one, zero],[zero, -one]])
    anti = finite_sp_matrix(esl*J)._weil_representation_case_2()
    X = antiunitary_normal_form(anti)
    return np.matrix(X[1][0]) @ np.matrix(X[1][0]).H

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

def unitary_from_sl(sl, p, n):
    sp = sl.sl_matrix_to_prime_matrix(p,n)
    unitary = np.matrix(sp.weil_representation())
    return unitary

def state_from_collapse_to_labyrinth(p,n):
    density_matrix = zero_state(p,n)
    sl = finite_sp_matrix.get_element_of_sl_2_from_field_extension(p,n)
    unitary = unitary_from_sl(sl, p, n)
    return collapse_to_eigenstate(density_matrix, unitary, 1)

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
