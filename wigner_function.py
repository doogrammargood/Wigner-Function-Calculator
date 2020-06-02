import itertools
import numpy as np
import math
import traceback
from sl_class import *

def is_unitary(m):
    return np.allclose(np.eye(m.shape[0]), np.matmul(m.H, m) )

def is_hermitian(m):
    return np.allclose(np.zeros(m.shape[0]), m.H - m)

def parity_operator(p):
    array=[[1.0 if c== (-r)%p else 0.0 for c in range(p)]for r in range(p)]
    return np.matrix(array)

def weilX(x,p):
    array=[[ 1.0 if c==(r-x)%p else 0.0 for c in range(p)]for r in range(p)]
    return np.matrix(array)
def weilZ(z,p):
    if p ==2:
        if z%2 == 1:
            return np.matrix([[1,0],[0,-1]])
        else:
            return np.matrix([[1,0],[0,1]])
    else:
        array=[[np.e**(2*np.pi*1j*z*c/p ) if c ==r else 0.0 for c in range(p)]for r in range(p)]
    return np.matrix(array)


#def weilZ(z,p)

def weil_elementary( z, x, p ):
    #assumes that n == 1
    if p==2:
        return 1j**(x*z)*np.matmul(weilZ(z,p),weilX(x,p))
    x=x%p
    z=z%p
    return np.e**(-modinv(2,p)*2*np.pi*1j*x*z/ p)*np.matmul(weilZ(z,p),weilX(x,p))


#-Start of changes to finite field elements
def phase_pt_elementary(x,p):
    #returns phase point operator A(x) for the discrete Wigner function
    #assumes p is an odd prime.
    x0 = x[0]
    x1 = x[1]
    if (x0,x1)==(0,0):
        return parity_operator(p)
    else:

        return np.matmul(weil_elementary( 2*x0,2*x1 , p ),parity_operator(p))

def phase_ptA(point):
    assert isinstance(point, point_of_plane)
    x_coords = point.x.coordinates
    y_coords = point.y.coordinates
    p = point.x.p
    matrix = None
    for x,y in zip(x_coords, y_coords):
        if matrix is None:
            matrix = phase_pt_elementary( (x,y), p)
        else:
            matrix =np.kron( matrix, phase_pt_elementary( (x,y), p) )
    return matrix


def discrete_wig_fuct(x,mat):
    p=x.x.p
    n=x.x.n
    return np.real(np.trace(np.matmul(phase_ptA(x), mat))) * 1/(p**n)

def random_pure_state(p,n):
    # for line in traceback.format_stack():
    #     print(line.strip())
    vector1 = np.array(np.random.rand(p**n)) - np.array([0.5 for i in range(p**n)])
    vector2 = 1j*( np.array(np.random.rand(p**n)) - np.array([0.5 for i in range(p**n)]) )
    vector = vector1 + vector2
    vector = np.matrix(vector)
    vector = vector / np.linalg.norm(vector)
    #vector = np.matmul(weil_elementary(0,2,p),vector.T).T
    #print (vector)
    return np.matmul(vector.H,vector)

def zero_state(p,n):
    vector = np.matrix([1.0 if i ==0 else 0 for i in range(p**n)])
    return np.matmul(vector.H,vector)

def state_with_special_order(p, i=4):#assume n = 1 for now..
    s = sl_matrix.gen_with_order(p,1).__next__()
    unitary = unitary_from_sl(s,p)
    eig = np.linalg.eig(unitary)
    eig = eig[1][:,i].T
    mat = np.matrix(eig)
    mat = np.matmul(mat.H,mat)
    return mat



def unitary_from_sl(mat, p):
    #Assumes p is prime.
    alpha = int(mat.matrix[0][0])
    beta = int(mat.matrix[0][1])
    gamma = int(mat.matrix[1][0])
    epsilon = int(mat.matrix[1][1])
    omega = np.exp((2*np.pi*1j)/float(p))
    tau = omega ** modinv(2,p)
    if p==2:
        tau = 1j
    array = []
    if beta == 0:
        array = [[tau**(alpha*gamma*r**2) if (alpha*r)%p==c else 0 for c in range(p)] for r in range(p)]
        return np.matrix(array)
    else:
        for j in range(p):
            row = []
            for k in range(p):
                #print (alpha*k**2 - 2*j*k + epsilon*j**2)%p *(beta**-1)
                row.append(p**-0.5 * tau** ( (alpha*k**2 - 2*j*k + epsilon*j**2) *(beta**-1) ))
            array.append(row)
        return np.matrix(array)

def test_functions():
    p,n = 5,1
    r= random_pure_state(p,n)
    #print(phase_pt_elementary((0,1),p))
    pt = point_of_plane( (finite_field_element([0,1],5,1),finite_field_element.zero(5,1) ))
    print (discrete_wig_fuct(pt,r))
#test_functions()
def example_from_gross():
    p = 3
    vector = np.matrix([0. ,1. ,-1. ])
    vector = 2**(-.5)*vector

    vectors = [vector]
    vectors.append(np.matmul(weil_elementary(-1, 0,p),vector.T).T)
    vectors.append(np.matmul(weil_elementary(-1, -1,p),vector.T).T)

    dense_matrix = 1/3. * np.matmul(vectors[0].H, vectors[0]) + 1/3. * np.matmul(vectors[1].H, vectors[1]) +1/3. * np.matmul(vectors[2].H, vectors[2])
    for x in finite_field_element.list_elements(3,1):
        for y in finite_field_element.list_elements(3,1):
            pt = point_of_plane((x,y))
            print(discrete_wig_fuct(pt, dense_matrix))
        print (" ")
#example_from_gross()
#test_functions()
