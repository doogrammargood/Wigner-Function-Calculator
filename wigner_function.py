import itertools
import numpy as np
import math
import traceback
from sl_class import *
from finite_matrix_class import *
#from wigner_function import *

def parity_operator(p):
    array=[[1.0 if c== (-r)%p else 0.0 for c in range(p)]for r in range(p)]
    return np.matrix(array)

def parity_operator_general(p,n):
    array=[[1.0 if c== -r else 0.0 for c in finite_field_element.list_elements(p,n)]for r in finite_field_element.list_elements(p,n)]
    return np.matrix(array)

def weilX(x,p):
    array=[[ 1.0 if c==(r-x)%p else 0.0 for c in range(p)]for r in range(p)]
    return np.matrix(array)
def weilZ(z,p):
    array=[[np.e**(2*np.pi*1j*z*c/p ) if c ==r else 0.0 for c in range(p)]for r in range(p)]
    return np.matrix(array)

def weil_elementary( z, x, p ):
    #assumes that n == 1
    x=x%p
    z=z%p
    return np.e**(-modinv(2,p)*2*np.pi*1j*x*z/ p)* ( weilZ(z,p) @ weilX(x,p) )

def weil_general(pt):
    p = pt.x.p
    mat = None
    for a,b in zip(pt.x.coordinates, pt.y.coordinates):
        new_mat = weil_elementary(a,b,p)
        if mat is None:
            mat = new_mat
        else:
            mat = np.kron(new_mat, mat)
    return mat


def phase_pt_elementary(x,p):
    #returns phase point operator A(x) for the discrete Wigner function
    #assumes p is an odd prime.
    x0 = x[0]
    x1 = x[1]
    if (x0,x1)==(0,0):
        return parity_operator(p)
    else:
        return weil_elementary( 2*x0,2*x1 , p ) @ parity_operator(p)

def phase_pt_general(pt):
    if pt.x.n >1:
        point = point_of_plane ( (finite_matrix.convert_to_dual(pt.x), pt.y) )
    else:
        point = pt
    new_pt = point_of_plane((point.x + point.x, point.y + point.y))
    return weil_general(new_pt) @ parity_operator_general(point.x.p,point.x.n)

def discrete_wig_fuct(x,mat):
    p=x.x.p
    n=x.x.n
    return np.real(np.trace( phase_pt_general(x) @ mat)) * (p**-n)

def test_functions():
    p,n = 5,1
    r= random_pure_state(p,n)
    #print(phase_pt_elementary((0,1),p))
    pt = point_of_plane( (finite_field_element([0,1],5,1),finite_field_element.zero(5,1) ))
    print (discrete_wig_fuct(pt,r))

def debug_test():
    dim =5
    a = np.matrix(weil_elementary(2,2,dim))
    b = np.matrix(weil_elementary(2,2,dim))
    c = np.matmul(a,b)
    def symp(a,b):
        return (a[0]*b[1]-b[0]*a[1])%dim
    d = np.e**(modinv(2,dim)*(symp((2,2),(2,2)) )*2*np.pi*1j / dim)* np.matrix(weil_elementary(4,4,dim))
    assert np.allclose(d,c)
    print(c)
    print(d)
#debug_test()
#kron_test()
#test_functions()

#example_from_gross()
#test_functions()
