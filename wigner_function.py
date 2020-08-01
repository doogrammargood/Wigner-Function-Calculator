import itertools
import numpy as np
import math
import traceback
from affine_plane_class import *
from finite_matrix_class import *

#This file contains the wigner function, and associated functions. Functions described as 'general' apply to all finite fields.

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

def phase_pt_general(pt, multiparticle = False):
    #Multiparticle is False when we calculate the Wigner Function of a single d^n dimensional particle.
    #True when we calculate the Wigner Function of n d-dimensional particles.
    if pt.x.n >1:
        if not multiparticle:
            point = point_of_plane ( (finite_matrix.convert_to_dual(pt.x), pt.y) )
        else:
            point = point_of_plane ( (pt.x, pt.y) )
    else:
        point = pt
    new_pt = point_of_plane((point.x + point.x, point.y + point.y))
    return weil_general(new_pt) @ parity_operator_general(point.x.p,point.x.n)

def discrete_wig_fuct(x,mat,multiparticle = False):
    p=x.x.p
    n=x.x.n
    return np.real(np.trace( phase_pt_general(x, multiparticle) @ mat)) * (p**-n)

def discrete_wig_fuct_pure_state(x,vect,multiparticle = False):
    vect = np.array(vect)[0]
    p=x.x.p
    n=x.x.n
    d=len(vect)
    if not multiparticle and n > 1:
        pt = point_of_plane ( (finite_matrix.convert_to_dual(x.x), x.y) )
    else:
        pt = x
    total = 0
    two_inv = finite_field_element([(p+1)//2],p,1)
    x_vec = pt.x.to_vector()
    y_vec = pt.y.to_vector()
    def vec_to_int(vect):
        tot = 0
        p = vect[0].p
        for i,v in enumerate(vect):
            tot += int(v)*p**i
        return tot
    for zeta in finite_field_element.list_elements(p,n):
        zeta_vec = zeta.to_vector()
        phase = np.e**(-2*np.pi*1j*int(reduce(lambda a,b: a+b, [a*b for a,b in zip(zeta_vec, x_vec)]))/p )
        factor1 = np.conj(vect[vec_to_int([a+two_inv*b for a,b in zip(y_vec, zeta_vec)])])
        factor2 = vect[vec_to_int([a-two_inv*b for a,b in zip(y_vec, zeta_vec)])]
        total += phase * factor1 *factor2
    to_return = total/d
    #assert np.allclose(to_return, np.real(to_return))
    return np.real(to_return)
def pure_state_from_density_matrix(mat):
    assert np.allclose(mat @ mat, mat)
    w,v = np.linalg.eigh(mat)
    i = len(v)-1
    vect= v[:,i].H
    assert np.allclose(np.matrix(vect).H @ np.matrix(vect), mat)
    return vect

def test_functions():
    p,n = 5,1
    r= random_pure_state(p,n)
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
