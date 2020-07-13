from finite_sp_matrix_class import *
from density_matrix_functions import *
from wigner_function import *
from grid_class import *
finite_matrix.load_dual_basis_matrices()
def sl_matrix_to_prime_matrix(sl,p,n):
    #takes an sl_2 matrix and produces the equivalent prime matrix.
    prime_matrix = sl.to_prime_field_matrix()
    convert = finite_matrix.dual_basis_conversion[(p,n)]
    I = finite_matrix.identity(n,p,1)
    O = finite_matrix.zero(n,p,1)
    factor1 = finite_matrix([[convert,O],[O,I]])
    factor2 = finite_matrix([[convert.inverse(),O],[O,I]])
    prime_matrix = factor1 * prime_matrix * factor2
    prime_matrix = finite_sp_matrix(prime_matrix)
    return prime_matrix
    #unitary = prime_matrix.weil_representation()


def test_representation():
    for s1 in finite_sp_matrix.list_sl_2(3,2):
        for s2 in finite_sp_matrix.list_sl_2(3,2):
            s1p =  sl_matrix_to_prime_matrix(s1,3,2)
            s2p = sl_matrix_to_prime_matrix(s2,3,2)

            unitary1 = s1p.weil_representation()*s2p.weil_representation()
            unitary2 = (s1p*s2p).weil_representation()
            assert np.allclose(unitary1,unitary2)
            #print("accurate")
        print("s2 exhausted")


def test_clifford_covariance(sp, density_matrix, p, n):
    #print(sp)
    unitary = sp.weil_representation()
    #unitary = sp.C_type_weil()
    #unitary = sp.A_type_weil()
    #print(unitary)
    new_density = unitary @ density_matrix @ unitary.H
    grid = grid_element(density_matrix, p, n, multiparticle = True)
    new_grid = grid_element(new_density, p, n, multiparticle = True)

    for x in finite_field_element.list_elements(p,n):
        for y in finite_field_element.list_elements(p,n):
            pt = point_of_plane((x,y))
            original_value = grid.get_value(pt)
            new_pt = sp.inverse() * pt
            val1 = new_grid.get_value(pt)
            val2 = grid.get_value(new_pt)
            # print (original_value)
            # print (new_value)
            # print ("--")
            assert np.allclose(val1,val2)

test_representation()
# for s in finite_sp_matrix.list_sl_2(3,2):
#     prime_matrix = sl_matrix_to_prime_matrix(s,3,2)
#     density_matrix = random_pure_state(3,2)
#     test_clifford_covariance(prime_matrix,density_matrix,3,2)
#
# n = 2
# for A in finite_matrix.list_symmetric_matrices(n,3,1):
#     prime_matrix= finite_sp_matrix.A_type_matrix(A)
#     #prime_matrix = sl_matrix_to_prime_matrix(s,3,1)
#     print(prime_matrix)
#     density_matrix = random_pure_state(3,n)
#     print(density_matrix)
#     test_clifford_covariance(prime_matrix,density_matrix,3,n)
# for C in finite_matrix.list_invertible_matrices(2,5,1):
#     prime_matrix= finite_sp_matrix.C_type_matrix(C)
#     #prime_matrix = sl_matrix_to_prime_matrix(s,3,1)
#     print(prime_matrix)
#     density_matrix = random_pure_state(5,2)
#     #print(density_matrix)
#     test_clifford_covariance(prime_matrix,density_matrix,5,2)
