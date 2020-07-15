from finite_sp_matrix_class import *
from density_matrix_functions import *
from wigner_function import *
from grid_class import *
from covariance_functions import *
finite_matrix.load_dual_basis_matrices()
# def sl_matrix_to_prime_matrix(sl,p,n):
#     #takes an sl_2 matrix and produces the equivalent prime matrix.
#     prime_matrix = sl.to_prime_field_matrix()
#     convert = finite_matrix.dual_basis_conversion[(p,n)]
#     I = finite_matrix.identity(n,p,1)
#     O = finite_matrix.zero(n,p,1)
#     factor1 = finite_matrix([[convert,O],[O,I]])
#     factor2 = finite_matrix([[convert.inverse(),O],[O,I]])
#     prime_matrix = factor1 * prime_matrix * factor2
#     prime_matrix = finite_sp_matrix(prime_matrix)
#     return prime_matrix
    #unitary = prime_matrix.weil_representation()


def test_representation():# This test doesn't do much..
    p,n = 7,1
    I = finite_field_element.one(p,n)
    O = finite_field_element.zero(p,n)
    #s1 = finite_sp_matrix([[I,O],[I,I]])
    for s1 in finite_sp_matrix.list_sl_2(p,n):
        for s2 in finite_sp_matrix.list_sl_2(p,n):
            s1p = s1
            s2p = s2
            unitary1 = s1p.weil_representation() @ s2p.weil_representation()
            unitary2 = (s1p*s2p).weil_representation()
            # print(unitary1)
            # print(unitary2)
            print(s1p*s2p)
            assert np.allclose(s2p.weil_representation() @ s2p.weil_representation(), (s2p*s2p).weil_representation())
            assert np.allclose(unitary1, unitary2)
    print("s2 exhausted")

def debugging_test():
    p,n = 5,1
    I = finite_field_element.one(p,n)
    O = finite_field_element.zero(p,n)
    T = finite_field_element([2],p,n)
    bad = finite_sp_matrix([[I,I],[I,T]])
    print(bad.determinant())
    print( (bad*bad).determinant())
    #print(bad.weil_representation())
    assert np.allclose(bad.weil_representation() @ bad.weil_representation(), (bad*bad).weil_representation())

def debug_test2():
    J = finite_sp_matrix(finite_matrix.symplectic_form(2,5,1))
    for C in finite_matrix.list_invertible_matrices(1,5,1):
        C_mat = finite_sp_matrix.C_type_matrix(C)
        assert np.allclose(J.B_type_weil() @ C_mat.C_type_weil() @ J.B_type_weil(),  (J * C_mat * J).C_type_weil())
        print("ok")
    for A in finite_matrix.list_symmetric_matrices(1,5,1):
        A_mat = finite_sp_matrix.A_type_matrix(A)
        assert np.allclose( J.B_type_weil() @ A_mat.A_type_weil() @ J.B_type_weil(), (J * A_mat * J).weil_representation() )
#debug_test2()
#debugging_test()
#test_representation()


#----   TESTS RELATED TO CLIFFORD COVARIANCE ------#
def test_clifford_covariance_prime_field(sp, density_matrix, p, n):
    unitary = sp.weil_representation()
    new_density = unitary @ density_matrix @ unitary.H
    grid = grid_element(density_matrix, p, n, multiparticle = True)
    new_grid = grid_element(new_density, p, n, multiparticle = True)

    for x in finite_field_element.list_elements(p,n):
        for y in finite_field_element.list_elements(p,n):
            pt = point_of_plane((x,y))
            original_value = grid.get_value(pt)
            new_pt = sp * pt
            val1 = new_grid.get_value(new_pt)
            val2 = grid.get_value(pt)
            assert np.allclose(val1,val2)
def test_clifford_covariance_sl_2(sl, density_matrix,p,n):
    grid = grid_element(density_matrix, p, n)
    sp = sl.sl_matrix_to_prime_matrix(p,n)
    unitary = sp.weil_representation()
    new_density = unitary @ density_matrix @ unitary.H
    new_grid = grid_element(new_density, p, n)
    for x in finite_field_element.list_elements(p,n):
        for y in finite_field_element.list_elements(p,n):
            pt = point_of_plane((x,y))
            original_value = grid.get_value(pt)
            new_pt = sl * pt
            val1 = new_grid.get_value(new_pt)
            val2 = grid.get_value(pt)
            assert np.allclose(val1,val2)

def test_clifford_covariance():
    for s in finite_sp_matrix.list_sl_2(3,2):
        density_matrix = random_pure_state(3,2)
        test_clifford_covariance_sl_2(s,density_matrix,3,2)
        print(s)
    print("covariance sl2 passed")
    for s in finite_sp_matrix.list_sl_2(3,2):
        prime_matrix = sl_matrix_to_prime_matrix(s,3,2)
        density_matrix = random_pure_state(3,2)
        test_clifford_covariance_prime_field(prime_matrix,density_matrix,3,2)

def test_generate_state_from_sl():
    p,n = 3,2
    for sl in [s for s in finite_sp_matrix.list_sl_2(p,n) if s.order()==p**n+1]:
        sp = sl.sl_matrix_to_prime_matrix(p,n)
        U = sp.weil_representation()
        for i in range(p**n):
            density_matrix = state_from_sl(sl,i,p,n)
            assert np.allclose(density_matrix, U @ density_matrix @ U.H)
            grid = grid_element(density_matrix, p, n)
            for x in finite_field_element.list_elements(p,n):
                for y in finite_field_element.list_elements(p,n):
                    pt = point_of_plane((x,y))
                    print(grid.get_value(pt))
                    print(grid.get_value(sl*pt))
                    assert np.allclose(grid.get_value(pt), grid.get_value(sl*pt))
test_generate_state_from_sl()
# p,n = 5,1
# for sl in [s for s in finite_sp_matrix.list_sl_2(p,n) if s.order()==p**n+1]:
#     density_matrix = state_from_sl(sl,0,p,n)
#     test_clifford_covariance_prime_field(sl, density_matrix, p, n)
test_clifford_covariance()


#----
