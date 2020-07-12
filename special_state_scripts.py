from finite_sp_matrix_class import *
finite_matrix.load_dual_basis_matrices()
def produce_state(sl,p,n):
    prime_matrix = sl.to_prime_field_matrix()
    convert = finite_matrix.dual_basis_conversion[(p,n)]
    I = finite_matrix.identity(n,p,1)
    O = finite_matrix.zero(n,p,1)
    factor1 = finite_matrix([[I,O],[O,convert]])
    factor2 = finite_matrix([[I,O],[O,convert.inverse()]])
    prime_matrix = factor1 * prime_matrix * factor2
    prime_matrix = finite_sp_matrix(prime_matrix)
    return prime_matrix
    #unitary = prime_matrix.weil_representation()

for s in finite_sp_matrix.list_sl_2(3,2):
    print(produce_state(s,3,2))
