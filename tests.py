import unittest
from unittest.mock import Mock, patch
from finite_sp_matrix_class import finite_sp_matrix
from finite_matrix_class import *
from functional_class import *
finite_matrix.load_dual_basis_matrices()

class SimpleFiniteSpMatrixClass(unittest.TestCase):

    def setUp(self):
        self.p,self.n = 5,1
        self.one = finite_field_element.one(self.p, self.n)
        self.zero = finite_field_element.zero(self.p, self.n)
        pass

    def test_finite_sp_matrix_init(self):
        self.one = finite_field_element.one(3, 1)
        self.zero = finite_field_element.zero(3, 1)
        test_finite_sp_matrix = finite_sp_matrix([[self.zero, -self.one],[self.one, self.zero]])
        self.assertEqual(test_finite_sp_matrix.p, 3)
        self.assertEqual(test_finite_sp_matrix.n, 1)

    def test_finite_sp_matrix_init_fail(self):
        self.one = finite_field_element.one(3, 1)
        self.zero = finite_field_element.zero(3, 1)
        with self.assertRaises(Exception):
            finite_sp_matrix([[self.zero, self.one],[self.one, self.zero]])


    def test_block_matrix_decompose(self):
        self.one = finite_field_element.one(3, 1)
        self.zero = finite_field_element.zero(3, 1)
        self.finite_sp_matrix = finite_sp_matrix([[self.zero, -self.one],[self.one, self.zero]])
        test_ret = self.finite_sp_matrix.block_matrix_decompose()
        self.assertEqual(len(test_ret), 2)
        self.assertEqual(len(test_ret[0]), 2)
        self.assertEqual(len(test_ret[1]), 2)

    # @patch('finite_sp_matrix_class.finite_sp_matrix.B_type_weil')# I just don't understand the need for mocks- Just check the required equations are true?
    # @patch('finite_matrix_class.finite_matrix.symplectic_form')
    # @patch('finite_sp_matrix_class.finite_sp_matrix.factorize_lower')
    # @patch('finite_sp_matrix_class.finite_sp_matrix.factorize_upper')
    # @patch('finite_sp_matrix_class.finite_sp_matrix.block_matrix_decompose')
    # @patch('finite_sp_matrix_class.finite_sp_matrix.factorize')
    # def test_weil_representation(self, mock_factorize, mock_block_matrix_decompose, mock_factorize_upper, mock_factorize_lower, mock_symplectic, mock_B_type_weil):
    #     # you can set more mock return values for factorize, block_matrix_decompose, etc. as you see fit.
    #     self.one = finite_field_element.one(3, 1)
    #     self.zero = finite_field_element.zero(3, 1)
    #     self.finite_sp_matrix = finite_sp_matrix([[self.zero, -self.one],[self.one, self.zero]])
    #     mock_symplectic.return_value = [[self.zero, -self.one],[self.one, self.zero]]
    #     ret_val = self.finite_sp_matrix.weil_representation()
    #     print(ret_val)
    #
    # def test_weil_represenation(self):
    #     p,n = 3,2
    #     unit_eq = None
    #     for s1 in finite_sp_matrix.list_sl_2(p,n):
    #         for s2 in finite_sp_matrix.list_sl_2(p,n):
    #             self.assertTrue(np.allclose((s1*s2).weil_representation(),s1.weil_representation()@s2.weil_representation()))
    #
    # def test_sl_from_extension(self):
    #     p_n_pairs = [(3,2), (5,2), (7,2), (11,2), (3,3), (5,3), (7,3), (3,4)]
    #     percents = [0.001, 0.05, 0.1, 0.2, 0.4, 0.5, 0.75]
    #     print("p, n, i, percent, q_val, c_val, ratio, numlines")
    #     for p,n in p_n_pairs:
    #         sl = finite_sp_matrix.get_element_of_sl_2_from_field_extension(p,n)
    #         for i in range(p**n):
    #             density_matrix = state_from_sl(sl,i,p,n)
    #             grid = grid_element(density_matrix, p, n)
    #             total_negativity = grid.total_negativity()
    #             for percent in percents:
    #                 functional = functional_on_grid.functional_from_sl(sl, grid, percent)
    #                 q_val = functional * grid
    #                 c_val = functional.evaluate_classical(grid)[1]
    #                 num_lines = len(list(functional.list_lines()))
    #                 print([p, n, i, percent, q_val, c_val, q_val/c_val, num_lines, total_negativity])

    # def test_functional_count_positive_points(self):
    #     p_n_pairs = [(3,2), (5,2), (7,2), (11,2), (3,3), (5,3), (7,3), (3,4)]
    #     print("p, n, i, q_val, c_val, ratio, numlines")
    #     for p,n in p_n_pairs:
    #         sl = finite_sp_matrix.get_element_of_sl_2_from_field_extension(p,n)
    #         for i in range(p**n):
    #             density_matrix = state_from_sl(sl,i,p,n)
    #             grid = grid_element(density_matrix, p, n)
    #             total_negativity = grid.total_negativity()
    #             functional = functional_on_grid.functional_counting_positive_points(grid)
    #             q_val = (functional * grid)
    #             c_val = functional.evaluate_classical(grid)[1]
    #             c_max = functional.point_count
    #             num_lines = len(list(functional.list_lines()))
    #             print([p, n, i, q_val, c_val, q_val/c_val, c_max, total_negativity])
    def test_antiunitary_normal_form(self):
        def special_sqrt(omega):
            x = omega**(.5)
            if np.imag(x) * np.imag(omega) < 0:
                x = -x
            return x
        p,n = self.p, self.n
        one, zero = self.one, self.zero
        m = finite_field_element([1 if i == 1 else 0 for i in range(2*n)], p, 2*n)
        mat = finite_matrix.from_finite_field_element(m, new_n = n)
        esl = mat ** ((p**n-1)/ 2)
        J = finite_matrix([[one, zero],[zero, -one]])
        anti = finite_sp_matrix(esl*J)._weil_representation_case_2()
        X = antiunitary_normal_form(anti)
        vectors = reduce(lambda a,b: a+b, [X[key] for key in X.keys()], [])
        for i in vectors:
            assert np.allclose(1, np.linalg.norm(i))
        for i,j in itertools.combinations(vectors,2):
            assert np.allclose(np.vdot(i,j), 0)
        for key in X.keys():
            for e, vect in enumerate(X[key]):
                conj_key = [k for k in X.keys() if np.allclose(k, np.conj(key))]
                assert len(conj_key)==1
                conj_key = conj_key[0]
                pair_vect = special_sqrt(key)*anti@np.conj(vect)
                assert np.allclose(X[conj_key][e], pair_vect)
if __name__ == '__main__':
    unittest.main()
