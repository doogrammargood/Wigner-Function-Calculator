from finite_matrix_class import *
from numpy_helpers import *
from affine_plane_class import *
from wigner_function import *
#finite_matrix.load_dual_basis_matrices()
#This class is a little bloated. It contains all the information for finite symplectic matrices,
#It also includes SL(2,p) as a special case.
class finite_sp_matrix(finite_matrix):
    def __init__(self, list_representation):
        if isinstance(list_representation, finite_matrix):
            #promote finite matrix to sp finite_matrix
            #assert list_representation.is_symplectic()
            self.elements = list_representation.elements
            self.p=list_representation.p
            self.n=list_representation.n
        else:
            super().__init__(list_representation)
            assert self.is_symplectic()

    def __pow__(self,exp):
        return finite_sp_matrix(super().__pow__(exp))

    def __mul__(self, other):
        if isinstance(other, point_of_plane):
            x = other.x
            y = other.y
            if len(self.elements[0])==  2: #when self is an element of sl_2
                col_matrix = finite_matrix([[x],[y]])
                new_col_matrix = self * col_matrix
                new_row_matrix = new_col_matrix.transpose()
                new_x = new_row_matrix.elements[0][0]
                new_y = new_row_matrix.elements[0][1]
                return point_of_plane((new_x, new_y))
            else: #for a prime field matrix.
                vector = x.to_vector() + y.to_vector()
                col_matrix = finite_matrix(vector)
                col_matrix = self * col_matrix #This time, we are multiplying matrices!
                row_matrix = col_matrix.transpose()
                new_vector = row_matrix.elements[0]
                new_x_vector = new_vector[:len(new_vector)//2]
                new_y_vector = new_vector[len(new_vector)//2:]
                new_x = finite_field_element.from_vector(new_x_vector)
                new_y = finite_field_element.from_vector(new_y_vector)
                return point_of_plane((new_x,new_y))
        elif isinstance(other, line_of_plane):
            #TODO: rewrite this in terms of the coefficients of line.
            iter = other.gen_points()
            pt1 = next(iter)
            pt2 = next(iter)
            new_pt1 = self * pt1
            new_pt2 = self * pt2
            return new_pt1.line_to(new_pt2)



        elif isinstance(other, finite_sp_matrix):
            return finite_sp_matrix(super().__mul__(other))
        else: #when other is a finite_matrix.
            return super().__mul__(other)
    def __neg__(self):
        return finite_sp_matrix(super().__neg__())
    def transpose(self):
        return finite_sp_matrix(super().transpose())
    def inverse(self):
        return finite_sp_matrix(super().inverse())

    def orbit(self, other):
        #assert isinstance(other, point_of_plane) or isinstance(other, line_of_plane) or isinstance(other, finite_sp_matrix)
        #Yay polymorphism! But maybe this should be in finite_matrix.
        to_return = [other]
        current = self*other
        count =0
        while not current == other:
            to_return.append(current)
            current = self * current
        return to_return

    def block_matrix_decompose(self):
        #returns a list-of-lists of matrices [[A,B][C,D]]
        size = len(self.elements)
        block_size = size//2
        A_elements = [[self.elements[j][i] for i in range(block_size) ]for j in range(block_size)]
        B_elements = [[self.elements[j][i] for i in range(block_size,size)] for j in range(block_size)]
        C_elements = [[self.elements[j][i] for i in range(block_size)] for j in range(block_size,size)]
        D_elements = [[self.elements[j][i] for i in range(block_size,size)] for j in range(block_size,size)]
        return [[finite_matrix(A_elements), finite_matrix(B_elements)],[finite_matrix(C_elements), finite_matrix(D_elements)]]

    def sl_matrix_to_prime_matrix(self,p,n):
        #takes an sl_2 matrix and produces the equivalent prime matrix.
        if n == 1:
            return self
        prime_matrix = self.to_prime_field_matrix()
        convert = finite_matrix.dual_basis_conversion[(p,n)]
        I = finite_matrix.identity(n,p,1)
        O = finite_matrix.zero(n,p,1)
        factor1 = finite_matrix([[convert,O],[O,I]])
        factor2 = finite_matrix([[convert.inverse(),O],[O,I]])
        prime_matrix = factor1 * prime_matrix * factor2
        prime_matrix = finite_sp_matrix(prime_matrix)
        return prime_matrix

    def factorize(self):
        #TODO: This method is not totally general. Maybe use the symplectic Gaussian Elimination?
        blocks = self.block_matrix_decompose()
        A=blocks[0][0]
        B=blocks[0][1]
        C=blocks[1][0]
        D=blocks[1][1]
        O = finite_matrix.zero(len(self.elements)//2,self.p,self.n)
        if D.determinant().is_zero():
            J = finite_sp_matrix.B_type_matrix(finite_matrix.identity(len(self.elements)//2,self.p,self.n))
            temp = self*J
            if C.determinant().is_zero():
                print("problem")
            fact = temp.factorize()
            return fact[0], fact[1], J.inverse()
        else:
            factor_1 = finite_sp_matrix.C_type_matrix(D.transpose())
            temp1 = factor_1 * self
            temp_blocks = temp1.block_matrix_decompose()
            Y = -temp_blocks[0][1]
            assert Y.transpose()==Y
            factor_2 = finite_sp_matrix.A_type_matrix(Y)
            lower = factor_2 * temp1
            final_blocks = lower.block_matrix_decompose()
            assert final_blocks[0][1].is_zero()
            upper = (factor_2*factor_1).inverse()
            #upper is triangular, lower is strictly triangular, and the last matrix is either J inverse or the identity
            return upper, lower, finite_matrix.identity(len(self.elements),self.p,self.n)

    def factorize_upper(self):
        #assume self is a block upper triangular matrix
        blocks = self.block_matrix_decompose()
        #assert blocks[1][0].is_zero()
        M = blocks[0][0]
        A = blocks[0][1]
        #assert not M.determinant().is_zero()
        factor1 = finite_sp_matrix.C_type_matrix(M)
        factor2 = finite_sp_matrix.A_type_matrix(M.inverse()*A)
        assert factor1 * factor2 == self
        return factor1, factor2

    def factorize_lower(self):
        #assume self is a strictly block lower triangular matrix.
        blocks = self.block_matrix_decompose()
        assert blocks[0][1].is_zero()
        assert blocks[0][0]==blocks[1][1] #these should both be the identity.
        A = blocks[1][0]
        factor = finite_sp_matrix.A_type_matrix(-A)
        J = finite_sp_matrix.B_type_matrix(finite_matrix.identity(len(self.elements)//2,self.p,self.n))
        assert J * factor * J.inverse() == self
        return J, factor, J.inverse()



    def A_type_weil(self):
        blocks = self.block_matrix_decompose()
        A = blocks[0][1]
        #assert self == finite_sp_matrix.A_type_matrix(A)
        size = len(blocks[0][0].elements)
        diag = []
        two_inv = finite_field_element([2],self.p,1).inverse()
        for ys in itertools.product(finite_field_element.list_elements(self.p,1), repeat = size):
            yps = list(ys)
            yps.reverse()
            col_matrix = finite_matrix([[y] for y in yps])
            row_matrix = col_matrix.transpose()
            result = finite_matrix([[two_inv]]) * (row_matrix * A * col_matrix)
            diag.append(result.character())
        return np.matrix(np.diag(diag))

    def C_type_weil(self):
        blocks = self.block_matrix_decompose()
        #assert blocks[0][1] == blocks[1][0]
        C = blocks[0][0]
        C_inv_trans = C.inverse().transpose()
        permutation = [] #the result should be a permutation matrix.
        size = len(blocks[0][0].elements)
        def col_matrix_to_int(col_matrix):
            vector = col_matrix.transpose().elements[0]
            ffe = finite_field_element.from_vector(vector)
            return int(ffe)

        for ys in itertools.product(finite_field_element.list_elements(self.p,1), repeat = size):
            yps = list(ys)
            yps.reverse()
            col_matrix = finite_matrix([[y] for y in yps])
            result = C_inv_trans * col_matrix
            permutation.append( col_matrix_to_int(result) )
        result = np.matrix([[1 if row == permutation[column] else 0 for column in range(self.p**size)] for row in range(self.p**size)])

        multiplier = legendre_symbol(int(C.determinant()), self.p)
        #multiplier = C.determinant().square_class()
        #print(multiplier)
        return multiplier*result

    def B_type_weil(self):
        #only defined for the symplectic form over a prime field
        #assert self == finite_matrix.symplectic_form(len(self.elements),self.p,1)
        neg_two_inv = -(finite_field_element([2],self.p,1).inverse())
        to_return = []
        guass_sum = 0
        size = len(self.elements)//2
        B = self.block_matrix_decompose()[0][1]
        #size = len(blocks[0][0].elements)
        for ys in itertools.product(finite_field_element.list_elements(self.p,1), repeat = size):
            current_col = []
            yps = list(ys)
            yps.reverse()
            col_y = finite_matrix([[y] for y in yps])
            row_y = col_y.transpose()
            guass_sum += (finite_matrix([[neg_two_inv]]) * (row_y * B * col_y) ).character()
            for xs in itertools.product(finite_field_element.list_elements(self.p,1), repeat = size):
                xps = list(xs)
                xps.reverse()
                col_x = finite_matrix([[x] for x in xps])
                #val = (finite_matrix([[neg_two_inv]]) * (row_y * B * col_x) ).character()
                val = ( (row_y * B * col_x) ).character()
                current_col.append(val)
            to_return.append(current_col)
        to_return = np.matrix(to_return).T
        return (1/guass_sum) * to_return

    def weil_representation(self):

        if len(self.elements) ==2:
            return self._weil_representation_case_2()
        #This is sketchy. For SL matrices, use https://arxiv.org/pdf/0909.5233.pdf
        upper, lower, last = self.factorize()
        upper_blocks = upper.block_matrix_decompose()

        up_factor1, up_factor2 = upper.factorize_upper()
        assert up_factor1 * up_factor2 == upper
        low_factor1, low_factor2, low_factor3 = lower.factorize_lower()
        assert up_factor1 * up_factor2 * low_factor1 * low_factor2 * low_factor3 * last == self
        J = finite_sp_matrix(finite_matrix.symplectic_form(len(self.elements),self.p,1))
        assert low_factor1 == J
        assert low_factor3 == J.inverse()
        weil_of_j = J.B_type_weil()
        #weil_of_j_inv = np.linalg.inv(weil_of_j)
        weil_of_j_inv = weil_of_j.H
        #assert np.allclose(np.eye(5,5), weil_of_j @ weil_of_j_inv)
        assert np.allclose((J*J).C_type_weil(),weil_of_j @ weil_of_j)
        #assert np.allclose(np.eye(5,5),  weil_of_j @ weil_of_j @ weil_of_j @ weil_of_j)
        assert is_unitary(up_factor1.C_type_weil())
        assert is_unitary(up_factor2.A_type_weil())
        assert is_unitary(weil_of_j)
        assert is_unitary(low_factor2.A_type_weil())
        #assert is_unitary()
        if last == J.inverse():
            to_return = up_factor1.C_type_weil() @ up_factor2.A_type_weil() @ weil_of_j @ low_factor2.A_type_weil() @ weil_of_j_inv @ weil_of_j_inv
        else:
            to_return = up_factor1.C_type_weil() @ up_factor2.A_type_weil() @ weil_of_j @ low_factor2.A_type_weil() @ weil_of_j_inv
        #assert is_unitary(to_return)
        return to_return

    def _weil_representation_case_2(self):
        #Following https://arxiv.org/pdf/0909.5233.pdf see equation 76.
        if len(self.elements) != 2:
            print ("Only call for 2x2 matrices.")
        p, n = self.p, self.n
        d = p**n

        two = finite_field_element.from_int(2,p,n)
        two_inv = int(two.inverse())
        tau = np.e**(np.pi *1j * (p+1)/p)

        one = finite_field_element.one(p,n)
        zero = finite_field_element.zero(p,n)
        K = finite_matrix([[zero, one], [one, zero]]) #Appleby's displacement operators are XZ, so we have to flip. TODO: prove this is correct
        other = K*self*K
        alpha, beta = other.elements[0]
        gamma, delta = other.elements[1]
        if beta.is_zero():
            factor = alpha.square_class()
            to_return = factor * np.array([[tau** (int((alpha*gamma*col**2).trace())) if row == alpha*col else 0 for col in finite_field_element.list_elements(p,n)]for row in finite_field_element.list_elements(p,n)])
        else:
            l_tilde = -(1j**-((n*(p+3)/2)))*(-beta).square_class()
            to_return = l_tilde/(d**0.5) * np.array([[ tau**(int(((alpha*col**2-two*row*col+delta*row**2)/beta).trace())) for col in finite_field_element.list_elements(p,n)]for row in finite_field_element.list_elements(p,n)])
        return to_return


    @classmethod
    def identity(cls, size,p,n):
        return finite_sp_matrix(finite_matrix.identity(size,p,n))

    #The following generators of the symplectic group come from https://www.math.wisc.edu/~shamgar/Small-Representations-Howe70th-Proceedings.pdf
    @classmethod
    def A_type_matrix(cls, A):
        #needs A to be symmetric
        assert A == A.transpose()
        size = len(A.elements)
        I = finite_matrix.identity(size, A.p, A.n)
        O = finite_matrix.zero(size, A.p, A.n)
        return finite_sp_matrix([[I,A],[O,I]])

    @classmethod
    def B_type_matrix(cls, B):
        #needs B to be symmetric and invertible
        assert B == B.transpose()
        assert not B.determinant().is_zero()
        size = len(B.elements)
        O = finite_matrix.zero(size, B.p, B.n)
        return finite_sp_matrix([[O,B],[-(B.inverse()),O]])

    @classmethod
    def C_type_matrix(cls, C):
        #assert not C.determinant().is_zero()
        size = len(C.elements)
        O = finite_matrix.zero(size, C.p, C.n)
        return finite_sp_matrix([[C,O],[O, C.inverse().transpose()]])

    @classmethod
    def list_sl_2(cls,p,n):
        I = finite_field_element.one(p,n)
        O = finite_field_element.zero(p,n)
        def part1():
            for A,B,C in itertools.product(finite_field_element.list_elements(p,n), finite_field_element.list_nonzero_elements(p,n), finite_field_element.list_elements(p,n)):
                yield finite_sp_matrix([[I,O],[A,I]])*finite_sp_matrix([[B,O],[O,B.inverse()]])*finite_sp_matrix([[I,C],[O,I]])
        def part2():
            for A,B in itertools.product(finite_field_element.list_elements(p,n), finite_field_element.list_nonzero_elements(p,n)):
                yield finite_sp_matrix([[O,-B.inverse()],[B,A]])
        return itertools.chain(part1(),part2())

    def check_order(self, n):
        if (self ** n).is_identity():
            for x in maximal_proper_factors(n):
                if (self**x).is_identity():
                    return False
            return True
        else:
            return False
    def is_identity(self):
        return self == finite_sp_matrix.identity(len(self.elements), self.p, self.n)
    @classmethod
    def get_element_of_sl_2_from_field_extension(cls,p,n):
        #returns an element of sl_2 with order p**n+1.
        m = finite_field_element([1 if i == 1 else 0 for i in range(2*n)], p, 2*n)
        mat = finite_matrix.from_finite_field_element(m, new_n = n)
        return finite_sp_matrix(mat ** (p**n-1))

    @classmethod
    def list_sp(cls, size, p, n):
        #TODO: finish this method.
        #assert size % 2 == 0
        for A0 in finite_matrix.list_symmetric_matrices(size//2,p,n):
            for C in finite_matrix.list_invertible_matrices(size//2,p,n):
                pass


#_weil_representation_case_2()
# p,n =3,2
# count = 0
# for x in finite_sp_matrix.list_sl_2(p,n):
#     if x.check_order(p**n + 1):
#         assert x.order()==p**n+1
#         count += 1
# print(count)

# for x in finite_sp_matrix.list_sl_2_from_field_extension(3,2):
#     print(x)
# for y in finite_matrix.list_invertible_matrices(4,3,2):
#     print(y)
#     assert not y.determinant().is_zero()

# for y in finite_sp_matrix.list_sl_2(3,1):
#     y.weil_representation()
