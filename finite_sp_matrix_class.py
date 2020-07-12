from finite_matrix_class import *
class finite_sp_matrix(finite_matrix):
    def __init__(self, list_representation):
        if isinstance(list_representation, finite_matrix):
            #promote finite matrix to sp finite_matrix
            assert list_representation.is_symplectic()
            self.elements = list_representation.elements
            self.p=list_representation.p
            self.n=list_representation.n
        else:
            super().__init__(list_representation)
            assert self.is_symplectic()

    def __mul__(self, other):
        return finite_sp_matrix(super().__mul__(other))
    def __neg__(self):
        return finite_sp_matrix(super().__neg__())
    def transpose(self):
        return finite_sp_matrix(super().transpose())
    def inverse(self):
        return finite_sp_matrix(super().inverse())

    def block_matrix_decompose(self):
        #returns a list-of-lists of matrices [[A,B][C,D]]
        size = len(self.elements)
        block_size = size//2
        A_elements = [[self.elements[j][i] for i in range(block_size) ]for j in range(block_size)]
        B_elements = [[self.elements[j][i] for i in range(block_size,size)] for j in range(block_size)]
        C_elements = [[self.elements[j][i] for i in range(block_size)] for j in range(block_size,size)]
        D_elements = [[self.elements[j][i] for i in range(block_size,size)] for j in range(block_size,size)]
        return [[finite_matrix(A_elements), finite_matrix(B_elements)],[finite_matrix(C_elements), finite_matrix(D_elements)]]


    def factorize(self):
        #TODO: Is this method totally general? Is it necessary that either C or D is invertible?
        #Maybe this can be deduced from the block-matrix inverse.

        blocks = self.block_matrix_decompose()
        A=blocks[0][0]
        B=blocks[0][1]
        C=blocks[1][0]
        D=blocks[1][1]
        O = finite_matrix.zero(len(self.elements)//2,self.p,self.n)
        if D.determinant().is_zero():
            J = finite_sp_matrix.B_type_matrix(finite_matrix.identity(len(self.elements)//2,self.p,self.n))
            temp = self*J
            if B.determinant().is_zero():
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
        assert blocks[1][0].is_zero()
        M = blocks[0][0]
        A = blocks[0][1]
        assert not M.determinant().is_zero()
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
        assert self == finite_sp_matrix.A_type_matrix(A)
        size = len(blocks[0][0].elements)
        diag = []
        two_inv = finite_field_element([2],self.p,1).inverse()
        for ys in itertools.product(finite_field_element.list_elements(self.p,1), repeat = size):
            col_matrix = finite_matrix([[y] for y in ys])
            row_matrix = col_matrix.transpose()
            result = finite_matrix([[two_inv]]) * (row_matrix * A * col_matrix)
            diag.append(result.character())
        return np.matrix(np.diag(diag))

    def C_type_weil(self):
        blocks = self.block_matrix_decompose()
        assert blocks[0][1] == blocks[1][0]
        C = blocks[0][0]
        C_inv_trans = C.inverse().transpose()
        permutation = [] #the result should be a permutation matrix.
        size = len(blocks[0][0].elements)
        def col_matrix_to_int(col_matrix):
            to_return = 0
            vector = col_matrix.transpose().elements[0]
            vector.reverse() #itertools counts by incrementing the last element most frequently.
            ffe = finite_field_element.from_vector(vector)
            return int(ffe)

        for ys in itertools.product(finite_field_element.list_elements(self.p,1), repeat = size):
            col_matrix = finite_matrix([[y] for y in ys])
            # print(size)
            # print(len(col_matrix.elements))
            # print (len(C_inv_trans.elements[0]))
            result = C_inv_trans * col_matrix
            permutation.append( col_matrix_to_int(col_matrix) )
        result = np.matrix([[1 if i == permutation[x] else 0 for x in range(self.p**size)] for i in range(self.p**size)])
        multiplier = legendre_symbol(int(C.determinant()), self.p)
        return multiplier*result

    def B_type_weil(self):
        #only defined for the symplectic form over a prime field
        assert self == finite_matrix.symplectic_form(len(self.elements),self.p,1)
        neg_two_inv = -finite_field_element([2],self.p,1).inverse()
        to_return = []
        guass_sum = 0
        size = len(self.elements)//2
        #size = len(blocks[0][0].elements)
        for ys in itertools.product(finite_field_element.list_elements(self.p,1), repeat = size):
            current_col = []
            row_y = finite_matrix(ys).transpose()
            guass_sum += (finite_matrix([[neg_two_inv]]) * (row_y * row_y.transpose()) ).character()
            for yps in itertools.product(finite_field_element.list_elements(self.p,1), repeat = size):
                col_yp = finite_matrix(yps)
                val = (finite_matrix([[neg_two_inv]]) * (row_y * col_yp) ).character()
                current_col.append(val)
            to_return.append(current_col)
        to_return = np.matrix(to_return).T
        return guass_sum * to_return

    def weil_representation(self):
        print(self)
        upper, lower, last = self.factorize()
        upper_blocks = upper.block_matrix_decompose()

        up_factor1, up_factor2 = upper.factorize_upper()
        low_factor1, low_factor2, low_factor3 = lower.factorize_lower()
        assert up_factor1 * up_factor2 * low_factor1 * low_factor2 * low_factor3 * last == self
        weil_of_j = finite_sp_matrix(finite_matrix.symplectic_form(len(self.elements),self.p,1)).B_type_weil()
        weil_of_j_inv = np.linalg.inv(weil_of_j)
        if last == finite_matrix.symplectic_form(len(self.elements),self.p,1).inverse():
            return up_factor1.C_type_weil() * up_factor2.A_type_weil() * weil_of_j * low_factor2.A_type_weil() * weil_of_j_inv * weil_of_j_inv
        else:
            return up_factor1.C_type_weil() * up_factor2.A_type_weil() * weil_of_j * low_factor2.A_type_weil() * weil_of_j_inv
        # for f in factors:
        #     print(f)
        print("---")

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
        assert not C.determinant().is_zero()
        size = len(C.elements)
        O = finite_matrix.zero(size, C.p, C.n)
        return finite_sp_matrix([[C,O],[O, C.inverse().transpose()]])

    @classmethod
    def list_sl_2(cls,p,n):
        I = finite_field_element.one(p,n)
        O = finite_field_element.zero(p,n)
        def part1():
            for A,B,C in itertools.product(finite_field_element.list_elements(p,n), finite_field_element.list_nonzero_elements(p,n), finite_field_element.list_elements(p,n)):
                yield finite_sp_matrix([[I,O],[A,I]])*finite_matrix([[B,O],[O,B.inverse()]])*finite_matrix([[I,C],[O,I]])
        def part2():
            for A,B in itertools.product(finite_field_element.list_elements(p,n), finite_field_element.list_nonzero_elements(p,n)):
                yield finite_sp_matrix([[O,-B.inverse()],[B,A]])
        return itertools.chain(part1(),part2())


    @classmethod
    def list_sp(cls, size, p, n):
        #TODO: finish this method.
        assert size % 2 == 0
        for A0 in finite_matrix.list_symmetric_matrices(size//2,p,n):
            for C in finite_matrix.list_invertible_matrices(size//2,p,n):
                pass

# for y in finite_matrix.list_invertible_matrices(4,3,2):
#     print(y)
#     assert not y.determinant().is_zero()

for y in finite_sp_matrix.list_sl_2(3,1):
    y.weil_representation()
