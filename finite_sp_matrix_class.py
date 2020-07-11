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
            return upper, lower, finite_matrix.identity(len(self.elements),self.p,self.n)
        #returns

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
                yield finite_matrix([[I,O],[A,I]])*finite_matrix([[B,O],[O,B.inverse()]])*finite_matrix([[I,C],[O,I]])
        def part2():
            for A,B in itertools.product(finite_field_element.list_elements(p,n), finite_field_element.list_nonzero_elements(p,n)):
                yield finite_matrix([[O,-B.inverse()],[B,A]])
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
    print(y)
    assert y.determinant() == finite_field_element.one(3,1)
