import pickle
from finite_field_class import *
from functools import reduce
#finite_field_element.load_dual_vectors()
class finite_matrix(object):
    # turns a tuple of finite field elements into another.
    dual_basis_conversion = {}
    def __init__(self, list_representation):
        #assume list_representation is a list of lists of finite field elements.
        if isinstance(list_representation, finite_field_element):
            list_representation = list_representation.to_vector()
        if isinstance(list_representation[0], finite_field_element):
            #convert to column matrix
            self.elements = []
            for item in list_representation:
                self.elements.append([item])
            self.p = list_representation[0].p
            self.n = list_representation[0].p
        else:
            self.elements = list_representation
            self.p = list_representation[0][0].p
            self.n = list_representation[0][0].n

    def copy(self):
        return finite_matrix(self.elements)

    def transpose(self):
        rows = len(self.elements)
        cols = len(self.elements[0])
        return finite_matrix( [[ self.elements[c][r] for c in range(rows)] for r in range(cols)] )

    def __mul__(self, other):
        w =len(self.elements[0])
        assert w == len(other.elements)
        to_return = finite_matrix([[ reduce( (lambda x,y: x+y), [self.elements[c][k]*other.elements[k][r] for k in range(w)]) for c in range(len(self.elements))]for r in range(len(other.elements[0]))])
        return to_return.transpose()

    def __eq__(self, other):
        if not self.p == other.p:
            return False
        elif not self.n == other.n:
            return False
        elif not len(self.elements) == len(other.elements):
            return False
        elif not len(self.elements[0])==len(other.elements[0]):
            return False
        for i in range(len(self.elements)):
            for j in range(len(self.elements[0])):
                if not self.elements[i][j] == other.elements[i][j]:
                    return False
        return True

    def __str__(self):
        s = "["
        for row in range(len(self.elements)):
            s+="["
            for col in range(len(self.elements[0])):
                s += (str(self.elements[row][col]) + ",")
            s+=("],")
            s+="\n"
        s+=("]")
        return s

    def column_matrix_to_vector(self):
        assert len(self.elements[0])==1
        return self.transpose().elements[0]
    @classmethod
    def identity(cls,size,p,n):
        return finite_matrix( [[ finite_field_element.one(p,n) if c ==r else finite_field_element.zero(p,n) for c in range(size)]for r in range(size)] )

    @classmethod
    def load_dual_basis_matrices(cls):
        file = open("dual_basis_transformation_matrices.pickle", "rb")
        cls.dual_basis_conversion = pickle.load(file)

    @classmethod
    def convert_to_dual(cls, ffe):
        col_mat = cls.dual_basis_conversion[(ffe.p,ffe.n)] * finite_matrix(ffe)
        #print(finite_matrix(ffe))
        #print(col_mat)
        return finite_field_element.from_vector(col_mat.column_matrix_to_vector())
    def inverse(self):
        elts = self.elements.copy()
        pivot = None
        h = len(elts)
        w = len(elts[0])
        assert w==h
        size = w
        iden = finite_matrix.identity(size,self.p,self.n).elements
        top = 0 #No columns have been cleared
        index_y = 0 #the currently selected index
        def find_nonzero(matrix, lower, upper, col):
            #returns the coordinates of the first nonzero entry
            pivot = None
            for y in range(lower, upper):
                if not matrix[y][col].is_zero():
                    pivot = y
                    break
            return pivot

        def scale_row(matrix, row, factor):
            scaled_row = [factor*matrix[row][c] for c in range(w)]
            return [scaled_row if r ==row else [matrix[r][c] for c in range(w)]for r in range(h)]

        def swap_row(matrix,row1, row2):
            rw1 = [matrix[row1][c] for c in range(w)]
            rw2 = [matrix[row2][c] for c in range(w)]
            return [rw2 if r == row1 else rw1 if r==row2 else [matrix[r][c]for c in range(w)] for r in range(h)]

        def subtract_row(matrix, subtractee, subtractor_row,factor):
            #new_row = [piv_row[c]-factor* (subtractor_row[c]) for c in range(w)]
            new_row = [matrix[subtractee][c] - factor*matrix[subtractor_row][c] for c in range(w)]
            return [[matrix[r][c] for c in range(w)] if not r == subtractee else new_row for r in range(h)]

        def clear_column(matrix1, matrix2, y_val, x_val):
            for y in range(h):
                if y != y_val:
                    factor = matrix1[y][x_val]
                    if not factor.is_zero():
                        matrix1 = subtract_row(matrix1, y, y_val, factor)
                        matrix2 = subtract_row(matrix2, y, y_val, factor)
            return matrix1, matrix2

        for x_index in range(w):
            #print(elts)
            #print("ooh")
            pivot = find_nonzero(elts,top,h,x_index)
            elts = swap_row(elts, top, pivot)
            iden = swap_row(iden, top, pivot)

            elts = scale_row(elts, top, elts[top][x_index].inverse() )
            iden = scale_row(iden, top, elts[top][x_index].inverse() )

            elts, iden = clear_column(elts, iden, top, x_index)
            top += 1
        assert finite_matrix(elts) == finite_matrix.identity(size, self.p,self.n)
        return finite_matrix(iden)
# one = finite_field_element.one(5,2)
# zero = finite_field_element.zero(5,2)
# M = finite_matrix([ [one, zero, one],
#                     [zero, zero, one],
#                     [zero, one, one]])

def process_dual_matrices(dict, p, n):
    vectors = []
    for i in range(n):
        vectors.append(dict[(i,p,n)].to_vector())
    M = finite_matrix(vectors)
    M = M.transpose()
    M = M.inverse()
    return M
def record_dual_basis_transform():
    dict = finite_field_element.dual_basis
    to_store={}
    for p in [p for p in range(3,50) if is_prime(p)]:
        if p < 10:
            for n in range(2,5):
                to_store[(p,n)]=process_dual_matrices(dict,p,n)
        elif p < 20:
            for n in range(2,4):
                to_store[(p,n)]=process_dual_matrices(dict,p,n)
        else:
            for n in range(2,3):
                to_store[(p,n)]=process_dual_matrices(dict,p,n)
    file = open("dual_basis_transformation_matrices.pickle", "wb")
    pickle.dump(to_store, file)
#record_dual_basis_transform()
#finite_matrix.load_dual_basis_matrices()
def check_dual_basis_transform():
    x = finite_field_element([1,0,3], 5,2)
    print(finite_matrix.convert_to_dual(x))
#check_dual_basis_transform()
