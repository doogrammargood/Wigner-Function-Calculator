import pickle
from finite_field_class import *
from functools import reduce
import numpy as np
#finite_field_element.load_dual_vectors()
def swap_positions(list, pos1, pos2):
    #from Geeks for Geeks
    list[pos1], list[pos2] = list[pos2], list[pos1]
    return list
def find_nonzero(matrix, lower, upper, col):
    #returns the coordinates of the first nonzero entry
    pivot = None
    for y in range(lower, upper):
        if not matrix[y][col].is_zero():
            pivot = y
            break
    return pivot

def scale_row(matrix, row, factor):
    w = len(matrix[0])
    h = len(matrix)
    scaled_row = [factor*matrix[row][c] for c in range(w)]
    return [scaled_row if r ==row else [matrix[r][c] for c in range(w)]for r in range(h)]

def swap_row(matrix,row1, row2):
    w = len(matrix[0])
    h = len(matrix)
    rw1 = [matrix[row1][c] for c in range(w)]
    rw2 = [matrix[row2][c] for c in range(w)]
    return [rw2 if r == row1 else rw1 if r==row2 else [matrix[r][c]for c in range(w)] for r in range(h)]

def subtract_row(matrix, subtractee, subtractor_row,factor):
    #new_row = [piv_row[c]-factor* (subtractor_row[c]) for c in range(w)]
    w = len(matrix[0])
    h = len(matrix)
    new_row = [matrix[subtractee][c] - factor*matrix[subtractor_row][c] for c in range(w)]
    return [[matrix[r][c] for c in range(w)] if not r == subtractee else new_row for r in range(h)]

def clear_column(matrix1, matrix2, y_val, x_val):
    #for use in finding inverse.
    w = len(matrix1[0])
    h = len(matrix1)
    for y in range(h):
        if y != y_val:
            factor = matrix1[y][x_val]
            if not factor.is_zero():
                matrix1 = subtract_row(matrix1, y, y_val, factor)
                matrix2 = subtract_row(matrix2, y, y_val, factor)
    return matrix1, matrix2

def clear_column_down(matrix, y_val, x_val):
    #for use in finding determinant
    w = len(matrix[0])
    h = len(matrix)
    for y in range(y_val+1,h):
        if y != y_val:
            factor = matrix[y][x_val]
            if not factor.is_zero():
                matrix = subtract_row(matrix, y, y_val, factor)
    return matrix


class finite_matrix(object):
    # turns a tuple of finite field elements into another.
    dual_basis_conversion = {}
    def __init__(self, list_representation):
        if isinstance(list_representation, finite_field_element):
            #initialize from a finite field element.
            #defines a vector over the prime field. The vector will become a column matrix
            list_representation = list_representation.to_vector()

        if isinstance(list_representation[0], finite_field_element):
            #convert to column matrix
            self.elements = []
            for item in list_representation:
                self.elements.append([item])
            self.p = list_representation[0].p
            self.n = list_representation[0].n
        elif isinstance(list_representation[0][0], finite_field_element):
            #standard initialization: assume list_representation is a list of lists of finite field elements.
            self.elements = list_representation
            self.p = list_representation[0][0].p
            self.n = list_representation[0][0].n
        elif isinstance(list_representation[0][0], finite_matrix):
            #block matrix initializiation
            #assume all blocks are the have the same types of entries
            self.p = list_representation[0][0].p
            self.n = list_representation[0][0].n

            # for row in range(len(list_representation) - 1):
            #     for col in range(len(list_representation[0]) - 1):
                    #check shape of blocks
                    #assert len(list_representation[row][col].elements[0]) == len(list_representation[row+1][col].elements[0])
                    #assert len(list_representation[row][col].elements) == len(list_representation[row][col+1].elements)
            self.elements = []
            new_row = []
            for row in range(len(list_representation) ):
                for inner_row in range(len(list_representation[row][0].elements)):
                    for col in range(len(list_representation[0])):
                        new_row.extend( list_representation[row][col].elements[inner_row] )
                    self.elements.append(new_row)
                    new_row = []


    def copy(self):
        return finite_matrix(self.elements)

    def transpose(self):
        rows = len(self.elements)
        cols = len(self.elements[0])
        return finite_matrix( [[ self.elements[c][r] for c in range(rows)] for r in range(cols)] )

    def __mul__(self, other):
        w =len(self.elements[0])
        #assert w == len(other.elements)
        to_return = finite_matrix([[ reduce( (lambda x,y: x+y), [self.elements[c][k]*other.elements[k][r] for k in range(w)]) for c in range(len(self.elements))]for r in range(len(other.elements[0]))])
        return to_return.transpose()

    def __neg__(self):
        new_elements = [[-self.elements[r][c] for c in range(len(self.elements[0]))] for r in range(len(self.elements))]
        return finite_matrix(new_elements)

    def order(self):
        counter = 1
        curr = self.copy()
        #assert len(self.elements)==len(self.elements[0])
        I = finite_matrix.identity(len(self.elements),self.p,self.n)
        while(not curr == I):
            curr = curr*self
            counter += 1
        return counter

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
        #assert len(self.elements[0])==1
        return self.transpose().elements[0]

    def is_symplectic(self):
        if self.transpose()*finite_matrix.symplectic_form(len(self.elements),self.p,self.n)*self == finite_matrix.symplectic_form(len(self.elements),self.p,self.n):
            return True
        else:
            return False

    def is_zero(self):
        #checks to see if all elements of the matrix are zero
        for row in self.elements:
            for item in row:
                if not item.is_zero():
                    return False
        return True

    def to_prime_field_matrix(self):
        return finite_matrix([ [finite_matrix.from_finite_field_element(self.elements[i][j]) for j in range(len(self.elements[0]))] for i in range(len(self.elements))])


    def character(self):
        #from https://arxiv.org/pdf/1906.07230.pdf
        #assert len(self.elements)==len(self.elements[0])==1
        #Assumes the matrix is 1x1
        lamb = self.elements[0][0]
        #assert lamb.n == 1
        p = lamb.p
        return np.e**(1j*2*np.pi*int(lamb.trace())/p)
    @classmethod
    def identity(cls,size,p,n):
        return finite_matrix( [[ finite_field_element.one(p,n) if c ==r else finite_field_element.zero(p,n) for c in range(size)]for r in range(size)] )

    @classmethod
    def zero(cls, size, p, n):
        return (finite_matrix([[finite_field_element.zero(p,n) for c in range(size)] for r in range(size)]))

    @classmethod
    def symplectic_form(cls,size,p,n):
        #assert size%2==0
        I = finite_matrix.identity(size//2,p,n)
        O = finite_matrix.zero(size//2,p,n)
        to_return = finite_matrix([[O,I],[-I,O]])
        return to_return

    @classmethod
    def from_finite_field_element(cls, ffe):
        #defines a square matrix over the prime field.
        #As opposed to the constructor method in __init__ which creates a column matrix
        n=ffe.n
        p=ffe.p
        basis_elements = [finite_field_element([1 if i ==j else 0 for i in range(n)],p,n) for j in range(n)]
        vects = []
        for b in basis_elements:
            vects.append( (ffe*b).to_vector() )
        M = finite_matrix(vects)
        return M.transpose()

    @classmethod
    def list_square_matrices(cls, m,p,n):
        for elements in itertools.product(finite_field_element.list_elements(p,n), repeat = m**2):
            yield finite_matrix( [ [elements[i*m + j] for i in range(m)] for j in range(m)] )

    @classmethod
    def list_invertible_matrices(cls, m,p,n):

        #Uses partial Guassian Elimination as suggested here: https://stackoverflow.com/questions/62766721/how-to-list-all-invertible-matrices-over-a-finite-field/62767443#62767443

        current_elements = []
        list_of_row_iterators = [ itertools.chain(itertools.product(finite_field_element.list_elements(p,n), repeat = m), [None]) for r in range(m)]

        current_index = 0
        completed = False

        while not completed:
            new_row = next(list_of_row_iterators[current_index])
            if not new_row is None:
                if current_index >= len(current_elements):
                    current_elements.append(new_row)
                else:
                    current_elements[current_index]=new_row

                reduced_list = finite_matrix([c for c in current_elements[:current_index+1]]).row_eschelon().elements
                if not all([e.is_zero() for e in reduced_list[-1]]):
                    current_index += 1
                    if current_index==m:
                        current_index = m-1
                        yield finite_matrix(current_elements)


            else:
                list_of_row_iterators[current_index] = itertools.chain(itertools.product(finite_field_element.list_elements(p,n), repeat = m), [None])
                temp_elements = [c for e,c in enumerate(current_elements) if (e != current_index) and (e != current_index-1)]
                current_index = current_index - 1
                if current_index == -1:
                    completed = True




    @classmethod
    def list_symmetric_matrices(cls,m,p,n):
        for elements in itertools.product(finite_field_element.list_elements(p,n), repeat = m*(m+1)//2):
            yield finite_matrix( [ [elements[i*(i+1)//2 + j] if i>=j else elements[j*(j+1)//2 + i]  for i in range(m)] for j in range(m)] )

    @classmethod
    def load_dual_basis_matrices(cls):
        file = open("dual_basis_transformation_matrices.pickle", "rb")
        cls.dual_basis_conversion = pickle.load(file)

    @classmethod
    def convert_to_dual(cls, ffe):
        # if ffe.n == 1:
        #     return ffe
        col_mat = cls.dual_basis_conversion[(ffe.p,ffe.n)] * finite_matrix(ffe)
        return finite_field_element.from_vector(col_mat.column_matrix_to_vector())

    #Guassian Elimination Methods:


    def inverse(self):
        elts = self.elements.copy()
        pivot = None
        h = len(elts)
        w = len(elts[0])
        #assert w==h
        size = w
        iden = finite_matrix.identity(size,self.p,self.n).elements
        x_index = 0 #No columns have been cleared
        index_y = 0 #the currently selected index
        top = 0

        for x_index in range(w):
            #print(elts)
            #print("ooh")
            pivot = find_nonzero(elts,top,h,x_index)

            elts = swap_row(elts, top, pivot)
            iden = swap_row(iden, top, pivot)
            scale = elts[top][x_index].inverse()
            elts = scale_row(elts, top, scale )
            iden = scale_row(iden, top, scale )

            elts, iden = clear_column(elts, iden, top, x_index)
            top += 1
        #assert finite_matrix(elts) == finite_matrix.identity(size, self.p,self.n)
        return finite_matrix(iden)

    def determinant(self):
        #computes a row eschelon form.
        elts = self.elements.copy()
        pivot = None
        size = len(elts)
        #top = 0
        det = finite_field_element.one(self.p,self.n)
        #print(size)
        #assert len(elts[0])==size #determinant is only defined for square matrices
        for x_index in range(size):
            pivot = find_nonzero(elts,x_index,size,x_index)
            if pivot is None:
                return finite_field_element.zero(self.p,self.n)
            if x_index != pivot:
                det = -det
                elts = swap_row(elts, x_index, pivot)
            scale = elts[x_index][x_index]
            det = det*scale
            elts = scale_row(elts, x_index, scale.inverse() )
            elts = clear_column_down(elts, x_index, x_index)
        return det

    def row_eschelon(self, permutation = False):
        #returns the row-eschelon form of the matrix.
        #If permutation == True, then it will also return the permutation from the input matrix to its rre form.
        elts = self.elements.copy()
        pivot = None
        size = len(elts)
        #assert len(elts[0])>=size
        current_row = 0
        permute = [i for i in range(len(elts))] #initially, the identity.
        for x_index in range(len(elts[0])):
            pivot = find_nonzero(elts,current_row,size,x_index)
            if pivot is None:
                continue
                #return finite_field_element.zero(self.p,self.n)
            if current_row != pivot:
                elts = swap_row(elts, current_row, pivot)

                if permutation:
                    swap_positions(permute, current_row, pivot)


            scale = elts[current_row][x_index]
            elts = scale_row(elts, current_row, scale.inverse() )
            elts = clear_column_down(elts, current_row, x_index)
            current_row +=1
        if permutation:
            #invert the permutation:
            #permute =[x[0] for x in sorted(enumerate(permute), key = lambda x: x[1])]
            return finite_matrix(elts), permute
        else:
            return finite_matrix(elts)
# one = finite_field_element.one(5,2)
# zero = finite_field_element.zero(5,2)
# M = finite_matrix([ [one, zero, one],
#                     [zero, zero, one],
#                     [zero, one, one]])

def process_dual_matrices(dict, p, n):
    vectors = []
    for i in range(n):
        vectors.append(dict[(i,p,n)].to_vector())
    #print(M)
    M = finite_matrix(vectors)
    M = M.transpose()
    M = M.inverse()
    assert M*(M.inverse() )==finite_matrix.identity(n,p,1)
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
    #file = open("dual_basis_transformation_matrices.pickle", "wb")
    #pickle.dump(to_store, file)

#record_dual_basis_transform()
#finite_matrix.load_dual_basis_matrices()
def check_dual_basis_transform():
    finite_matrix.load_dual_basis_matrices()
    x = finite_field_element([1,0,3], 5,3)
    y = finite_matrix.convert_to_dual(x)
    print(y)
    for i in range(3):
        print(finite_field_element.dual_basis[(i,3,3)])

# def check_block_matrix_initialization():
#     I = finite_matrix.identity(3,3,1)
#     O = finite_matrix.zero(3,3,1)
#     to_check = finite_matrix([[O,I],[I,O]])
#     print(to_check)
#check_block_matrix_initialization()
#check_dual_basis_transform()
#record_dual_basis_transform()
# a = finite_field_element([2,1],3,2)
# i = finite_field_element.one(3,2)
# o = finite_field_element.zero(3,2)
# A = finite_matrix([[a,a,a],[a,a,a],[a,a,a]])
# print(A.determinant())
#for m in [ m for m in finite_matrix.list_symmetric_matrices(2,3,1) if not m.determinant().is_zero()]:
#for m in finite_matrix.list_symmetric_matrices(2,3,1):
# for m in finite_matrix.list_square_matrices(2,3,1):
#     if m.is_symplectic():
#         print(m)
#         print(m.determinant())
    # print(m)
    # print(m.determinant())
# print("---")
# for m in [ m for m in finite_matrix.list_square_matrices(2,3,1) if not m.determinant().is_zero()]:
#     print(m)
