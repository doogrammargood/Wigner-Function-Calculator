from affine_plane_class import *
import math
class sl_matrix(object):
    def __init__(self, matrix):
        #assume matrix is a 2x2 nested list of finite_field_elements
        self.matrix = matrix
        self.p = matrix[0][0].p
        self.n = matrix[0][0].n
        assert self.determinant() == finite_field_element.one(self.p,self.n)


    def determinant(self):
        return self.matrix[0][0]*self.matrix[1][1] - self.matrix[1][0]*self.matrix[0][1]

    def inverse(self):
        sl = self.matrix
        a = sl[0][0]
        b = sl[0][1]
        c = sl[1][0]
        d = sl[1][1]
        return sl_matrix([[d, -b],[-c,a]])

    def transpose(self):
        sl = self.matrix
        a = sl[0][0]
        b = sl[0][1]
        c = sl[1][0]
        d = sl[1][1]
        return [[a,c],[b,d]]

    def __mul__(self, other):
        assert self.p == other.p
        assert self.n == other.n
        m1 = self.matrix
        if isinstance(other,sl_matrix):
            m2 = other.matrix
            entry00 = m1[0][0]*m2[0][0] + m1[0][1]*m2[1][0]
            entry01 = m1[0][0]*m2[0][1] + m1[0][1]*m2[1][1]
            entry10 = m1[1][0]*m2[0][0] + m1[1][1]*m2[1][0]
            entry11 = m1[1][0]*m2[0][1] + m1[1][1]*m2[1][1]
            return sl_matrix([[entry00, entry01],[entry10, entry11]])
        elif isinstance(other,point_of_plane):

            x, y = other.x, other.y
            xp = m1[0][0]*x + m1[0][1]*y
            yp = m1[1][0]*x + m1[1][1]*y
            return (point_of_plane( (xp, yp) ))
        elif isinstance(other,line_of_plane):
            pass

    def __eq__(self,other):
        A = self.matrix
        B = other.matrix
        return A[0][0]==B[0][0] and A[1][0]==B[1][0] and A[0][1]==B[0][1] and A[1][1]==B[1][1]

    def __str__(self):
        M = self.matrix
        return str([ [ str(M[0][0]), str(M[0][1]) ], [str(M[1][0]), str(M[1][1]) ]])

    def order(self):
        M = self.matrix
        copy = sl_matrix(M)
        count = 1
        while not copy == sl_matrix.identity(self.p, self.n):
            copy = self*copy
            count +=1
        return count

    def orbit(self,pt):
        def orbit_to_create(self,pt):
            to_yield = point_of_plane(pt.coordinates)
            yield to_yield
            to_yield = self*to_yield
            while to_yield != pt:
                yield to_yield
                to_yield = self*to_yield
        return orbit_to_create(self,pt)

    @classmethod
    def identity(cls, p, n):
        one = finite_field_element.one(p,n)
        zero = finite_field_element.zero(p,n)
        return sl_matrix([[one, zero],[zero, one]])

    @classmethod
    def basis_preserving(cls,a,g):
        a_inv = a.inverse()
        return sl_matrix([[a,finite_field_element.zero(a.p,a.n)],[g,a_inv]])
    @classmethod
    def coset_rep(cls, b, p, n):
        if not isinstance(b,finite_field_element):
            two = finite_field_element.generator(p,n)
            return sl_matrix([[two, finite_field_element.one(p,n)], [-finite_field_element.one(p,n),finite_field_element.zero(p,n)]])
        else:
            return sl_matrix([[finite_field_element.one(p,n),b],[finite_field_element.zero(p,n),finite_field_element.one(p,n)]])
    @classmethod
    def gen_sl(cls,p,n):
        coset_reps = itertools.chain(finite_field_element.list_elements(p,n) , [float('inf')])
        basis_pres = itertools.product(finite_field_element.list_nonzero_elements(p,n), finite_field_element.list_elements(p,n))
        for cr,bp in itertools.product(coset_reps, basis_pres):
            F = sl_matrix.coset_rep(cr, p, n)
            B = sl_matrix.basis_preserving(bp[0],bp[1])
            yield F*B

    @classmethod
    def gen_with_order(cls, p, n, order = 'natural'):
        if order == 'natural':
            order = p**n + 1
        return (s for s in sl_matrix.gen_sl(p, n) if s.order() == order )

def test_gen_sl():
    total_sl = []
    for x in sl_matrix.gen_sl(5,2):
        if x in total_sl:
            print ("Duplicate!")
            print (x)
        total_sl.append(x)
    print("done getting total sl")
    print(len(total_sl))
    for a,b,c,d in itertools.product(finite_field_element.list_elements(5,2), repeat =4):
        if (a*d - b*c).is_one():
            M = sl_matrix([[a,b],[c,d]])
            if not M in total_sl:
                print("Missing M")
                print(M)
    print("All done")

def test_gen_with_order(p,n):
    S = sl_matrix.gen_with_order(p,n)
    one = finite_field_element.one(p,n)
    s=next(S)
    print ("ss")
    R = s.orbit(point_of_plane((one,one)))
    for r in R:
        print (r)

#test_gen_with_order(5,2)

#test_gen_with_order()
#test_gen_sl()
