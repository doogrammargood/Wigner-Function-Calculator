#from sage.all import *
import itertools


def modinv(a, m):
    a = a % m;
    for x in range(1, m) :
        if ((a * x) % m == 1) :
            return x
    return 1

class polynomial_element(object):
    #n is an upper bound on the degree, p is the prime field of coefficients.
    def __init__(self,coordinates,p,n):
        self.p=p
        self.n=n
        self.coordinates = [c%p for c in coordinates] #coordinate list should have length n

    def __add__(self, other):
        new_n = max(self.n, other.n)
        temp1 = [self.coordinates[index] if index < len(self.coordinates) else 0 for index in range(new_n)]
        temp2 = [other.coordinates[index] if index < len(other.coordinates) else 0 for index in range(new_n)]
        new_coordinates = [(a+b)%self.p for a,b in zip(temp1, temp2)]
        return polynomial_element(new_coordinates, self.p, max(self.n, other.n))

    def copy(self):
        return polynomial_element(self.coordinates, self.p, self.n)

    def degree(self):
        temp = [e for e in enumerate(self.coordinates) if e[1] != 0 ]
        if len(temp)==0:
            return 0
        else: return max(temp, key= lambda x: x[0])[0]

    def trim_to(self,n):
        self.n = n
        self.coordinates = [self.coordinates[index] for index in range(n)]

    def is_zero(self):
        return all([c==0 for c in self.coordinates])

    def _add_inverse(self):
        neg_coords = [(-c)%self.p for c in self.coordinates]
        return polynomial_element(neg_coords, self.p, self.n)

    def __sub__(self,other):
        return self.__add__(other._add_inverse())

    def __eq__(self,other):
        if other is None:
            return False
        return all([c==d for c,d in zip(self.coordinates, other.coordinates)])

    def __mul__(self, other):
        list1 = self.coordinates
        list2 = other.coordinates
        new_length = len(list1) + len(list2) - 1
        new_list = [sum([list1[index]*list2[new_dex-index] for index in range(len(list1)) if index <=new_dex and index + len(list2) > new_dex]) for new_dex in range(new_length) ]
        to_return = polynomial_element(new_list, self.p, self.n + other.n -1)
        return to_return

    def leading_term(self):
        #returns (a,b), where b x^a is the leading term.
        terms = [e for e in enumerate(self.coordinates) if e[1] != 0]
        if len(terms)==0:
            return((0,0))
        else:
            return max(terms, key= lambda x: x[0])
    def remove_leading_term(self):
        pass

    def __str__(self):
        return str(self.coordinates)

    def __truediv__(self, other):
        #implements the division algorithm self = a*other + r; returns (a,r)
        self_copy = polynomial_element(self.coordinates, self.p, self.n)
        degree_diff = self_copy.degree()-other.degree()
        quotient = polynomial_element([0]*(degree_diff+1), self.p, degree_diff + 1)
        if other.is_zero():
            print('DIVISION BY ZERO')
            return None

        def reduce(self_copy, other, quotient):

            leading_self = self_copy.leading_term()
            leading_other = other.leading_term()
            shift = leading_self[0] - leading_other[0]
            mul = (leading_self[1]*modinv(leading_other[1],self.p) )%self.p
            pol = polynomial_element([mul if index == shift else 0 for index in range(quotient.n)],quotient.p,quotient.n)
            quotient = quotient + pol
            self_copy = self - other*quotient
            return self_copy, quotient

        while degree_diff >=0 and not self_copy.is_zero():
            self_copy, quotient=reduce(self_copy, other, quotient)
            degree_diff = self_copy.degree()-other.degree()
        self_copy.trim_to(self.n)
        return (quotient, self_copy)

    def euclidean_algorithm(self,other):
        #Stright from the Wvikerpedias
        assert self.p == other.p
        r0 = self.copy()
        r1 = other.copy()
        one = polynomial_element([1 if i ==0 else 0 for i in range(self.n)], self.p, self.n)
        zero = polynomial_element([0]*self.n, self.p, self.n)
        s0, s = one, zero
        t0, t = zero, one
        while (not r1.is_zero()):
            temp_r = r1.copy()
            temp_s = s.copy()
            temp_t = t.copy()
            q, r1 = r0/r1
            s = s0 - q*s
            t = t0 - q*t
            r0 = temp_r
            s0 = temp_s
            t0 = temp_t
        return r0, s0, t0


class finite_field_element(polynomial_element):
    #These are polynomial_elements identified up to their conway polynomials.
    conway_polynomial={(3,3): polynomial_element([1,2,0,1],3,4)}

    @classmethod
    def load_conway_polynomials(cls):
        f = open('conway_poly_lists.txt', 'r')
        lines = f.readlines()
        for l in lines:
            ln = l.split("%")
            #print(eval(ln[0]))
            key=(eval(ln[0]), eval(ln[1]))
            #print (key)
            list = eval( ln[2] )

            value = polynomial_element(list,key[0],key[1])
            cls.conway_polynomial[key]=value
        #print(finite_field_element.conway_polynomial)

    def as_poly(self):
        return polynomial_element(self.coordinates, self.p, self.n)

    def is_zero(self):
        return self.as_poly().is_zero()

    def is_one(self):
        return self == finite_field_element.one(self.p,self.n)

    def __init__(self,*args, **kwargs):
        if len(args)==3:
            coordinates, p, n= args
            super().__init__(coordinates,p,n)
            conway_polynomial = finite_field_element.conway_polynomial[(self.p,self.n)]
            new_poly = (self.as_poly() / conway_polynomial)[1]
            self.coordinates = new_poly.coordinates

        elif len(args)==1:
            poly = args[0]
            conway_polynomial = finite_field_element.conway_polynomial[(poly.p,poly.n)]
            self.coordinates = ((poly / conway_polynomial)[1] ).coordinates
            self.p = poly.p
            self.n = poly.n


    def __mul__(self,other):
        temp = self.as_poly() * other
        temp.n = self.n
        return finite_field_element(temp)

    def __add__(self, other):
        return finite_field_element(self.as_poly() + other)

    def __sub__(self,other):
        return finite_field_element(self.as_poly() - other)

    def __neg__(self):
        return finite_field_element.zero(self.p,self.n)-self

    def inverse(self):
        if self.degree() == 0:
            x = modinv(self.coordinates[0], self.p)
            return finite_field_element([x if i==0 else 0 for i in range(self.n)], self.p, self.n)

        r, s, t = self.as_poly().euclidean_algorithm(finite_field_element.conway_polynomial[(self.p,self.n)])
        #r = finite_field_element(r.coordinates,self.p, self.n)
        assert r.degree()==0
        assert not r.is_zero()
        x = r.coordinates[0]
        x = modinv(x,self.p)

        s = finite_field_element([i*x for i in s.coordinates],self.p, self.n)
        return s

    def __truediv__(self,other):
        return self * other.inverse()

    def __str__(self):
        return str(self.coordinates)

    def __int__(self):
        cur_sum = 0
        for i in range(self.n):
            cur_sum += self.coordinates[i]*self.p**i
        return cur_sum

    def order(self):
        if self.is_one():
            return 1

        x = self.copy()
        x = self * x
        count = 2
        while not x == self:
            x = x*self
            count += 1
        return count
    @classmethod
    def list_elements(cls,prime, power):
        return ( (finite_field_element(list(reversed(l)),prime,power) for l in itertools.product(range(prime), repeat = power) ) )

    @classmethod
    def list_nonzero_elements(cls,prime, power):
        return ( (finite_field_element(list(reversed(l)),prime,power) for l in itertools.product(range(prime), repeat = power) if not finite_field_element(l,prime,power).is_zero() ) )
    @classmethod
    def one(cls, p, n):
        return (finite_field_element([1 if i == 0 else 0 for i in range(n)], p, n))
    @classmethod
    def zero(cls, p, n):
        return (finite_field_element([0]*n, p, n))

    @classmethod
    def generator(cls,p,n):
        #returns a generator.
        for x in finite_field_element.list_elements(p,n):
            if x.order() == p**n:
                return x

def test_finite_fields():
    #Tests inverses, associativity of addition, distributivity, commutativity.
    finite_field_element.load_conway_polynomials()
    test_primes = [3,5]
    test_powers = [1,2]

    for prime, power in itertools.product(test_primes, test_powers):
        lister = finite_field_element.list_elements(prime,power)
        for e1, e2 in itertools.product(lister, repeat =2):
            if not e1.is_zero():
                assert e1*e1.inverse()*e2 == e2
                assert e1*e2 == e2*e1
        lister = finite_field_element.list_elements(prime,power)
        for e1, e2, e3 in itertools.product(lister, repeat = 3):
            assert (e1+e2)+e3==e1+(e2+e3)
            assert (e1+e2)*e3==e1*e3 + e2*e3

#test_finite_fields()

finite_field_element.load_conway_polynomials()
#print (finite_field_element.generator(5,1))
# def test_division():
#
#     primes = [3,5,7]
#     ns = [2,3]
#     for p in primes:
#         for n in ns:
#             for list1, list2 in itertools.combinations(itertools.product(range(p), repeat = n),2):
#                 a = polynomial_element(list1, p, n)
#                 b = polynomial_element(list2, p, n)
#                 q,r = a/b
#
#                 if not b.is_zero():
#                     assert b*q+r==a
#                     print("good")

#test_division()
# print(finite_field_element.conway_polynomial[(5,2)])
#x = finite_field_element([1,1,1],3,3)
#y = finite_field_element([1,2,0],3,3)


#r, s, t = y.euclidean_algorithm(x)
# print(s)
# print (t)
# print(r)
# print(x)
# print(y)
#assert r == s*y + t*x
# a = finite_field_element(x)
# a = finite_field_element([1,1],3,2)
# b = finite_field_element([2,1],3,2)
#print (a- a._add_inverse()).coordinates
#print (a+b).coordinates
#print (a*b).coordinates


#rint a._polynomial_remainder([3,3,1],[2,1,0],3)
