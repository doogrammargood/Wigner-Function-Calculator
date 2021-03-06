from finite_field_class import *
from object_modified_class import *
class point_of_plane(object_modified):
    def __init__(self, coordinates):
        #Expect coordinates to be a pair (x,y) where x and y are finite_field_elements.
        if coordinates == None:
            self.coordinates = (None,None)
            self.p = None
            self.n = None
            self.x = None
            self.y = None
        else:
            assert coordinates[0].p == coordinates[1].p
            assert coordinates[0].n == coordinates[1].n
            self.coordinates = coordinates
            self.p = coordinates[0].p
            self.n = coordinates[0].n
            self.x = self.coordinates[0]
            self.y = self.coordinates[1]
    def line_to(self, other):
        #returns the line from self to other.
        if other.isNone() or other == self:
            return line_of_plane(None)
        dx = self.x-other.x
        dy = self.y-other.y
        if not dy.is_zero():
            a = finite_field_element.one(self.p, self.n)
            b = -dx/dy
            return line_of_plane((a, b, a*self.x + b*self.y))
        else:
            a = finite_field_element.zero(self.p,self.n)
            b = finite_field_element.one(self.p,self.n)
            c = dy
            return line_of_plane((a, b, a*self.x + b*self.y))
    def __str__(self):
        return(str(self.x)+", "+str(self.y))
        #return(str(self.coordinates))

    def __eq__(self, other):
        if other.isNone():
            return False
        return self.x == other.x and self.y == other.y

    def isNone(self):
        if self.p is None:
            return True
        else:
            return False

    def is_on_line(self,line):
        if self.isNone() or line.isNone():
            return False
        a,b,c = line.coefficients
        return ( (a*self.x+b*self.y-c).is_zero() )

    def gen_lines(self):
        if self.isNone():
            return([])
        #returns a generator for the lines through the point.
        iter1 = ( line_of_plane( (finite_field_element.one(self.p, self.n), b, self.x +b*self.y) ) for b in finite_field_element.list_elements(self.p, self.n) )
        iter2 = ( line_of_plane( (finite_field_element.zero(self.p, self.n),  finite_field_element.one(self.p, self.n), self.y) ) for i in range(1))
        return itertools.chain( iter1, iter2 )

    @classmethod
    def from_string(cls, string,p,n):
        if string[0] == 'N':
            return point_of_plane(None)
        coords = string.split("], ")
        x = eval(coords[0] + "]")
        y = eval(coords[1])
        x = finite_field_element(x,p,n)
        y = finite_field_element(y,p,n)
        return point_of_plane((x,y))

    @classmethod
    def origin(cls,p,n):
        return point_of_plane((finite_field_element.zero(p,n),finite_field_element.zero(p,n)))


class line_of_plane(object_modified):
    def __init__(self, coefficients):
        if coefficients is None:
            self.coefficients = [None, None, None]
            self.p = None
            self.n = None
        else:
            a,b,c = coefficients
            assert a.p==b.p==c.p
            assert a.n==b.n==c.n
            self.coefficients = (a,b,c)
            self.p = a.p
            self.n = a.n
            self.coefficients = self.canonize_line()

    def __eq__(self,other):
        if all([self.coefficients[i]==other.coefficients[i] for i in range(3)]):
            return True
        else:
            return False
    def __neg__(self):
        return line_of_plane((self.coefficients[0], self.coefficients[1], -self.coefficients[2]))
    def canonize_line(self):
        a, b, c, = self.coefficients
        if not a.is_zero():
            ainv = a.inverse()
            a, b, c = finite_field_element.one(a.p, a.n), ainv*b, ainv*c
        else:
            binv = b.inverse()
            b, c = finite_field_element.one(a.p,a.n), binv*c
        return (a,b,c)

    def parallel_through(self,pt):
        #returns the line parallel to self through pt.
        a, b, c, = self.coefficients
        c = a*pt.x + b*pt.y
        return line_of_plane((a,b,c))

    def __str__(self):
        a,b,c = self.coefficients
        return str(a) + str(b) + str(c)
    @classmethod
    def from_string(cls,string,p,n):
        str_coefficients = string.split(']')[:-1]
        str_coefficients = [s + ']' for s in str_coefficients]
        coefficients = [ finite_field_element(eval(s), p,n) for s in str_coefficients]
        if string[0]=='N':
            return line_of_plane(None)
        else:
            str1, str2, str3 = str_coefficients
            return line_of_plane(coefficients)

    def gen_points(self):
        #returns a generator for the points on the line
        if self.isNone():
            return([])
        a,b,c = self.coefficients
        field_elements = finite_field_element.list_elements(self.p,self.n)
        if not b.is_zero():
            for x in field_elements:
                yield( point_of_plane( (x, (c-(a*x)) /b ) ) )
        else:
            for y in field_elements:
                yield( point_of_plane((c/a, y) )  )
            #ax+by=c => y=(c-ax)/b
    def gen_parallel_lines(self):
        if self.isNone():
            return([])
        a,b,c = self.coefficients
        return ( (line_of_plane((a,b,x)) for x in finite_field_element.list_elements(self.p, self.n)) )

    def is_parallel(self, other):
        a,b,c = self.coefficients
        d,e,f = other.coefficients
        if a==d and b==e:
            return True
        else:
            return False

    def intersect_line(self, other):
        if self.is_parallel(other):
            return None #Note weird behavior for intersecting with self.
        else:
            a,b,c = self.coefficients
            d,e,f = other.coefficients
            if not b.is_zero():
                x= (f - ((e/b)*c) ) / (d - a*(e/b))
                y = (c-a*x)/b
            else:
                a,b,c = other.coefficients
                d,e,f = self.coefficients
                x= (f - ((e/b)*c) ) / (d - a*(e/b))
                y = (c-a*x)/b
            return point_of_plane((x,y))


def test_line_to():
    #o = finite_field_element.zero(5,1)
    #origin = point_of_plane((o,o))
    x1 =finite_field_element([4],5,1)
    y1 = finite_field_element([0],5,1)
    x2 = finite_field_element([2],5,1)
    y2 = finite_field_element([1],5,1)
    pt1 = point_of_plane((x1,y1))
    pt2 = point_of_plane((x2,y2))

    l = pt1.line_to(pt2)
    assert pt1.is_on_line(l)
    assert pt2.is_on_line(l)
    count = 0



# x = finite_field_element([2,1],3,2)
# p = point_of_plane((x,x))
# #o = finite_field_element.zero
# o =point_of_plane( (finite_field_element.zero(3,2), finite_field_element.zero(3,2) ) )
# l = p.line_to(o)
# #print(o)
# count = 0
# for line1, line2 in itertools.product(p.gen_lines(), repeat =2):
#     if line1.intersect_line(line2) is not None:
#         assert p==line1.intersect_line(line2)
    #print(a)
    #pass
#print(p.return_line())
