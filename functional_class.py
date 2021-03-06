from grid_class import *
from finite_sp_matrix_class import *
import profile


#This file defines the class of functionals on grids. These are a collection of lines on the grid to be counted.
#For multidimensional particles, it should really be a collection of isotropic spaces.
def top_lines_from_direction(grid, direction, percent, from_sl = False):
    #I now realize this method for generating functionals is doomed to fail, due to the entropic uncertainty principle.
    p,n = direction.coefficients[0].p, direction.coefficients[0].n
    percent_copy = percent
    a,b = direction.coefficients[0], direction.coefficients[1]
    marginal = grid.marginalize_grid(direction)
    if from_sl:
        marginal = enumerate(marginal)
        representative_marginal = []
        for c, value in marginal:
            if (not c in [x[0] for x in representative_marginal]) and (c not in [int(-finite_field_element.from_int(x[0],p,n)) for x in representative_marginal]):
                if c!= 0:
                    representative_marginal.append((c,value))
        marginal = [x[1] for x in representative_marginal]
        percent_copy = percent / 2

    sorted_decorated_marginal = sorted(enumerate(marginal), key = lambda x: x[1], reverse = True) #sort decreasing by value
    small_sum = 0
    lines_in_this_direction = []
    current_index = 0
    while small_sum < percent_copy and current_index < len(sorted_decorated_marginal):
        lines_in_this_direction.append(sorted_decorated_marginal[current_index][0])
        small_sum += sorted_decorated_marginal[current_index][1]
        current_index += 1
    if from_sl:
        lines_in_this_direction = [line_of_plane((a,b, finite_field_element.from_int(representative_marginal[index][0],p,n))) for index in lines_in_this_direction] + [line_of_plane((a,b, -finite_field_element.from_int(representative_marginal[index][0],p,n))) for index in lines_in_this_direction]
        lines_in_this_direction.sort(key = lambda x: int(x.coefficients[2]))
    else:
        lines_in_this_direction = [line_of_plane((a,b, finite_field_element.from_int(index,p,n))) for index in lines_in_this_direction]


    return lines_in_this_direction

class functional_on_grid(object_modified):
    def __init__(self, dict_of_lines, p, n, multiparticle = False):
        #assume dict of lines is a dictonary from str(direction) to lists of  pairs: (lines, weight) where weights are real.
        self.dictionary_of_lines = {}
        self.point_count = None #counts the number of negative points.
        self.p = p
        self.n = n
        for direction in point_of_plane.origin(p,n).gen_lines():
            if str(direction) in dict_of_lines.keys():
                if isinstance(dict_of_lines[str(direction)][0],tuple):
                    self.dictionary_of_lines[str(direction)] = dict_of_lines[str(direction)]
                else:
                    self.dictionary_of_lines[str(direction)] = zip(dict_of_lines[str(direction)], [1]*len(dict_of_lines[str(direction)]))
            else:
                self.dictionary_of_lines[str(direction)] = []
    def list_lines(self):
        for key in self.dictionary_of_lines.keys():
            for item in self.dictionary_of_lines[key]:
                yield item

    def __mul__(self, other):
        assert isinstance(other, grid_element)
        return self.evaluate(other)

    def evaluate(self, grid):
        #evaluates the functional on the grid.
        assert grid.p == self.p
        assert grid.n == self.n
        sum = 0
        for direction in point_of_plane.origin(self.p,self.n).gen_lines():
            for item in self.dictionary_of_lines[str(direction)]:
                #sum += grid.marginalize_grid(direction)[int(item.coefficients[2])]
                sum += item[1]*grid.sum_line(item[0])
        return sum

    def evaluate_classical_pt(self, grid, pt):
        counter = 0
        for l in self.list_lines():
            if pt.is_on_line(l[0]):
                counter += l[1]
        return counter

    def evaluate_classical(self,grid):
        #finds the point in the grid which lies on the most lines.
        #I wish we could overload '@' for this.
        assert grid.p == self.p
        assert grid.n == self.n
        best_val = -1
        for x in finite_field_element.list_elements(self.p,self.n):
            for y in finite_field_element.list_elements(self.p,self.n):
                pt = point_of_plane((x,y))
                new_val = self.evaluate_classical_pt(grid, pt)
                if new_val > best_val:
                    best_pt = pt
                    best_val = new_val
        return pt, best_val

    @classmethod
    def top_outcomes(cls, grid, percent):
        #returns a functional which counts the most occuring outcomes whose total probability is at least percent within that marginal.
        p = grid.p
        n = grid.n
        dict_of_lines = {}
        assert 0 < percent < 1
        for direction in point_of_plane.origin(p,n).gen_lines():
            lines_in_this_direction = top_lines_from_direction(grid, direction, percent)
            dict_of_lines[str(direction)] = lines_in_this_direction
        return functional_on_grid(dict_of_lines,p,n)

    @classmethod
    def functional_from_sl(cls,sl,grid,percent):
        #input is a grid, SL element, and percent. Returns a functional so that the lines in each direction sum to probability at least 'percent'
        #I've since found out that this method of creating functionals is doomed.
        p = grid.p
        n = grid.n
        dict_of_lines = {}
        rep_direction1 = None #first cannonical direction.
        rep_direction2 = None
        orbit = []
        for direction in point_of_plane.origin(p,n).gen_lines():
            assert np.allclose(grid.sum_line(direction), grid.sum_line(sl.inverse()*direction))
            if not direction in orbit:
                if len(orbit) == 0:
                    orbit = sl.orbit(direction)
                    lines_in_direction1 = top_lines_from_direction(grid, direction, percent, from_sl = True)
                    rep_direction1 = direction
                else:
                    lines_in_direction2 = top_lines_from_direction(grid, direction, percent, from_sl = True)
                    rep_direction2 = direction
                    break
        group_element = finite_sp_matrix.identity(2,p,n)
        while not group_element == -finite_sp_matrix.identity(2,p,n):
            dict_of_lines[str(group_element * rep_direction1)] = [group_element * l for l in lines_in_direction1]
            dict_of_lines[str(group_element * rep_direction2)] = [group_element * l for l in lines_in_direction2]
            group_element = sl * group_element
        return functional_on_grid(dict_of_lines,p,n)

    @classmethod
    def functional_counting_positive_points(cls, grid):
        p = grid.p
        n = grid.n
        dict_of_lines = {}
        for direction in point_of_plane.origin(p,n).gen_lines():
            lines_in_this_direction = []
            total_positive_count = 0
            for line in direction.gen_parallel_lines():
                positive_count = len([ 1 for pt in line.gen_points() if grid.get_value(pt) < 0])
                if positive_count > 0:
                    lines_in_this_direction.append((line, positive_count))
                    total_positive_count += positive_count
            dict_of_lines[str(direction)] = lines_in_this_direction
        f = functional_on_grid(dict_of_lines,p,n)
        f.point_count = total_positive_count
        return f

def sandbox_test():
    p,n = 7,2
    x = finite_field_element([0,1,0,0], p, 2*n)
    #y = finite_field_element([0,1], 3, 2)
    pow = (p**(2*n) -1) /(p**n - 1)
    y=x**pow
    print(pow)
    print(y)

def sl_decompose_test():
    #try to find a general pattern for the special sl matrices.
    p = 3
    n = 2
    I = finite_field_element.one(p,n)
    O = finite_field_element.zero(p,n)
    for A,B,C in itertools.product(finite_field_element.list_elements(p,n), finite_field_element.list_nonzero_elements(p,n), finite_field_element.list_elements(p,n)):
        mat = finite_sp_matrix([[I,O],[A,I]])*finite_sp_matrix([[B,O],[O,B.inverse()]])*finite_sp_matrix([[I,C],[O,I]])
        if mat.order()== p**n + 1:
            print(finite_sp_matrix([[I,O],[A,I]]).order())
            print(finite_sp_matrix([[B,O],[O,B.inverse()]]).order())
            print(finite_sp_matrix([[I,C],[O,I]]).order())
            #print(mat.order())
            print("---")
    for A,B in itertools.product(finite_field_element.list_elements(p,n), finite_field_element.list_nonzero_elements(p,n)):
        mat = finite_sp_matrix([[O,-B.inverse()],[B,A]])
        if mat.order()== p**n + 1:
            print(A)
            print(B)
            print("---")
        #functional = functional_on_grid.functional_from_sl(sl, grid, 0.1)
    #I = finite_matrix.identity(n,p,1)
#sl_decompose_test()
#sandbox_test()
#test_sl_from_extension()
#profile.run('sandbox_test()')
