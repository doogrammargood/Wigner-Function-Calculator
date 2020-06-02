from wigner_function import *
class grid_element(object_modified):
    def __init__(self, matrix, p, n):
        assert len(matrix) == p**n
        self.p = p
        self.n = n
        self.values = [[ discrete_wig_fuct(point_of_plane((col, row)), matrix) for col in finite_field_element.list_elements(p,n)]for row in finite_field_element.list_elements(p,n)]

    def get_value(self, pt):
        if pt is None:
            return None
        return self.values[int(pt.y)][int(pt.x)]

    def marginalize_grid(self, line):
        return None
        if line is None:
            return None
        lines = line_of_plane.gen_parallel_lines(line)
        marginal = []
        for l in lines:
            marginal.append(sum([self.get_value(pt) for pt in line.gen_points()]) )
        return marginal

    def sum_line(self, line):
        if line is None:
            return None
        #returns the sum of the values over the line
        return sum([self.get_value(pt) for pt in line.gen_points()])

    def total_negativity(self):
        p,n =self.p, self.n
        return sum([ sum([ self.values[row][col] if self.values[row][col]<0 else 0 for col in range(p**n)]) for row in range(p**n)])
    def most_neg_pt(self):
        #returns a point where the value is minimized.
        val = float('inf')
        for x in finite_field_element.list_elements(self.p,self.n):
            for y in finite_field_element.list_elements(self.p,self.n):
                if val > self.get_value(point_of_plane((x,y))):
                    val = self.get_value(point_of_plane((x,y)))
                    current = point_of_plane((x,y))
        return current

def test_grid():
    p=5
    n=2
    matrix = random_pure_state(p,n)
    G=grid_element(matrix, p, n)
    a=5+3
    b=2
    A = finite_field_element([3,1],5,2)
    B = finite_field_element([2,0],5,2)
    pt = point_of_plane((A,B))
    assert G.get_value(pt) == G.values[b][a]
    assert G.get_value(pt) == discrete_wig_fuct(pt, matrix)
    #for l in line_of_plane.list_lines(p,n)
    zero = finite_field_element.zero(p,n)
    origin = point_of_plane((zero,zero))
    for l in origin.gen_lines():
        print(G.marginalize_grid(l))
#test_grid()
