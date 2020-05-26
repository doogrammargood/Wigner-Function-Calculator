from wigner_function import *
class grid_element(object):
    def __init__(self, matrix, p, n):
        assert len(matrix) == p**n
        self.p = p
        self.n = n
        self.values = [[ discrete_wig_fuct(point_of_plane((col, row)), matrix) for col in finite_field_element.list_elements(p,n)]for row in finite_field_element.list_elements(p,n)]

    def get_value(self, pt):
        return self.values[int(pt.y)][int(pt.x)]

    def marginalize_grid(self, line):
        lines = line_of_plane.gen_parallel_lines(line)
        marginal = []
        for l in lines:
            marginal.append(sum([self.get_value(pt) for pt in line.gen_points()]) )
        return marginal


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
test_grid()
