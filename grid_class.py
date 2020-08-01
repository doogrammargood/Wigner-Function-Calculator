from wigner_function import *
import multiprocessing
from density_matrix_functions import *
import profile
#finite_matrix.load_dual_basis_matrices()
class grid_element(object_modified):
    def __init__(self, matrix, p, n, compute_fresh = True, multiparticle = False):
        pure = np.allclose(matrix, matrix @ matrix)
        if compute_fresh:
            assert len(matrix) == p**n
            self.p = p
            self.n = n
            if pure:
                vect = pure_state_from_density_matrix(matrix)
                self.values = [[ discrete_wig_fuct_pure_state(point_of_plane((col, row)), vect, multiparticle = multiparticle) for col in finite_field_element.list_elements(p,n)]for row in finite_field_element.list_elements(p,n)]
            else:
                self.values = [[ discrete_wig_fuct(point_of_plane((col, row)), matrix, multiparticle = multiparticle) for col in finite_field_element.list_elements(p,n)]for row in finite_field_element.list_elements(p,n)]
            # with multiprocessing.Pool(processes = 4) as poo:
            #     self.values= poo.starmap(discrete_wig_fuct, [(point_of_plane((col,row)),matrix) for col, row in itertools.product(finite_field_element.list_elements(p,n), repeat = 2)] )
            # self.values = [self.values[i*(p**n):(i+1)*(p**n)] for i in range(p**n)]
            # self.values = [list(x) for x in zip(*self.values)]

            self.marginals = {}
            self.record_marginals()
            #assert np.allclose(mat_in, matrix)
        else:
            pass
    def get_value(self, pt):
        if pt is None:
            return None
        return self.values[int(pt.y)][int(pt.x)]

    def record_marginals(self):
        for line in point_of_plane.origin(self.p, self.n).gen_lines():
            self.marginals[str(line)] = self._marginalize_grid(line)
    def _marginalize_grid(self, line):
        epsilon = 10**-7
        if line is None:
            return None
        marginal = []
        for l in line.gen_parallel_lines():
            marginal.append(sum([self.get_value(pt) for pt in l.gen_points()]) )
        #assert(abs (sum(marginal) -1) < epsilon)
        return marginal

    def marginalize_grid(self, line):
        if line.isNone():
            return []
        cached_line = line.parallel_through(point_of_plane.origin(self.p,self.n))
        return self.marginals[str(cached_line)]

    def sum_line(self, line):
        if line.isNone():
            return 0
        cached_line = line.parallel_through(point_of_plane.origin(self.p,self.n))
        if not line.coefficients[0].is_zero():
            x = int(line.coefficients[2]/line.coefficients[0])
        else:
            x = int(line.coefficients[2]/line.coefficients[1])
        #assert self.marginals[str(cached_line)][x] == sum([self.get_value(pt) for pt in line.gen_points()])
        return self.marginals[str(cached_line)][x]

    def total_negativity(self):
        p,n =self.p, self.n
        return np.sum([ np.sum([ self.values[row][col] if self.values[row][col]<0 else 0 for col in range(p**n)]) for row in range(p**n)])
    def most_neg_pt(self):
        #returns a point where the value is minimized.
        val = float('inf')
        for x in finite_field_element.list_elements(self.p,self.n):
            for y in finite_field_element.list_elements(self.p,self.n):
                if val > self.get_value(point_of_plane((x,y))):
                    val = self.get_value(point_of_plane((x,y)))
                    current = point_of_plane((x,y))
        return current

    def dictionary_state(self):
        #returns a dictionary with its information. For saving.
        dict = {}
        dict['p']=self.p
        dict['n']=self.n
        dict['values']=self.values
        dict['marginals']=self.marginals
        return dict

    @classmethod
    def from_dictonary_state(cls, dict):
        to_return = grid_element(compute_fresh=False)
        to_return.p=dict['p']
        to_return.n=dict['n']
        to_return.value = dict['values']
        to_rreturn.marginals = dict['marginals']
        return to_return

def test_grid():
    p=5
    n=1
    matrix = super_position_state_negatives(p,n)
    G=grid_element(matrix, p, n)
    origin = point_of_plane.origin(p,n)
    for direction in origin.gen_lines():
        for line in direction.gen_parallel_lines():
            tot = np.zeros((p**n,p**n))
            for pt in line.gen_points():
                tot = tot + phase_pt_general(pt)
            print(line)
            print(np.linalg.matrix_rank(tot, tol = 10**-3))
            assert is_hermitian(tot)
            #assert np.linalg.matrix_rank(tot) == 1
def test_wig_fuct_pure_state():
    p=5
    n=2
    matrix = super_position_state_negatives(p,n)
    #G2=grid_element(matrix, p, n, pure = True)
    G1=grid_element(matrix, p, n)
    # for x in finite_field_element.list_elements(p,n):
    #     for y in finite_field_element.list_elements(p,n):
    #         pt = point_of_plane((x,y))
    #         assert np.allclose(G1.get_value(pt), G2.get_value(pt))
#profile.run('test_wig_fuct_pure_state()')
