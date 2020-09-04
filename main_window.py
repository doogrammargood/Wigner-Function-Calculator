from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
#from PyQt5.QtCore import *
from main import *
from density_matrix_functions import *
class MyMainWindow(QMainWindow):
    def __init__(self, *args, **kwargs):
        super(QMainWindow, self).__init__(*args, **kwargs)
        self.p = 19
        self.n = 1
        #self.density_matrix = super_position_state_negatives(self.p, self.n)
        self.density_matrix = random_pure_state(self.p,self.n)
        self.grid = grid_element(self.density_matrix, self.p, self.n) #this is potentially confusing: grid is not a layout.
        self.entropy_value = 1
        self.clear_data()
    def clear_data(self):
        self.pt1 = point_of_plane(None)
        self.pt2 = point_of_plane(None)
        self.pos1 = Pos(None)
        self.pos2 = Pos(None)
        self.line = line_of_plane(None)

    def _init2(self):
        #to be called after the ui has been set up.
        self.wig = self.findChild(QWidget, "wignerWidget")
        self.local_view = self.findChild(QWidget,"local_view_widget")
        self.p_spinbox = self.findChild(QWidget,"p_spinbox")
        self.n_spinbox = self.findChild(QWidget,"n_spinbox")
        self.p_spinbox.setMinimum(3)
        self.p_spinbox.setMaximum(97)
        self.p_spinbox_value = self.p
        self.p_spinbox.setValue(self.p)
        self.n_spinbox.setMinimum(1)
        self.n_spinbox.setMaximum(5)
        self.n_spinbox.setValue(self.n)

    def connect_signals(self, reconnect_wig = False):
        #establishes the connection with
        p,n = self.wig.p, self.wig.n
        for x in finite_field_element.list_elements(p,n):
            for y in finite_field_element.list_elements(p,n):
                w = self.wig.grid.itemAtPosition(int(x),int(y)).widget()
                w.hovered.connect(self.handle_hover)
                w.clicked.connect(self.handle_click)
                w.right_clicked.connect(self.handle_right_click)

        if not reconnect_wig:
            random_button = self.findChild(QPushButton,"button_random")
            random_button.clicked.connect(self.handle_random_button)
            zero_button = self.findChild(QPushButton,"button_zero")
            zero_button.clicked.connect(self.handle_zero_button)
            p_spin = self.findChild(QSpinBox,"p_spinbox")
            p_spin.valueChanged.connect(self.handle_p_spinbox_change)
            change_size_button = self.findChild(QPushButton, "button_change_size")
            change_size_button.clicked.connect(self.handle_change_size)
            maximally_mixed_button = self.findChild(QPushButton, "button_maximally_mixed")
            maximally_mixed_button.clicked.connect(self.handle_maximally_mixed_button)
            superposition_button = self.findChild(QPushButton, "button_superposition")
            superposition_button.clicked.connect(self.handle_superpositon_button)
            superposition_alternating_button = self.findChild(QPushButton, "button_superpositon_alternating")
            superposition_alternating_button.clicked.connect(self.handle_superposition_alternating_button)
            superposition_zero_negative_button = self.findChild(QPushButton, "button_superpositon_zero_negative")
            superposition_zero_negative_button.clicked.connect(self.handle_superposition_zero_negative_button)
            cat_state_button = self.findChild(QPushButton, "button_cat_state")
            cat_state_button.clicked.connect(self.handle_cat_state_button)

    def set_labels(self):
        tot_neg = self.findChild(QLabel, "tot_neg_label")
        most_neg_pt = self.findChild(QLabel, "most_neg_pt")
        state_info = self.findChild(QLabel, "label_state_info")
        tot_neg.setText("Total Negativity = " + str(self.grid.total_negativity()))
        most_neg_pt.setText("Most Negative Point = " +str(self.grid.most_neg_pt()))
        state_info.setText(str(self.density_matrix))

    def change_matrix(self, new_matrix):
        self.density_matrix = new_matrix
        self.grid = grid_element(self.density_matrix, self.p, self.n)
        self.wig.set_values_from_grid(self.grid)
        self.set_labels()

    def change_size(self, p, n):
        self.clear_data()
        self.p=p
        self.n=n
        self.density_matrix = random_pure_state(self.p,self.n)
        self.grid = grid_element(self.density_matrix, self.p, self.n)
        self.flagged_pos = []

        current = self.wig
        layout = current.parent().layout()
        for i in reversed(range(layout.count())):
            layout.itemAt(i).widget().setParent(None)
        w = WignerWidget(current.parent(), p=self.p, n=self.n, grid =self.grid)
        w.setObjectName("wignerWidget")
        self.wig = w
        layout.addWidget(w)
        #w.show()

    def keyPressEvent(self, e):
        if e.key() == QtCore.Qt.Key_Escape:
            self.close()
        elif e.key() == QtCore.Qt.Key_Space:
            #self.change_matrix(state_with_special_order(self.p))
            self.change_matrix(labyrinth_state(self.p,self.n))
            #self.change_matrix(None)
        elif e.key() == QtCore.Qt.Key_A:
            #mat = cat_state(self.p,self.n)
            mat = super_position_state_negatives(self.p, self.n)
            self.change_matrix(mat)
        elif e.key() == QtCore.Qt.Key_S:
            self.save("test_file.wig")
        elif e.key() == QtCore.Qt.Key_L:
            self.load("test_file.wig")
        elif e.key() == QtCore.Qt.Key_U:
            tab = self.findChild(QTabWidget,"tabWidget")
            tab.setCurrentIndex(1)
        elif e.key() == QtCore.Qt.Key_O:
            tab = self.findChild(QTabWidget,"tabWidget")
            tab.setCurrentIndex(0)
        elif e.key()== QtCore.Qt.Key_P:
            matrix = weil_elementary(-1,0,self.p) @ self.density_matrix @ weil_elementary(1,0,self.p)
            self.change_matrix(matrix)
        elif e.key()== QtCore.Qt.Key_Y:
            matrix = weil_elementary(-1,-1,self.p) @ self.density_matrix @ weil_elementary(1,1,self.p)
            self.change_matrix(matrix)

    def handle_click(self,pt):
        #This function is ugly. N
        pos = self.wig.grid.itemAtPosition(int(pt.x),int(pt.y)).widget()
        if self.pt1.isNone():
            self.pt1 = pt
            self.wig.set_decorators('flagged',[self.pt1, self.pt2])
            self.pos1.copy_data(pos)
        elif pt == self.pt1:
            self.pt1 = self.pt2
            self.pt2 = point_of_plane(None)
            self.wig.set_decorators('flagged',[self.pt1, self.pt2])
            self.pos1.copy_data(self.pos2)
            self.pos2.copy_data(Pos(None))
            self.pos1.set_decorator('marked',val=False)
            self.line = line_of_plane(None)
        elif pt == self.pt2:
            self.pt2 = point_of_plane(None)
            self.wig.set_decorators('flagged',[self.pt1, self.pt2])
            self.pos2.copy_data(Pos(None))
            self.pos1.set_decorator('marked',val=False)
            self.line = line_of_plane(None)
        else:
            self.pt2 = pt
            self.wig.set_decorators('flagged',[self.pt1, self.pt2])
            self.pos2.copy_data(pos)
            self.pos1.set_decorator('marked')
            self.line = self.pt1.line_to(self.pt2)
        self.update_views(pt,pos)

    def handle_right_click(self,pt):
        pos = self.wig.grid.itemAtPosition(int(pt.x),int(pt.y)).widget()
        decorated = self.wig.decorators['highlighted']
        if pt not in decorated:
            self.wig.set_decorators('highlighted',decorated+[pt])
        else:
            self.wig.set_decorators('highlighted',[p for p in decorated if not p == pt])
        #self.update_views(pt,pos)

    def handle_hover(self,pt):
        pos = self.wig.grid.itemAtPosition(int(pt.x),int(pt.y)).widget()
        self.update_views(pt,pos)

    def handle_change_size(self):
        self.change_size(self.p_spinbox.value(),self.n_spinbox.value())
        self.connect_signals(reconnect_wig = True)
    def handle_random_button(self):
        r = random_pure_state(self.p, self.n)
        self.change_matrix(r)
    def handle_zero_button(self):
        z = zero_state(self.p,self.n)
        self.change_matrix(z)
    def handle_maximally_mixed_button(self):
        m = maximally_mixed_state(self.p,self.n)
        self.change_matrix(m)
    def handle_superpositon_button(self):
        s = super_position_state(self.p,self.n)
        self.change_matrix(s)
    def handle_superposition_alternating_button(self):
        s = super_position_state_negatives(self.p,self.n)
        self.change_matrix(s)
    def handle_superposition_zero_negative_button(self):
        s = super_position_zero_negative(self.p,self.n)
        self.change_matrix(s)
    def handle_cat_state_button(self):
        c = cat_state(self.p,self.n)
        self.change_matrix(c)

    def handle_p_spinbox_change(self,new_val):
        if new_val > self.p_spinbox_value:
            new_val = next_prime(new_val)
            self.p_spinbox.setValue(new_val)
            self.p_spinbox_value = new_val
        else:
            new_val = previous_prime(new_val)
            self.p_spinbox.setValue(new_val)
            self.p_spinbox_value = new_val

    def update_views(self,pt,pos):
        #determines what the local view controller should see.

        #self.pt1 = pt
        if self.pt1.isNone():
            self.wig.set_decorators('marked',[])
            value1 = self.grid.get_value(pt)
            self.local_view.set_values(pt, value1, pos, point_of_plane(None), None, Pos(None), line_of_plane(None), None, None, None)
        elif self.pt2.isNone():
            value1 = self.grid.get_value(self.pt1)
            value2 = self.grid.get_value(pt)
            line = self.pt1.line_to(pt)
            valuel = self.grid.sum_line(line)
            marginal = self.grid.marginalize_grid(line)
            if self.pt1.is_on_line(line):
                self.pos1.set_decorator('marked')
            else:
                self.pos1.set_decorator('marked', val = False)
            pos.set_decorator('marked')
            entropy = self.grid.entropy(2, line)
            self.local_view.set_values(self.pt1, value1, self.pos1, pt, value2, pos, line, valuel, marginal, (self.entropy_value, entropy) )
            self.wig.set_decorators('marked',list(line.gen_points()))
        else:
            value1 = self.grid.get_value(self.pt1)
            value2 = self.grid.get_value(self.pt2)
            line = self.line.parallel_through(pt)
            valuel = self.grid.sum_line(line)
            marginal = self.grid.marginalize_grid(line)
            if self.pt1.is_on_line(line):
                self.pos1.set_decorator('mark')
            else:
                self.pos1.set_decorator('mark',val = False)
            if self.pt2.is_on_line(line):
                self.pos2.set_decorator('mark')
            else:
                self.pos2.set_decorator('mark',val=False)
            entropy = self.grid.entropy(self.entropy_value, line)
            self.local_view.set_values(self.pt1, value1, self.pos1, self.pt2, value2, self.pos2, line, valuel, marginal, (self.entropy_value, entropy))
            self.wig.set_decorators('marked',list(line.gen_points()))

    def save(self,filename):
        file = open(filename, "wb")
        pickle.dump(self.dictionary_state(), file)

    def load(self,filename):
        file = open(filename, "rb")
        dict = pickle.load(file)
        self.update_from_dict(dict)
        self.update()

    def dictionary_state(self):
        #produces a dictionary which represents the state of the program.
        dict = {'p': self.p, 'n': self.n, 'density_matrix': self.density_matrix,
                'pt1': str(self.pt1), 'pt2': str(self.pt2), 'grid': self.grid.dictionary_state()}
        return dict

    def update_from_dict(self,dict):
        self.p = int(dict['p'])
        self.n = int(dict['n'])
        self.change_size(self.p, self.n)
        self.density_matrix = dict['density_matrix']
        self.change_matrix(self.density_matrix)
        self.grid = grid_element(self.density_matrix,self.p,self.n)
        self.pt1 = point_of_plane.from_string( dict['pt1'], self.p, self.n)
        self.pt2 = point_of_plane.from_string( dict['pt2'], self.p, self.n)
        pos1 = self.get_pos(self.pt1)
        pos1.flag()
        self.pos1.flag()
        self.pos1.copy_data(pos1)
        pos2 = self.get_pos(self.pt2)
        pos2.flag()
        self.pos2.copy_data(pos2)
        self.pos2.flag()


    def get_pos(self,pt):
        if pt.isNone():
            return Pos(None)
        else:
            return self.wig.grid.itemAtPosition(int(x),int(y)).widget()
