#From the helpful moonsweeper tutorial: https://www.learnpyqt.com/examples/moonsweeper/

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
#from sl_class import *
from wigner_function import *
from grid_class import *
import math
IMG_FLAG = QImage("./images/flag.png")

def find_color(value, min, max):
    epsilon = 0.001
    assert value <= max + epsilon and value >= min - epsilon
    interval = max - min
    moved_val = value - min
    ratio = moved_val / interval
    ratio = ratio * math.pi / 2

    return QColor(255*(math.cos(ratio)), 255*math.cos(ratio), 255*math.sin(ratio))


class Pos(QWidget):
    clicked = pyqtSignal(point_of_plane)
    hovered = pyqtSignal(point_of_plane)
    unhovered = pyqtSignal()

    def __init__(self, pt, *args, **kwargs):
        super(QWidget, self).__init__(*args, **kwargs)

        self.is_flagged = False# flagged points are the ones that have been clicked
        self.is_marked = False# Marked to show when points are on a line.
        self.value = 0 #determines the shading.
        if pt is None:
            self.dim = 5**2
            self.magnify = True
            self.setFixedSize(QSize(60, 60))
            self.pt = None
            self.x = None
            self.y = None
        else:
            self.pt = pt
            self.x = pt.x
            self.y = pt.y
            self.dim = self.x.p**self.x.n
            self.magnify = False
            if self.dim < 15:
                self.setFixedSize(QSize(30, 30))
            elif self.dim < 30:
                self.setFixedSize(QSize(20, 20))
            elif self.dim < 60:
                self.setFixedSize(QSize(10, 10))
            else:
                self.setFixedSize(QSize(5, 5))
            self.installEventFilter(self)

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)

        inner = find_color(self.value,-1/self.dim, 1/self.dim)
        if self.is_marked:
            outer = Qt.black
        else:
            outer = Qt.gray
        r = event.rect()
        p.fillRect(r, QBrush(inner))
        pen = QPen(outer)
        if self.is_marked:
            pen.setWidth(7)
        else:
            pen.setWidth(1)
        p.setPen(pen)
        p.drawRect(r)
        if self.is_flagged:
            p.drawPixmap(r, QPixmap(IMG_FLAG))

    def copy_data(self, other):
        self.is_flagged = other.is_flagged
        self.is_marked = other.is_marked
        self.value = other.value
        self.pt = other.pt
        self.x = other.x
        self.y = other.y
        self.dim = other.dim
        self.update()

    def flag(self):
        self.is_flagged = True
        self.update()

    def unflag(self):
        self.is_flagged = False
        self.update()

    def set_mark(self):
        self.is_marked = True
        self.update()
    def set_unmark(self):
        self.is_marked = False
        self.update()
    def eventFilter(self, object, event):
        if event.type() == QEvent.MouseButtonPress:
            self.clicked.emit(self.pt)
        elif event.type() == QEvent.Enter:
            self.hovered.emit(self.pt)
        elif event.type() == QEvent.Leave:
            self.unhovered.emit()
#print ("left ")
        return False


class WignerWidget(QWidget):
    def __init__(self, *args, **kwargs):
        super(QWidget, self).__init__(*args)
        self.flagged_pts=[]
        self.marked_pts=[]
        #self.b_size=25
        #w = QWidget()
        hb = QHBoxLayout()
        vb = QVBoxLayout()
        vb.addLayout(hb)

        self.grid = QGridLayout()
        self.grid.setSpacing(2)

        vb.addLayout(self.grid)
        self.setLayout(vb)
        #w.setLayout(vb)
        #self.setCentralWidget(w)
        self.p=kwargs['p']
        self.n=kwargs['n']
        if self.p**self.n > 60:
            self.grid.setSpacing(1)
        self.init_map(self.p,self.n)

        grid = kwargs['grid']
        self.set_values_from_grid(grid)
        #self.show()

    def init_map(self,p,n):
        # Add positions to the map
        for x in finite_field_element.list_elements(p,n):
            for y in finite_field_element.list_elements(p,n):
                w = Pos(point_of_plane((x,y)))
                self.grid.addWidget(w, int(x), int(y))
                # Connect signal to handle expansion.
                #w.clicked.connect(self.flags_changed)


    def set_flagged(self, flags):
        flagged =[f for f in flags if not f.isNone()]
        for pt in self.flagged_pts:
            w = self.grid.itemAtPosition(int(pt.x), int(pt.y)).widget()
            w.unflag()

        for pt in flagged:
            w = self.grid.itemAtPosition(int(pt.x), int(pt.y)).widget()
            w.flag()
        self.flagged_pts = flagged

    def set_markings(self, new_marks):
        new_markings = [m for m in new_marks if not m.isNone()]
        for pt in self.marked_pts:
            w = self.grid.itemAtPosition(int(pt.x), int(pt.y)).widget()
            w.set_unmark()
        for pt in new_markings:
            w = self.grid.itemAtPosition(int(pt.x), int(pt.y)).widget()
            w.set_mark()
        self.marked_pts = new_markings

    def set_flags_sl(self,pt):
        generator = sl_matrix.gen_with_order(self.p,self.n)
        sl = next(generator)
        positions_to_mark = sl.orbit(pt)
        for p in positions_to_mark:
            w = self.grid.itemAtPosition(int(p.x), int(p.y)).widget()
            w.flag()
            self.flagged_pos.append(p)
        self.update_markings()

    def set_values_from_grid(self, grid):
        p,n = self.p, self.n
        for x in finite_field_element.list_elements(p,n):
            for y in finite_field_element.list_elements(p,n):
                w = self.grid.itemAtPosition(int(x),int(y)).widget()
                w.value = grid.get_value(point_of_plane((x,y)))
                w.update()

    def values(self):
        #prints a list of lists for the values of points.
        to_return = []
        p,n = self.p, self.n
        for x in range(p**n):
            row = []
            for y in range(p**n):
                w = self.grid.itemAtPosition(x,y).widget()
                row.append(w.value)
            to_return.append(row)
        return to_return

class LocalView(QWidget):
    def __init__(self, *args, **kwargs):
        super(QWidget, self).__init__(*args)
        self.pt1 = point_of_plane(None)
        self.value1 = None
        self.pt2 = point_of_plane(None)
        self.value2 = None
        self.line = line_of_plane(None)
        self.valuel = None
        self.marginal = None
        self.pos1 = Pos(None)

        self.pos2 = Pos(None)
        vb = QVBoxLayout()
        self.setLayout(vb)
        #self.label_pt.setAlignment(Qt.AlignCenter)
        #self.label_pt.setWordWrap(True)
        #self.setGeometry(50,50,30,10)
        self.label_pt1 = QLabel()
        self.label_value1 = QLabel()
        self.label_pt2 = QLabel()
        self.label_value2 = QLabel()
        self.label_line = QLabel()
        self.label_valuel = QLabel()
        self.label_marginal = QLabel()
        h1 = QHBoxLayout()
        v1 = QVBoxLayout()
        vb.addLayout(h1)
        h1.addWidget(self.pos1)
        v1.addWidget(self.label_pt1)
        v1.addWidget(self.label_value1)
        h1.addLayout(v1)

        h2 = QHBoxLayout()
        v2 = QVBoxLayout()
        vb.addLayout(h2)
        h2.addWidget(self.pos2)

        v2.addWidget(self.label_pt2)
        v2.addWidget(self.label_value2)
        h2.addLayout(v2)

        vb.addWidget(self.label_line)
        vb.addWidget(self.label_valuel)
        vb.addWidget(self.label_marginal)

        self.set_labels()
        self.show()
    def set_values(self, pt1, value1, pos1, pt2, value2, pos2, line, valuel, marginal):
        self.pt1 = pt1
        self.value1 = value1
        self.pos1.copy_data(pos1)
        self.pt2 = pt2
        self.value2 = value2
        self.pos2.copy_data(pos2)
        self.line = line
        self.valuel = valuel
        self.marginal = marginal
        self.set_labels()
        self.update()

    def set_labels(self):
        self.label_pt1.setText("point 1 = (" + str(self.pt1) + ")")
        self.label_value1.setText("Value of pt1 ="+ str(self.value1) )
        self.label_pt2.setText("point 2 = (" + str(self.pt2) + ")")
        self.label_value2.setText("Value of pt2 = "+ str(self.value2) )
        self.label_line.setText("line = " + str(self.line))
        self.label_valuel.setText("Value of Line = " + str(self.valuel))
        self.label_marginal = self.draw_marginal()

    def draw_marginal(self):
        #TODO: find a way to represent the maginal distribution, like a bar graph.
        return QLabel()

    # def self.paint_event(self):
    #     p = QPainter(self)
    #     p.setRenderHint(QPainter.Antialiasing)
    #     r = event.rect()
    #
    #     inner = find_color(self.value,-1/self.dim, 1/self.dim)
    #     if self.is_marked:
    #         outer = Qt.black
    #     else:
    #         outer = Qt.gray
    #     p.fillRect(r, QBrush(inner))
    #     pen = QPen(outer)
    #     pen.setWidth(5)
    #
    #     p.setPen(pen)
    #     p.drawRect(r)
    #     if self.is_flagged:
    #         p.drawPixmap(r, QPixmap(IMG_FLAG))
