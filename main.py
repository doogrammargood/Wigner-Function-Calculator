#From the helpful moonsweeper tutorial: https://www.learnpyqt.com/examples/moonsweeper/

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
#from sl_class import *
from wigner_function import *

IMG_FLAG = QImage("./images/flag.png")

def find_color(value, min, max):
    #assume value between min and max.
    assert value <= max and value >= min
    interval = max - min
    moved_val = value - min
    ratio = moved_val / interval
    return QColor(255*(1-ratio), 0, 255*ratio)

NUM_COLORS = {
    0: Qt.lightGray,
    1: QColor('#f44336'),
    2: QColor('#9C27B0'),
    3: QColor('#3F51B5'),
    4: QColor('#03A9F4'),
    5: QColor('#00BCD4'),
    6: QColor('#4CAF50'),
    7: QColor('#E91E63'),
    8: QColor('#FF9800'),
    9: QColor('#FBF9800'),
}



class Pos(QWidget):
    clicked = pyqtSignal(point_of_plane)

    def __init__(self, pt, *args, **kwargs):
        super(Pos, self).__init__(*args, **kwargs)

        self.setFixedSize(QSize(20, 20))
        self.is_flagged = False# flagged points are the ones that have been clicked
        self.is_marked = False# Marked to show when points are on a line.
        self.value = 0 #determines the shading.

        self.pt = pt
        self.x = pt.x
        self.y = pt.y
        self.dim = self.x.p**self.x.n

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)

        r = event.rect()

        inner = find_color(self.value,-1/self.dim, 1/self.dim)
        if self.is_marked:
            outer = Qt.black
        else:
            outer = Qt.gray
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


    def flag(self):
        if self.is_flagged:
            self.is_flagged = False
        else:
            self.is_flagged = True
        self.update()

    def set_mark(self):
        self.is_marked = True
        self.update()
    def set_unmark(self):
        self.is_marked = False
        self.update()

    def mouseReleaseEvent(self, e):
        self.flag()
        self.clicked.emit(self.pt)

    def wig_function(self, mat):
        self.value = discrete_wig_fuct(self.pt,mat)
        if self.pt.x.is_zero():
            print(self.value)
        self.update()


class MainWindow(QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.flagged_pos=[]
        self.marked_pos=[]
        #self.b_size=25
        w = QWidget()
        hb = QHBoxLayout()


        vb = QVBoxLayout()
        vb.addLayout(hb)

        self.grid = QGridLayout()
        self.grid.setSpacing(5)

        vb.addLayout(self.grid)
        w.setLayout(vb)
        self.setCentralWidget(w)

        self.p=5
        self.n=2
        self.init_map(self.p,self.n)
        self.show()

    def init_map(self,p,n):
        # Add positions to the map
        for x in finite_field_element.list_elements(p,n):
            for y in finite_field_element.list_elements(p,n):
                w = Pos(point_of_plane((x,y)))
                self.grid.addWidget(w, int(x), int(y))
                # Connect signal to handle expansion.
                w.clicked.connect(self.flags_changed)


    def flags_changed(self, *args):
        pt = args[0]
        if not pt in self.flagged_pos:
            self.flagged_pos.append(pt)
        else:
            self.flagged_pos.remove(pt)
        self.update_markings()

    def update_markings(self):
        new_markings = []
        for p1,p2 in itertools.combinations(self.flagged_pos,2):
            line = p1.line_to(p2)
            for pt_to_mark in line.gen_points():
                new_markings.append(pt_to_mark)
        old_markings = self.marked_pos
        for pt in [o for o in old_markings if not o in new_markings]:
            w = self.grid.itemAtPosition(int(pt.x), int(pt.y)).widget()
            w.set_unmark()
        for pt in [n for n in new_markings if not n in old_markings]:
            w = self.grid.itemAtPosition(int(pt.x), int(pt.y)).widget()
            w.set_mark()
        self.marked_pos = new_markings

    def set_flags_sl(self,pt):
        generator = sl_matrix.gen_with_order(self.p,self.n)
        sl = next(generator)
        positions_to_mark = sl.orbit(pt)
        for p in positions_to_mark:
            w = self.grid.itemAtPosition(int(p.x), int(p.y)).widget()
            w.flag()
            self.flagged_pos.append(p)
        self.update_markings()

    def wig_function(self, matrix):
        p,n = self.p, self.n
        for x in range(p**n):
            for y in range(p**n):
                w = self.grid.itemAtPosition(x,y).widget()
                w.wig_function(matrix)

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
if __name__ == '__main__':
    app = QApplication([])
    window = MainWindow()
    one = finite_field_element.one(window.p,window.n)
    pt = point_of_plane((one,one))
    #window.set_flags_sl(pt)
    matrix = random_pure_state(window.p, window.n)
    window.wig_function(matrix)
    print(window.values())
    app.exec_()
