#From the helpful moonsweeper tutorial: https://www.learnpyqt.com/examples/moonsweeper/

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from affine_plane_class import *

IMG_FLAG = QImage("./images/flag.png")

NUM_COLORS = {
    0: Qt.lightGray,
    1: QColor('#f44336'),
    2: QColor('#9C27B0'),
    3: QColor('#3F51B5'),
    4: QColor('#03A9F4'),
    5: QColor('#00BCD4'),
    6: QColor('#4CAF50'),
    7: QColor('#E91E63'),
    8: QColor('#FF9800')
}



class Pos(QWidget):
    clicked = pyqtSignal(point_of_plane)

    def __init__(self, pt, *args, **kwargs):
        super(Pos, self).__init__(*args, **kwargs)

        self.setFixedSize(QSize(20, 20))
        self.is_flagged = False# flagged points are the ones that have been clicked
        self.level = 0 #determines the shading.
        self.pt = pt
        self.x = pt.x
        self.y = pt.y

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)

        r = event.rect()

        inner = NUM_COLORS[self.level]
        outer = Qt.gray
        p.fillRect(r, QBrush(inner))
        pen = QPen(outer)
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

    def set_mark(self, value):
        self.level = value
        self.update()

    def mouseReleaseEvent(self, e):
        self.flag()
        self.clicked.emit(self.pt)
        #self.clicked.emit()


class MainWindow(QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.flagged_pos=[]
        self.marked_pos=[]#this is going to consist of pairs ((point), int)
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
        # for pts in self.flagged_pos:
        #     print(pts.x,pts.y)
        self.update_markings()

    def update_markings(self):
        new_markings = []
        for p1,p2 in itertools.combinations(self.flagged_pos,2):
            line = p1.line_to(p2)
            for pt_to_mark in line.gen_points():
                new_markings.append(pt_to_mark)
        new_markings = [(pt, new_markings.count(pt)) for pt in new_markings]
        old_markings = self.marked_pos
        for pt, val in old_markings:
            w = self.grid.itemAtPosition(int(pt.x), int(pt.y)).widget()
            w.set_mark(0)
        for pt, val in new_markings:
            w = self.grid.itemAtPosition(int(pt.x), int(pt.y)).widget()
            w.set_mark(val)
        self.marked_pos = new_markings



if __name__ == '__main__':
    app = QApplication([])
    window = MainWindow()
    app.exec_()
