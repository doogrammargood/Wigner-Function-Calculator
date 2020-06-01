2# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'wig_calc_interface.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from main import *
import pickle
class Ui_MainWindow(object):

    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1020, 689)
        MainWindow.setStatusTip("")
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setMovable(True)
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.tab)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.wigner_grid_frame = QtWidgets.QFrame(self.tab)
        self.wigner_grid_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.wigner_grid_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.wigner_grid_frame.setObjectName("wigner_grid_frame")
        self.verticalLayoutWidget_2 = QtWidgets.QWidget(self.wigner_grid_frame)
        #self.verticalLayoutWidget_2.setGeometry(QtCore.QRect(-1, -1, 491, 581))

        self.verticalLayoutWidget_2.setObjectName("verticalLayoutWidget_2")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_2)
        self.verticalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        #self.widget = QtWidgets.QWidget(self.verticalLayoutWidget_2)
        self.widget = WignerWidget(self.verticalLayoutWidget_2,p=MainWindow.p, n=MainWindow.n, grid = MainWindow.grid)

        self.widget.setObjectName("wignerWidget")
        self.verticalLayout_3.addWidget(self.widget)
        self.verticalLayoutWidget_2.raise_()
        self.gridLayout_2.addWidget(self.wigner_grid_frame, 0, 0, 2, 1)
        self.local_view_frame = QtWidgets.QFrame(self.tab)
        self.local_view_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.local_view_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.local_view_frame.setObjectName("local_view_frame")

        self.local_view_layout = QVBoxLayout(self.local_view_frame)
        self.local_view_widget = LocalView()
        self.local_view_layout.addWidget(self.local_view_widget)

        #self.gridLayout_2.addWidget(self.local_view_widget)
        self.local_view_widget.setObjectName("local_view_widget")
        self.gridLayout_2.addWidget(self.local_view_frame, 0, 1, 1, 1)
        self.frame = QtWidgets.QFrame(self.tab)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.frame)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setContentsMargins(-1, 15, -1, 15)
        self.verticalLayout.setSpacing(0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.most_neg_pt = QtWidgets.QLabel(self.frame)
        self.most_neg_pt.setObjectName("most_neg_pt")
        self.verticalLayout.addWidget(self.most_neg_pt)
        self.tot_neg_label = QtWidgets.QLabel(self.frame)
        self.tot_neg_label.setObjectName("tot_neg_label")
        self.verticalLayout.addWidget(self.tot_neg_label)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        self.gridLayout_2.addWidget(self.frame, 1, 1, 1, 1)
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.button_random = QtWidgets.QPushButton(self.tab_2)
        self.button_random.setGeometry(QtCore.QRect(40, 60, 211, 41))
        self.button_random.setObjectName("button_random")
        self.p_spinbox = QtWidgets.QSpinBox(self.tab_2)
        self.p_spinbox.setGeometry(QtCore.QRect(60, 30, 61, 27))
        self.p_spinbox.setObjectName("p_spinbox")
        self.n_spinbox = QtWidgets.QSpinBox(self.tab_2)
        self.n_spinbox.setGeometry(QtCore.QRect(160, 30, 61, 27))
        self.n_spinbox.setObjectName("n_spinbox")
        self.button_zero = QtWidgets.QPushButton(self.tab_2)
        self.button_zero.setGeometry(QtCore.QRect(40, 110, 211, 41))
        self.button_zero.setObjectName("button_zero")
        self.frame_2 = QtWidgets.QFrame(self.tab_2)
        self.frame_2.setGeometry(QtCore.QRect(49, 169, 201, 221))
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.label = QtWidgets.QLabel(self.frame_2)
        self.label.setGeometry(QtCore.QRect(10, 10, 171, 201))
        self.label.setObjectName("label")
        self.tabWidget.addTab(self.tab_2, "")
        self.gridLayout.addWidget(self.tabWidget, 0, 0, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1020, 25))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuEdit = QtWidgets.QMenu(self.menubar)
        self.menuEdit.setObjectName("menuEdit")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionSave = QtWidgets.QAction(MainWindow)
        self.actionSave.setObjectName("actionSave")
        self.menuFile.addAction(self.actionSave)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuEdit.menuAction())

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Wigner Function Calculator"))
        self.most_neg_pt.setText(_translate("MainWindow", "Most Negative Point: "))
        self.tot_neg_label.setText(_translate("MainWindow", "Total Negativity:"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MainWindow", "Wigner Function"))
        self.button_random.setText(_translate("MainWindow", "random pure state"))
        self.button_zero.setText(_translate("MainWindow", "zero state"))
        self.label.setText(_translate("MainWindow", "Current State"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MainWindow", "Density Matrix"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuEdit.setTitle(_translate("MainWindow", "Edit"))
        self.actionSave.setText(_translate("MainWindow", "Save"))

class MyMainWindow(QMainWindow):
    def __init__(self, *args, **kwargs):
        super(QMainWindow, self).__init__(*args, **kwargs)
        self.p = 3
        self.n = 2
        self.density_matrix = random_pure_state(self.p,self.n)
        self.grid = grid_element(self.density_matrix, self.p, self.n) #this is potentially confusing: grid is not a layout.
        self.pt1 = point_of_plane(None)
        self.pt2 = point_of_plane(None)
        self.pos1 = Pos(None)
        self.pos2 = Pos(None)
        self.line = line_of_plane(None)

    def dictionary_state(self):
        #produces a dictionary which represents the state of the program.
        dict = {'p': self.p, 'n': self.n, 'density_matrix': self.density_matrix,
                'pt1': str(self.pt1), 'pt2': str(self.pt2) }
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
        wig = self.findChild(QWidget, "wignerWidget")
        pos1 = self.get_pos(self.pt1)
        pos1.flag()
        self.pos1.flag()
        self.pos1.copy_data(pos1)
        pos2 = self.get_pos(self.pt2)
        pos2.flag()
        self.pos2.copy_data(pos2)
        self.pos2.flag()

    def change_matrix(self, new_matrix):
        self.density_matrix = new_matrix
        self.grid = grid_element(self.density_matrix, self.p, self.n)
        w = self.findChild(QWidget, "wignerWidget")
        w.set_values_from_grid(self.grid)
        self.set_labels()

    def change_size(self, p, n):
        self.p=p
        self.n=n
        self.density_matrix = random_pure_state(self.p,self.n)
        self.grid = grid_element(self.density_matrix, self.p, self.n)
        current = self.findChild(QWidget, "wignerWidget")
        layout = current.parent().layout()
        for i in reversed(range(layout.count())):
            layout.itemAt(i).widget().setParent(None)
        w = WignerWidget(p=self.p, n=self.n, grid =self.grid)
        w.setObjectName("wignerWidget")
        layout.addWidget(w)
        w.show()

    def keyPressEvent(self, e):
        if e.key() == QtCore.Qt.Key_Escape:
            self.close()
        elif e.key() == QtCore.Qt.Key_Space:
            self.change_matrix(None)
        elif e.key() == QtCore.Qt.Key_A:
            self.change_size(3, 2)
        elif e.key() == QtCore.Qt.Key_S:
            self.save("test_file.wig")
        elif e.key() == QtCore.Qt.Key_L:
            self.load("test_file.wig")

    def handle_click(self,pt):
        #This function is ugly. N
        wig = self.findChild(QWidget, "wignerWidget")
        pos = wig.grid.itemAtPosition(int(pt.x),int(pt.y)).widget()
        if self.pt1.isNone():
            self.pt1 = pt
            wig.set_flagged([self.pt1, self.pt2])
            self.pos1.copy_data(pos)
        elif pt == self.pt1:
            self.pt1 = self.pt2
            self.pt2 = point_of_plane(None)
            wig.set_flagged([self.pt1, self.pt2])
            self.pos1.copy_data(self.pos2)
            self.pos2.copy_data(Pos(None))
            self.pos1.set_unmark()
            self.line = line_of_plane(None)
        elif pt == self.pt2:
            self.pt2 = point_of_plane(None)
            wig.set_flagged([self.pt1, self.pt2])
            self.pos2.copy_data(Pos(None))
            self.pos1.set_unmark()
            self.line = line_of_plane(None)
        else:
            self.pt2 = pt
            wig.set_flagged([self.pt1, self.pt2])
            self.pos2.copy_data(pos)
            self.pos1.set_mark()
            self.line = self.pt1.line_to(self.pt2)
        self.update_views(pt,pos)


    def handle_hover(self,pt):
        wig = self.findChild(QWidget, "wignerWidget")
        pos = wig.grid.itemAtPosition(int(pt.x),int(pt.y)).widget()
        self.update_views(pt,pos)
        # if self.pt2.isNone():
        #     self.update_views(pt,pos)

    def set_labels(self):
        tot_neg = self.findChild(QLabel, "tot_neg_label")
        most_neg_pt = self.findChild(QLabel, "most_neg_pt")
        tot_neg.setText("Total Negativity = " + str(self.grid.total_negativity()))
        most_neg_pt.setText("Most Negative Point = " +str(self.grid.most_neg_pt()))

    def update_views(self,pt,pos):
        #determines what the local view controller should see.

        local_view = self.findChild(QWidget, "local_view_widget")
        wig = self.findChild(QWidget, "wignerWidget")
        #self.pt1 = pt
        if self.pt1.isNone():
            wig.set_markings([])
            value1 = self.grid.get_value(pt)
            local_view.set_values(pt, value1, self.pos1, point_of_plane(None), None, Pos(None), line_of_plane(None), None, None)
        elif self.pt2.isNone():
            value1 = self.grid.get_value(self.pt1)
            value2 = self.grid.get_value(pt)
            line = self.pt1.line_to(pt)
            valuel = self.grid.sum_line(line)
            marginal = self.grid.marginalize_grid(line)
            local_view.set_values(self.pt1, value1, self.pos1, pt, value2, pos, line, valuel, marginal)
            wig.set_markings(list(line.gen_points()))
        else:
            value1 = self.grid.get_value(self.pt1)
            value2 = self.grid.get_value(self.pt2)
            line = self.line
            valuel = self.grid.sum_line(line)
            marginal = self.grid.marginalize_grid(line)
            local_view.set_values(self.pt1, value1, self.pos1, self.pt2, value2, self.pos2, self.line, valuel, marginal)
            wig.set_markings(list(self.line.gen_points()))

    def connect_signals(self):
        #establishes the connection with
        local_view = self.findChild(QWidget, "local_view_widget")
        wignerWidget = self.findChild(QWidget, "wignerWidget")
        p,n = wignerWidget.p, wignerWidget.n
        for x in finite_field_element.list_elements(p,n):
            for y in finite_field_element.list_elements(p,n):
                w = wignerWidget.grid.itemAtPosition(int(x),int(y)).widget()
                w.hovered.connect(self.handle_hover)
                w.clicked.connect(self.handle_click)

    def save(self,filename):
        file = open(filename, "wb")
        pickle.dump(self.dictionary_state(), file)

    def load(self,filename):
        file = open(filename, "rb")
        dict = pickle.load(file)
        self.update_from_dict(dict)
        self.update()

    def get_pos(self,pt):
        wignerWidget = self.findChild(QWidget, "wignerWidget")
        if pt.isNone():
            return Pos(None)
        else:
            return wignerWidget.grid.itemAtPosition(int(x),int(y)).widget()

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    #MainWindow = QtWidgets.QMainWindow()
    MainWindow = MyMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.connect_signals()
    MainWindow.set_labels()
    MainWindow.show()
    sys.exit(app.exec_())