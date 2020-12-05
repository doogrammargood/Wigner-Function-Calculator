# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'wig_calc_interface.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from main_window import *
import pickle
from finite_matrix_class import *
finite_matrix.load_dual_basis_matrices()
class Ui_MainWindow(object):

    def setupUi(self, MainWindow):
        MainWindow.setWindowIcon(QtGui.QIcon('./images/rrg_icons/SVG/Icons_Large_-03.svg'))
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1320, 689)
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


        self.verticalLayoutWidget_2.setObjectName("verticalLayoutWidget_2")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_2)
        self.verticalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.widget = WignerWidget(self.verticalLayoutWidget_2,p=MainWindow.p, n=MainWindow.n, grid = MainWindow.grid)

        self.widget.setObjectName("wignerWidget")
        self.verticalLayout_3.addWidget(self.widget)
        self.verticalLayoutWidget_2.raise_()
        self.gridLayout_2.addWidget(self.wigner_grid_frame, 0, 0, 1, 1)
        self.local_view_frame = QtWidgets.QFrame(self.tab)
        self.local_view_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.local_view_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.local_view_frame.setObjectName("local_view_frame")

        self.local_view_layout = QVBoxLayout(self.local_view_frame)
        self.local_view_widget = LocalView()
        self.local_view_layout.addWidget(self.local_view_widget)

        #self.gridLayout_2.addWidget(self.local_view_widget)
        self.local_view_widget.setObjectName("local_view_widget")
        self.gridLayout_2.addWidget(self.local_view_frame, 1, 0, 1, 1)
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
        self.entropy_label = QtWidgets.QLabel(self.frame)
        self.entropy_label.setObjectName("entropy_label")
        self.l1_label = QtWidgets.QLabel(self.frame)
        self.l1_label.setObjectName("l1_label")
        self.mana_label = QtWidgets.QLabel(self.frame)
        self.mana_label.setObjectName("mana_label")
        self.verticalLayout.addWidget(self.tot_neg_label)
        self.verticalLayout.addWidget(self.l1_label)
        self.verticalLayout.addWidget(self.mana_label)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        self.gridLayout_2.addWidget(self.frame, 2, 0, 1, 1)
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        #self.horizontal_layoutTab2 = QtWidgets.QHBoxLayout(self.tab_2)
        #self.verticalLayoutTab2 =QtWidgets.QVBoxLayout()
        #self.horizontal_layoutTab2.addLayout(self.verticalLayoutTab2)
        self.vertical_layout_tab_2 = QVBoxLayout()
        self.grid_layout_tab_2 = QGridLayout()
        self.tab_2.setLayout(self.vertical_layout_tab_2)
        self.vertical_layout_tab_2.addLayout(self.grid_layout_tab_2)

        #self.button_random.setGeometry(QtCore.QRect(40, 60, 211, 41))

        self.spinbox_layout = QHBoxLayout()
        self.p_spinbox = QtWidgets.QSpinBox(self.tab_2)
        self.p_spinbox.setObjectName("p_spinbox")
        self.n_spinbox = QtWidgets.QSpinBox(self.tab_2)
        self.n_spinbox.setObjectName("n_spinbox")
        self.spinbox_layout.addWidget(self.p_spinbox)
        self.spinbox_layout.addWidget(self.n_spinbox)
        #self.p_spinbox.setGeometry(QtCore.QRect(60, 30, 61, 27))
        self.button_change_size = QtWidgets.QPushButton(self.tab_2)
        self.button_change_size.setObjectName("button_change_size")
        self.button_zero = QtWidgets.QPushButton(self.tab_2)
        self.button_zero.setObjectName("button_zero")
        self.button_random = QtWidgets.QPushButton(self.tab_2)
        self.button_random.setObjectName("button_random")

        self.grid_layout_tab_2.addLayout(self.spinbox_layout,0,0)
        self.grid_layout_tab_2.addWidget(self.button_change_size,1,0)
        self.grid_layout_tab_2.addWidget(self.button_zero,2,0)
        self.grid_layout_tab_2.addWidget(self.button_random,3,0)

        self.button_maximally_mixed = QtWidgets.QPushButton(self.tab_2)
        self.button_maximally_mixed.setObjectName("button_maximally_mixed")
        self.button_superposition = QtWidgets.QPushButton(self.tab_2)
        self.button_superposition.setObjectName("button_superposition")
        self.button_superpositon_alternating = QtWidgets.QPushButton(self.tab_2)
        self.button_superpositon_alternating.setObjectName("button_superpositon_alternating")
        self.button_superpositon_zero_negative = QtWidgets.QPushButton(self.tab_2)
        self.button_superpositon_zero_negative.setObjectName("button_superpositon_zero_negative")

        self.grid_layout_tab_2.addWidget(self.button_maximally_mixed,0,1)
        self.grid_layout_tab_2.addWidget(self.button_superposition,1,1)
        self.grid_layout_tab_2.addWidget(self.button_superpositon_alternating,2,1)
        self.grid_layout_tab_2.addWidget(self.button_superpositon_zero_negative,3,1)


        self.special_state_layout = QHBoxLayout()
        self.button_special_state_left = QtWidgets.QPushButton(self.tab_2)
        self.button_special_state_right = QtWidgets.QPushButton(self.tab_2)
        self.spinbox_special_state = QtWidgets.QSpinBox(self.tab_2)

        self.special_state_layout.addWidget(self.button_special_state_left)
        self.special_state_layout.addWidget(self.button_special_state_right)
        self.special_state_layout.addWidget(self.spinbox_special_state)

        self.button_cat_state = QtWidgets.QPushButton(self.tab_2)
        self.button_cat_state.setObjectName("button_cat_state")

        self.grid_layout_tab_2.addLayout(self.special_state_layout, 0,2)
        self.grid_layout_tab_2.addWidget(self.button_cat_state,1,2)


        #self.n_spinbox.setGeometry(QtCore.QRect(160, 30, 61, 27))
        #self.button_change_size.setGeometry(QtCore.QRect(40, 160, 211, 41))
        #self.button_zero.setGeometry(QtCore.QRect(40, 110, 211, 41))
        # self.frame_2 = QtWidgets.QFrame(self.tab_2)
        # self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        # self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        # self.frame_2.setObjectName("frame_2")
        # self.frame_2_v_layout=QVBoxLayout()
        # self.frame_2.setLayout(self.frame_2_v_layout)
        #self.state_scroll_area = QScrollArea(self.frame_2)#Why doesn't this do anything?
        #self.frame_2_v_layout.addWidget(self.state_scroll_area)
        #self.label_state_info = QtWidgets.QLabel(self.frame_2)
        #self.frame_2_v_layout.addWidget(self.label_state_info)
        #self.label_state_info.setObjectName("label_state_info")
        #self.label_state_info.setFixedWidth(1000)
        #self.vertical_layout_tab_2.addWidget(self.frame_2)
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
        self.button_change_size.setText(_translate("MainWindow", "change size"))
        self.button_zero.setText(_translate("MainWindow", "zero state"))
        self.button_random.setText(_translate("MainWindow", "random pure state"))
        self.button_maximally_mixed.setText(_translate("MainWindow", "maximally mixed state"))
        self.button_superposition.setText(_translate("MainWindow", "equiphase superposition state"))
        self.button_superpositon_alternating.setText(_translate("MainWindow", "superposition with alternating phase"))
        self.button_superpositon_zero_negative.setText(_translate("MainWindow", "superposition with zero negative"))
        self.button_cat_state.setText(_translate("MainWindow", "cat state"))
        self.button_special_state_left.setText(_translate("MainWindow", "<-"))
        self.button_special_state_right.setText(_translate("MainWindow", "->"))
        #self.label_state_info.setText(_translate("MainWindow", "Current State"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MainWindow", "Density Matrix"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuEdit.setTitle(_translate("MainWindow", "Edit"))
        self.actionSave.setText(_translate("MainWindow", "Save"))

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon('./images/rrg_icons/SVG/Icons_Large_-03.svg'))
    #MainWindow = QtWidgets.QMainWindow()
    MainWindow = MyMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow._init2()
    MainWindow.connect_signals()
    MainWindow.set_labels()
    MainWindow.show()
    sys.exit(app.exec_())
