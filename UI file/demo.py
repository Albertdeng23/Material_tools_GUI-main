# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'demo.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_AseAtomInput(object):
    def setupUi(self, AseAtomInput):
        AseAtomInput.setObjectName("AseAtomInput")
        AseAtomInput.resize(686, 424)
        self.centralwidget = QtWidgets.QWidget(AseAtomInput)
        self.centralwidget.setObjectName("centralwidget")
        self.graphicsView = QtWidgets.QGraphicsView(self.centralwidget)
        self.graphicsView.setGeometry(QtCore.QRect(380, 40, 261, 241))
        self.graphicsView.setObjectName("graphicsView")
        self.layoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.layoutWidget.setGeometry(QtCore.QRect(10, 330, 651, 30))
        self.layoutWidget.setObjectName("layoutWidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.layoutWidget)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.AddAtom_Buttom = QtWidgets.QPushButton(self.layoutWidget)
        self.AddAtom_Buttom.setObjectName("AddAtom_Buttom")
        self.horizontalLayout.addWidget(self.AddAtom_Buttom)
        self.DelAtom = QtWidgets.QPushButton(self.layoutWidget)
        self.DelAtom.setObjectName("DelAtom")
        self.horizontalLayout.addWidget(self.DelAtom)
        self.ViewStructure_Button = QtWidgets.QPushButton(self.layoutWidget)
        self.ViewStructure_Button.setObjectName("ViewStructure_Button")
        self.horizontalLayout.addWidget(self.ViewStructure_Button)
        self.SetCaluator_Button = QtWidgets.QPushButton(self.layoutWidget)
        self.SetCaluator_Button.setObjectName("SetCaluator_Button")
        self.horizontalLayout.addWidget(self.SetCaluator_Button)
        self.Process_Button = QtWidgets.QPushButton(self.layoutWidget)
        self.Process_Button.setObjectName("Process_Button")
        self.horizontalLayout.addWidget(self.Process_Button)
        self.Cancel_Buttom = QtWidgets.QPushButton(self.layoutWidget)
        self.Cancel_Buttom.setObjectName("Cancel_Buttom")
        self.horizontalLayout.addWidget(self.Cancel_Buttom)
        self.layoutWidget1 = QtWidgets.QWidget(self.centralwidget)
        self.layoutWidget1.setGeometry(QtCore.QRect(10, 10, 351, 287))
        self.layoutWidget1.setObjectName("layoutWidget1")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.layoutWidget1)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtWidgets.QLabel(self.layoutWidget1)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.AtomInPut = QtWidgets.QLineEdit(self.layoutWidget1)
        self.AtomInPut.setObjectName("AtomInPut")
        self.verticalLayout.addWidget(self.AtomInPut)
        self.label_2 = QtWidgets.QLabel(self.layoutWidget1)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.X_Input = QtWidgets.QLineEdit(self.layoutWidget1)
        self.X_Input.setText("")
        self.X_Input.setObjectName("X_Input")
        self.verticalLayout.addWidget(self.X_Input)
        self.label_3 = QtWidgets.QLabel(self.layoutWidget1)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.verticalLayout.addWidget(self.label_3)
        self.Y_Input = QtWidgets.QLineEdit(self.layoutWidget1)
        self.Y_Input.setText("")
        self.Y_Input.setObjectName("Y_Input")
        self.verticalLayout.addWidget(self.Y_Input)
        self.label_4 = QtWidgets.QLabel(self.layoutWidget1)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_4.setFont(font)
        self.label_4.setObjectName("label_4")
        self.verticalLayout.addWidget(self.label_4)
        self.Z_Input = QtWidgets.QLineEdit(self.layoutWidget1)
        self.Z_Input.setText("")
        self.Z_Input.setObjectName("Z_Input")
        self.verticalLayout.addWidget(self.Z_Input)
        self.label_5 = QtWidgets.QLabel(self.layoutWidget1)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.verticalLayout.addWidget(self.label_5)
        self.Caluator_Input = QtWidgets.QLineEdit(self.layoutWidget1)
        self.Caluator_Input.setText("")
        self.Caluator_Input.setObjectName("Caluator_Input")
        self.verticalLayout.addWidget(self.Caluator_Input)
        AseAtomInput.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(AseAtomInput)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 686, 26))
        self.menubar.setObjectName("menubar")
        self.menu = QtWidgets.QMenu(self.menubar)
        self.menu.setObjectName("menu")
        self.menuConnect_Server = QtWidgets.QMenu(self.menubar)
        self.menuConnect_Server.setObjectName("menuConnect_Server")
        self.menuAbout = QtWidgets.QMenu(self.menubar)
        self.menuAbout.setObjectName("menuAbout")
        self.menuPlot = QtWidgets.QMenu(self.menubar)
        self.menuPlot.setObjectName("menuPlot")
        AseAtomInput.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(AseAtomInput)
        self.statusbar.setObjectName("statusbar")
        AseAtomInput.setStatusBar(self.statusbar)
        self.actionFrom = QtWidgets.QAction(AseAtomInput)
        self.actionFrom.setObjectName("actionFrom")
        self.actionServer = QtWidgets.QAction(AseAtomInput)
        self.actionServer.setObjectName("actionServer")
        self.actionVision = QtWidgets.QAction(AseAtomInput)
        self.actionVision.setObjectName("actionVision")
        self.actionAbout_Author = QtWidgets.QAction(AseAtomInput)
        self.actionAbout_Author.setObjectName("actionAbout_Author")
        self.actionBand_Structure = QtWidgets.QAction(AseAtomInput)
        self.actionBand_Structure.setObjectName("actionBand_Structure")
        self.menu.addAction(self.actionFrom)
        self.menuConnect_Server.addAction(self.actionServer)
        self.menuAbout.addAction(self.actionVision)
        self.menuAbout.addAction(self.actionAbout_Author)
        self.menuPlot.addAction(self.actionBand_Structure)
        self.menubar.addAction(self.menu.menuAction())
        self.menubar.addAction(self.menuConnect_Server.menuAction())
        self.menubar.addAction(self.menuAbout.menuAction())
        self.menubar.addAction(self.menuPlot.menuAction())

        self.retranslateUi(AseAtomInput)
        QtCore.QMetaObject.connectSlotsByName(AseAtomInput)

    def retranslateUi(self, AseAtomInput):
        _translate = QtCore.QCoreApplication.translate
        AseAtomInput.setWindowTitle(_translate("AseAtomInput", "MainWindow"))
        self.AddAtom_Buttom.setText(_translate("AseAtomInput", "Add Atom"))
        self.DelAtom.setText(_translate("AseAtomInput", "Del Atom"))
        self.ViewStructure_Button.setText(_translate("AseAtomInput", "View Structure"))
        self.SetCaluator_Button.setText(_translate("AseAtomInput", "Set Calculator"))
        self.Process_Button.setText(_translate("AseAtomInput", "Process"))
        self.Cancel_Buttom.setText(_translate("AseAtomInput", "Cancel"))
        self.label.setText(_translate("AseAtomInput", "Atom"))
        self.label_2.setText(_translate("AseAtomInput", "X"))
        self.label_3.setText(_translate("AseAtomInput", "Y"))
        self.label_4.setText(_translate("AseAtomInput", "Z"))
        self.label_5.setText(_translate("AseAtomInput", "Calculator"))
        self.menu.setTitle(_translate("AseAtomInput", "File"))
        self.menuConnect_Server.setTitle(_translate("AseAtomInput", "VASP"))
        self.menuAbout.setTitle(_translate("AseAtomInput", "About"))
        self.menuPlot.setTitle(_translate("AseAtomInput", "Plot"))
        self.actionFrom.setText(_translate("AseAtomInput", "Import cif"))
        self.actionServer.setText(_translate("AseAtomInput", "Connect Server"))
        self.actionVision.setText(_translate("AseAtomInput", "Version"))
        self.actionAbout_Author.setText(_translate("AseAtomInput", "About Author"))
        self.actionBand_Structure.setText(_translate("AseAtomInput", "Band Structure"))