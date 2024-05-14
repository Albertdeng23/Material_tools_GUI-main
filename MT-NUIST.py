import paramiko
import sys
import os
import time
import numpy as np
from matplotlib.backends.backend_qt5agg import  FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from functools import partial
from PyQt5.QtWidgets import *
from ase.visualize import view
from ase import Atoms
from ase.calculators.emt import EMT
from ase.calculators.vasp import Vasp
from ase.io import read
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from CIFDownload import CIFDownloader
from DescriptorDesign import FileSelectionWidget


### TODO     
#   1.使用VASP生成POSCAR , INCAR 之类的输入文件

###----------------------------------------------------------------------------------------------------------####
##工具


def timer(func):
    def wrapper(*args, **kwargs):
        app = QApplication([])
        window = QWidget()
        layout = QVBoxLayout()
        
        label = QLabel()
        layout.addWidget(label)
        window.setLayout(layout)
        window.setWindowTitle("执行时间统计")
        window.setGeometry(100, 100, 200, 100)
        window.show()
        
        def update_label(elapsed_time):
            label.setText(f"执行时间: {elapsed_time:.4f} 秒")
        
        def execute_func():
            start_time = time.time()
            result = func(*args, **kwargs)
            elapsed_time = time.time() - start_time
            update_label(elapsed_time)
            return result
        
        # 创建一个子线程来执行函数，以不阻塞主线程的方式实时更新窗口
        from threading import Thread
        thread = Thread(target=execute_func)
        thread.start()
        
        app.exec_()
        
    return wrapper


### 绘制能带图
def plot_band_structure(eigenval_file, figure):
    # 读取EIGENVAL文件
    with open(eigenval_file, 'r') as f:
        data = f.readlines()

    # 获取能带和能级数据
    efermi = float(data[5].split()[0])  # 费米能级
    eigenval_raw = np.loadtxt(eigenval_file, skiprows=8, usecols=(1,)) - efermi

    # 绘制能带图
    nkpoints = len(eigenval_raw)
    kpoints = np.linspace(0, 1, nkpoints)

    ax = figure.add_subplot(111)
    ax.plot(kpoints, eigenval_raw, color='b')

    # 设置坐标轴和标签
    ax.axhline(y=0, color='k', ls='--')
    ax.set_xlabel('k-points')
    ax.set_ylabel('Energy (eV)')

    # 显示能级线
    for i in range(1, nkpoints - 1):
        if kpoints[i] % 0.25 == 0:
            ax.axvline(x=kpoints[i], color='k', ls='--')

    # 显示费米能级
    ax.axhline(y=0, color='r', linewidth=0.5, ls='-')

    # 添加标题
    ax.set_title('Band Structure')
###----------------------------------------------------------------------------------------------------------####

### 绘制能带图类
class Band_Structure(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("能带图绘制")
        self.setGeometry(100, 100, 400, 300)

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.create_main_layout()

    def create_main_layout(self):
    # Apply main style to the main widget
        self.setStyleSheet("""
            QWidget {
                font-size: 14px;
            }
            QPushButton {
                background-color: #5AA469;
                color: white;
                border-style: outset;
                border-width: 2px;
                border-radius: 10px;
                border-color: beige;
                font: bold 14px;
                padding: 6px;
                margin: 6px;
            }
            QPushButton:hover {
                background-color: #66c280;
            }
            QLabel {
                min-height: 3em;
                align: center;
                font: 14px;
            }
            QScrollArea {
                border: 0;
            }
        """)

        main_widget = QWidget(self)
        self.setCentralWidget(main_widget)
        main_layout = QVBoxLayout(main_widget)
        main_layout.setContentsMargins(20, 20, 20, 20)
        main_layout.setSpacing(10)

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        main_layout.addWidget(scroll_area)

        scroll_widget = QWidget()
        scroll_area.setWidget(scroll_widget)

        # Add the QVBoxLayout to the scrollable widget
        scroll_layout = QVBoxLayout(scroll_widget)
        self.image_label = QLabel()
        self.image_label.setAlignment(Qt.AlignCenter)
        self.image_label.setMinimumSize(400, 300)  # Minimum size for image display
        scroll_layout.addWidget(self.image_label)

        # Layout for buttons
        button_layout = QHBoxLayout()
        scroll_layout.addLayout(button_layout)

        self.file_button = QPushButton('选择文件')
        self.file_button.clicked.connect(partial(self.select_file, "EIGENVAL"))
        button_layout.addWidget(self.file_button)

        self.plot_button = QPushButton('开始绘制')
        self.plot_button.clicked.connect(self.plot_band_structure)
        button_layout.addWidget(self.plot_button)

        self.save_button = QPushButton('保存图片')
        self.save_button.clicked.connect(self.save_figure)
        button_layout.addWidget(self.save_button)

        # Label for showing messages or information
        self.label = QLabel()
        self.label.setAlignment(Qt.AlignCenter)
        scroll_layout.addWidget(self.label)

    def select_file(self, file_type):
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(self, '选择文件')
        if file_path:
            setattr(self, file_type, file_path)
            self.label.setText(f"{file_type}文件已选择：{file_path}")

    def plot_band_structure(self):
        if not hasattr(self, 'EIGENVAL'):
            self.label.setText("请先选择EIGENVAL文件！")
            return
            
        eigenval_file = self.EIGENVAL
        plot_band_structure(eigenval_file, self.figure)
        self.update_image_label()

    def update_image_label(self):
        self.canvas.draw()
        w, h = self.canvas.get_width_height()
        image = QImage(self.canvas.buffer_rgba(), w, h, QImage.Format_RGBA8888)
        self.image_label.setPixmap(QPixmap.fromImage(image))

    def save_figure(self):
        save_dialog = QFileDialog()
        save_path, _ = save_dialog.getSaveFileName(self, '保存图片', '', 'PNG (*.png);;JPEG (*.jpg *.jpeg)')
        if save_path:
            self.figure.savefig(save_path)
            self.label.setText(f"图片已保存：{save_path}")


### mainUI类
class Ui_AseAtomInput(object):
    def setupUi(self, AseAtomInput):
        AseAtomInput.setObjectName("AseAtomInput")
        AseAtomInput.resize(686, 418)
        styleSheet = """
            QWidget {
                font-size: 14px;
            }
            QPushButton {
                background-color: #5AA469;
                color: white;
                border-style: outset;
                border-width: 2px;
                border-radius: 10px;
                border-color: beige;
                font: bold 14px;
                padding: 6px;
                margin: 1px;
            }
            QPushButton:hover {
                background-color: #66c280;
            }
            QLineEdit {
                font: 14px Arial;
                border: 1px solid #ccc;
                padding: 3px;
                margin: 4px;
            }
            QLabel {
                font: bold italic;
                color: #555;
            }
            QGraphicsView {
                border-style: outset;
                border-width: 2px;
                border-radius: 5px;
                border-color: #ccc;
            }
            QMenuBar {
                background-color: #5AA469;
                color: white;
            }
            QMenuBar::item {
                background-color: #5AA469;
                color: white;
            }
            QMenuBar::item::selected {
                background-color: #66c280;
            }
            QMenu {
                background-color: white;
                color: #5AA469;
            }
            QMenu::item::selected {
                background-color: #66c280;
                color: white;
            }
        """
        AseAtomInput.setStyleSheet(styleSheet)
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
        self.menubar.setGeometry(QtCore.QRect(0, 0, 686, 22))
        self.menubar.setObjectName("menubar")
        self.menu = QtWidgets.QMenu(self.menubar)
        self.menu.setObjectName("menu")
        self.menuConnect_Server = QtWidgets.QMenu(self.menubar)
        self.menuConnect_Server.setObjectName("menuConnect_Server")
        self.menuAbout = QtWidgets.QMenu(self.menubar)
        self.menuAbout.setObjectName("menuAbout")
        self.menuPlot = QtWidgets.QMenu(self.menubar)
        self.menuPlot.setObjectName("menuPlot")
        self.menuMaterials_Project = QtWidgets.QMenu(self.menubar)
        self.menuMaterials_Project.setObjectName("menuMaterials_Project")
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
        self.actionCIF_Download = QtWidgets.QAction(AseAtomInput)
        self.actionCIF_Download.setObjectName("actionCIF_Download")
        self.actionDescriptor_Design = QtWidgets.QAction(AseAtomInput)
        self.actionDescriptor_Design.setObjectName("actionDescriptor_Design")
        self.menu.addAction(self.actionFrom)
        self.menuConnect_Server.addAction(self.actionServer)
        self.menuAbout.addAction(self.actionVision)
        self.menuAbout.addAction(self.actionAbout_Author)
        self.menuPlot.addAction(self.actionBand_Structure)
        self.menuMaterials_Project.addAction(self.actionCIF_Download)
        self.menuMaterials_Project.addAction(self.actionDescriptor_Design)
        self.menubar.addAction(self.menu.menuAction())
        self.menubar.addAction(self.menuConnect_Server.menuAction())
        self.menubar.addAction(self.menuAbout.menuAction())
        self.menubar.addAction(self.menuPlot.menuAction())
        self.menubar.addAction(self.menuMaterials_Project.menuAction())

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
        self.menuMaterials_Project.setTitle(_translate("AseAtomInput", "Materials Project"))
        self.actionFrom.setText(_translate("AseAtomInput", "Import cif"))
        self.actionServer.setText(_translate("AseAtomInput", "Connect Server"))
        self.actionVision.setText(_translate("AseAtomInput", "Version"))
        self.actionAbout_Author.setText(_translate("AseAtomInput", "About Author"))
        self.actionBand_Structure.setText(_translate("AseAtomInput", "Band Structure"))
        self.actionCIF_Download.setText(_translate("AseAtomInput", "CIF Download"))
        self.actionDescriptor_Design.setText(_translate("AseAtomInput", "Descriptor Design"))

### verison类
class MyApp(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        # 设置窗口基本属性
        self.setWindowTitle('Software Information')
        self.setGeometry(300, 300, 500, 350)
        self.setMinimumSize(500, 350)

        # 创建并应用背景
        self.setStyleSheet("background-color: #F5DEB3;")

        # 创建布局
        layout = QVBoxLayout()

        # 创建并设置版本标签
        version_label = QLabel('Version 1.2')
        version_label.setFont(QFont('Arial', 20, QFont.Bold))
        version_label.setStyleSheet("color: #FFFFFF; margin-bottom: 20px;")
        version_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(version_label)

        # 创建并设置GitHub标签
        github_label = QLabel('<a href="https://github.com/Albertdeng23/Material_tools_GUI-main">GitHub: https://github.com/Albertdeng23/Material_tools_GUI-main</a>')
        github_label.setFont(QFont('Arial', 16))
        github_label.setStyleSheet("color: #FFFFFF;")

        # 设置标签的背景颜色为透明，覆盖默认的背景颜色
        github_label.setAutoFillBackground(True)
        p = github_label.palette()
        p.setColor(github_label.backgroundRole(), Qt.transparent)
        github_label.setPalette(p)

        github_label.setOpenExternalLinks(True)  # 设置标签可以打开外部链接
        github_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(github_label)

        # 创建并设置Release标签
        release_label = QLabel('<a href="https://github.com/Albertdeng23/Material_tools_GUI-main/releases">Release: https://github.com/Albertdeng23/Material_tools_GUI-main/releases</a>')
        release_label.setFont(QFont('Arial', 16))
        release_label.setStyleSheet("color: #FFFFFF; margin-top: 10px;")

        # 设置标签的背景颜色为透明，覆盖默认的背景颜色
        release_label.setAutoFillBackground(True)
        p = release_label.palette()
        p.setColor(release_label.backgroundRole(), Qt.transparent)
        release_label.setPalette(p)

        release_label.setOpenExternalLinks(True)  # 设置标签可以打开外部链接
        release_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(release_label)

        # 应用布局
        self.setLayout(layout)

## about author类
class AboutAuthorDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('About Author')
        self.setModal(True)

        self.author_label = QLabel('Author: Albertdeng23')
        self.author_label.setFont(QFont('Arial', 14, QFont.Bold))
        self.author_label.setAlignment(Qt.AlignCenter)

        self.email_label = QLabel('Email: issicdeng@outlook.com')
        self.email_label.setFont(QFont('Arial', 12))
        self.email_label.setAlignment(Qt.AlignCenter)

        self.github_label = QLabel('<a href="https://github.com/Albertdeng23">GitHub: https://github.com/Albertdeng23</a>')
        self.github_label.setFont(QFont('Arial', 12))
        self.github_label.setAlignment(Qt.AlignCenter)
        self.github_label.setOpenExternalLinks(True)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok)
        button_box.accepted.connect(self.accept)

        layout = QVBoxLayout()
        layout.addWidget(self.author_label)
        layout.addWidget(self.email_label)
        layout.addWidget(self.github_label)
        layout.addWidget(button_box)
        self.setLayout(layout)

### 文件传输进度条
class ProgressUpdater(QObject):
    progress_signal = pyqtSignal(int)

    def __init__(self, total_files):
        super().__init__()
        self.total_files = total_files

    def update_progress(self):
        for index in range(1, self.total_files + 1):
            self.progress_signal.emit(index)
            time.sleep(0.01)  # 模拟传输延迟
class ThreadedProgressUpdater(QRunnable):
    def __init__(self, progress_updater):
        super().__init__()
        self.progress_updater = progress_updater

    def run(self):
        self.progress_updater.update_progress()

### 创建MainWindow
class ASE_ui(Ui_AseAtomInput,QMainWindow): 
    def __init__(self,) -> None:
        super().__init__()
        self.atom_calculator = AtomCalculator()
        self.myAppInstance = MyApp()
        self.about_author_dialog = AboutAuthorDialog(self)
        self.initUI()

    def initUI(self):
        self.setupUi(self)
        self.retranslateUi(self)

        self.DelAtom.clicked.connect(self.deleteAtom)
        self.AddAtom_Buttom.clicked.connect(self.addAtom)
        self.ViewStructure_Button.clicked.connect(self.showStructure)
        self.SetCaluator_Button.clicked.connect(self.setCalculator)
        self.actionFrom.triggered.connect(self.importCIF)
        self.actionServer.triggered.connect(self.connectRemoteServer)
        self.Process_Button.clicked.connect(self.startCalculation)
        self.actionVision.triggered.connect(self.versionButton)
        self.actionAbout_Author.triggered.connect(self.showAboutAuthorDialog)
        self.actionBand_Structure.triggered.connect(self.plot_Band)
        self.actionCIF_Download.triggered.connect(self.showCFIDownloader)
        self.actionDescriptor_Design.triggered.connect(self.showDescriptorDesign)
    def showAboutAuthorDialog(self):
        self.about_author_dialog.exec_()

    def addAtom(self):
        element = self.AtomInPut.text()
        x = float(self.X_Input.text())
        y = float(self.Y_Input.text())
        z = float(self.Z_Input.text())
        self.atom_calculator.addAtom(element, x, y, z)
        self.AtomInPut.clear()
        self.X_Input.clear()
        self.X_Input.clear()
        self.X_Input.clear()

    def versionButton(self):
        self.myAppInstance.show()

    ### 对应：“删除上一个原子”按钮
    def deleteAtom(self):
        self.atom_calculator.deleteAtom()

    ### 对应：“结束输入”
    def finishInput(self):
        self.close()

    ### 对应：“显示原子结构”
    def showStructure(self):
        if self.atom_calculator.atoms is not None:
            view(self.atom_calculator.atoms)

    ### 对应：“设置计算器”按钮
    def setCalculator(self):
        calculator_name = self.Caluator_Input.text()
        self.atom_calculator.setCalculator(calculator_name)
        if calculator_name.lower() == "vasp":
            self.showVaspParamDialog()
        else:
            show_Calculator = QErrorMessage(self)
            show_Calculator.showMessage(f"caulater: {calculator_name}")
    
    ### 对应：Band_structure 
    def plot_Band(self):
        self.Band_structure = Band_Structure()
        self.Band_structure.show()

    ### 点击“设置计算器”按钮之后如果是“VASP”计算器就显示输入计算参数
    def showVaspParamDialog(self):
        self.dialog = VaspParamDialog(self)
        if self.dialog.exec_() == QDialog.Accepted:
            self.vasp_params = self.dialog.getParams()
            show_Calculator = QErrorMessage(self)
            show_Calculator.showMessage(f"Caulater: vasp, paragram: {self.vasp_params}")


    def setOptimizer(self):
        optimizer_name = self.optimizer_input.text()
        self.atom_calculator.setOptimizer(optimizer_name)
        show_Optimizer = QErrorMessage(self)
        show_Optimizer.showMessage(f"优化器为: {optimizer_name}")

    ### 对应：“开始计算”按钮
    def startCalculation(self):
        self.atom_calculator.startCalculation()
        energy = self.atom_calculator.energy
        forces = self.atom_calculator.forces
        message = QErrorMessage(self)
        message.showMessage(f"energy: {energy}, force为: {forces}")

    def importCIF(self):
        dialog = CIFImportDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            atoms = dialog.getAtoms()
            if atoms is not None:
                self.atom_calculator.atoms = atoms
                show_Structure = QErrorMessage(self)
                show_Structure.showMessage('Successfully Improt Atom Info')
    
    def connectRemoteServer(self):
        dialog = RemoteServerDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            server_info = dialog.getServerInfo()
            if server_info is not None:
                self.hostname = server_info['hostname']
                self.username = server_info['username']
                self.password = server_info['password']
                self.file_path = server_info['folder']


                try:
                    self.ssh = paramiko.SSHClient()
                    self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
                    self.ssh.connect(self.hostname, username=self.username, password=self.password)

                    use_times = self.atom_calculator.run_vasp_calculation(self.ssh , self.file_path)
                    QMessageBox.information(self, f'successful', f'successfully connect to {self.hostname} \n Result file path : .\\Results \n Use Time: {use_times /60 : .2f} minute. ')



                except paramiko.AuthenticationException:
                    QMessageBox.information(self, f'Authentication Error','please check you account and password')
                    self.connectRemoteServer()

                except Exception as e:
                    error_dialog = QErrorMessage(self)
                    error_dialog.showMessage(f'connect faild：{str(e)}')

    def showCFIDownloader(self):
        self.CIFDownload = CIFDownloader()
        self.CIFDownload.show()

    def showDescriptorDesign(self):
        self.DescriptorDesigner = FileSelectionWidget()
        self.DescriptorDesigner.show()

### DONE 链接远程服务器类，用于获取远程服务器的连接信息
class RemoteServerDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Connect Remote Server')
        self.setModal(True)
        self.server_info = None

        self.hostname_label = QLabel('Hostname:')
        self.username_label = QLabel('Account:')
        self.password_label = QLabel('Password:')

        self.hostname_input = QLineEdit()
        self.username_input = QLineEdit()
        self.password_input = QLineEdit()
        self.password_input.setEchoMode(QLineEdit.Password)

        self.folder_button = QPushButton('Select Folder')
        self.folder_button.clicked.connect(self.selectFolder)

        # 创建新的文本框用于显示选择的路径
        self.selected_folder_input = QLineEdit()

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)

        layout = QVBoxLayout()
        layout.addWidget(self.hostname_label)
        layout.addWidget(self.hostname_input)
        layout.addWidget(self.username_label)
        layout.addWidget(self.username_input)
        layout.addWidget(self.password_label)
        layout.addWidget(self.password_input)
        layout.addWidget(self.folder_button)
        # 添加新的文本框到UI中
        layout.addWidget(self.selected_folder_input)
        layout.addWidget(button_box)
        self.setLayout(layout)

        self.folder_path = None  # 添加folder_path属性

    def selectFolder(self):
        folder = QFileDialog.getExistingDirectory(self, 'Select VASP File', '/')
        if folder:
            self.folder_path = folder
            # 将选择的路径显示在新的文本框中
            self.selected_folder_input.setText(folder)

    def getServerInfo(self):
        self.server_info = {
            'hostname': self.hostname_input.text(),
            'username': self.username_input.text(),
            'password': self.password_input.text(),
            'folder': self.folder_path,
        }
        return self.server_info

### 从CIF文件导入UI类
class CIFImportDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('From CIF')
        self.setModal(True)
        self.filepath = None

        self.filepath_label = QLabel('CIF Path:')
        self.filepath_input = QLineEdit()
        self.btn_browse = QPushButton('browse', self)
        self.btn_browse.clicked.connect(self.browseCIF)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)

        layout = QVBoxLayout()
        layout.addWidget(self.filepath_label)
        layout.addWidget(self.filepath_input)
        layout.addWidget(self.btn_browse)
        layout.addWidget(button_box)

        self.setLayout(layout)

    def browseCIF(self):
        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile)
        file_dialog.setNameFilter('CIF files (*.cif)')
        if file_dialog.exec_():
            self.filepath_input.setText(file_dialog.selectedFiles()[0])

    def getAtoms(self):
        self.filepath = self.filepath_input.text()
        try:
            atoms = read(self.filepath)
            return atoms
        except Exception as e:
            error_dialog = QErrorMessage(self)
            error_dialog.showMessage(f'Could not improt CIF file：{str(e)}')
            return None

### VASP output ui 类
class SSHReader(QThread):
    output_signal = pyqtSignal(str)

    def __init__(self, channel):
        super().__init__()
        self.channel = channel

    def run(self):
        while not self.channel.exit_status_ready():
            if self.channel.recv_ready():
                chunk = self.channel.recv(4096).decode(errors='replace')
                self.output_signal.emit(chunk)
                # yield to the event loop to allow the main thread to handle the data
                self.msleep(1)

         # after the process exited, there still can be some data to read
        while self.channel.recv_ready():
            chunk = self.channel.recv(4096).decode(errors='replace')
            self.output_signal.emit(chunk)
            self.msleep(1)
class VASPSyntaxHighlighter(QSyntaxHighlighter):
    def highlightBlock(self, text):
        yellow_format = QTextCharFormat()
        yellow_format.setForeground(QColor(255, 255, 0))

        light_blue_format = QTextCharFormat()
        light_blue_format.setForeground(QColor(173, 216, 230))

        green_format = QTextCharFormat()
        green_format.setForeground(QColor(0, 255, 0))

        if text.strip().startswith('DAV'):
            self.setFormat(0, len(text), yellow_format)
            dav_index = text.index('DAV')
            self.setFormat(dav_index, 3, light_blue_format)
        elif "N       E                     dE             d eps       ncg     rms          rms(c)" in text:
            self.setFormat(0, len(text), green_format)
class VASPOutputWidget(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()

        # Create a scroll area
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)

        # Create the text edit widget
        self.text_area = QTextEdit()
        
        self.highlighter = VASPSyntaxHighlighter(self.text_area.document())
        
        self.text_area.setReadOnly(True)
        self.text_area.setLineWrapMode(QTextEdit.NoWrap)  # Disable line wrapping

        # Set the font
        font = QFont("Courier New", 12)
        self.text_area.setFont(font)

        # Set the background color
        palette = self.text_area.palette()
        palette.setColor(QPalette.Base, QColor(50, 50, 50))
        palette.setColor(QPalette.Text, QColor(200, 200, 200))
        self.text_area.setPalette(palette)

        # Set the text edit widget as the scroll area's widget
        scroll_area.setWidget(self.text_area)

        layout.addWidget(scroll_area)
        self.setWindowTitle('VASP Output')
        self.setLayout(layout)
        self.resize(1280, 720)

    def append_text(self, text):
        self.text_area.moveCursor(QtGui.QTextCursor.End)  # move cursor to end
        self.text_area.insertPlainText(text)
        self.text_area.moveCursor(QtGui.QTextCursor.End)

### VASP 设置参数UI类
class VaspParamDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('set VASP parameter')
        self.setModal(True)

        self.xc_label = QLabel('xc:')
        self.xc_input = QLineEdit()

        self.kpts_label = QLabel('kpts:')
        self.kpts_input = QLineEdit()

        self.encut_label = QLabel('encut:')
        self.encut_input = QLineEdit()

        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)

        layout = QVBoxLayout()
        layout.addWidget(self.xc_label)
        layout.addWidget(self.xc_input)
        layout.addWidget(self.kpts_label)
        layout.addWidget(self.kpts_input)
        layout.addWidget(self.encut_label)
        layout.addWidget(self.encut_input)
        layout.addWidget(button_box)
        self.setLayout(layout)

    def getParams(self):
        xc = self.xc_input.text()
        kpts = eval(self.kpts_input.text())
        encut = int(self.encut_input.text())
        print()
        return {'xc': xc, 'kpts': kpts, 'encut': encut , 'command':f'""'} 

### 内核主控制模块
class AtomCalculator:
    def __init__(self):
        self.atoms_list = []
        self.atoms = None
        self.calculator_name = None
        self.optimizer_name = None
        self.energy = None
        self.forces = None

    def addAtom(self, element, x, y, z):
        atom_info = {'symbol': element, 'position': (x, y, z)}
        self.atoms_list.append(atom_info)
        self.updateStructure()

    def deleteAtom(self):
        if self.atoms_list:
            self.atoms_list.pop()
            self.updateStructure()

    def updateStructure(self):
        symbols = [atom['symbol'] for atom in self.atoms_list]
        positions = [atom['position'] for atom in self.atoms_list]
        self.atoms = Atoms(symbols=symbols, positions=positions)

    def setCalculator(self, calculator_name):
        self.calculator_name = calculator_name

    def setOptimizer(self, optimizer_name):
        self.optimizer_name = optimizer_name

    def startCalculation(self):
        if self.atoms is not None:
            if self.calculator_name == "emt":
                emt = EMT()
                self.atoms.set_calculator(emt)
                self.energy = self.atoms.get_potential_energy()
                self.forces = self.atoms.get_forces()

            elif self.calculator_name == "vasp":
                vasp_params = {
                    'xc': 'PBE',
                    'kpts': (4, 4, 4),
                    'encut': 400,
                }

                vasp_calculator = Vasp(txt = "vasp.out",**vasp_params)
                self.atoms.set_calculator(vasp_calculator)

                self.energy = self.atoms.get_potential_energy()
                self.forces = self.atoms.get_forces()

            elif self.calculator_name == "gaussian":
                pass

            else:
                self.energy = None
                self.forces = None

        else:
            self.energy = None
            self.forces = None

    @timer
    def run_vasp_calculation(self, ssh, file_path):
        try:
            # Initialize SSH client

            # Create a folder for VASP calculation
            print("Creat Work Folder")
            _, stdout, _ = ssh.exec_command('mkdir -p VASP_calculation')

            # Confirm folder creation
            stdout.channel.recv_exit_status()

            _, stdout, _ = ssh.exec_command('rm -rf ./VASP_calculation/*')

            # Confirm folder creation
            stdout.channel.recv_exit_status()

            # Transfer files to the remote server
            sftp = ssh.open_sftp()
            local_files = ['INCAR', 'POSCAR', 'POTCAR', 'KPOINTS']
            local_files_with_path = [os.path.join(file_path, file) for file in local_files]  # 添加文件夹路径
            remote_path = 'VASP_calculation/'

            # 获取总共需要传输的文件数量
            total_files = len(local_files_with_path)

            # 创建进度对话框
            progress_dialog = QProgressDialog("Transferring Files...", "Cancel", 0, total_files)
            progress_dialog.setWindowTitle("File Transfer")
            progress_dialog.setWindowModality(Qt.WindowModal)
            progress_dialog.setAutoReset(False)
            progress_dialog.setAutoClose(False)

            # 创建并启动进度更新器
            progress_updater = ProgressUpdater(total_files)
            progress_updater.progress_signal.connect(progress_dialog.setValue)
            progress_updater_thread = QThread()
            progress_updater.moveToThread(progress_updater_thread)
            progress_updater_thread.started.connect(progress_updater.update_progress)
            progress_updater_thread.start()

            # 运行文件传输
            for index, file in enumerate(local_files_with_path):
                sftp.put(file, os.path.join(remote_path, os.path.basename(file)))

            # 等待进度更新器完成工作
            progress_updater_thread.quit()
            progress_updater_thread.wait()
            progress_dialog.close()
            sftp.close()

            # Execute VASP calculation
            transport = ssh.get_transport()
            channel = transport.open_session()
            channel.exec_command(f'cd {remote_path} && vasp_std')

            widget = VASPOutputWidget()
            thread = SSHReader(channel)
            thread.output_signal.connect(widget.append_text)

            thread.start()
            widget.show()

            # Wait for the calculation
            while thread.isRunning():
                QApplication.processEvents()

            # Transfer results back to local Windows folder
            local_result_folder = 'Results'
            os.makedirs(local_result_folder, exist_ok=True)
            sftp = ssh.open_sftp()
            remote_files = sftp.listdir(remote_path)

            # 重置进度对话框，用于显示下载结果的进度
            progress_dialog.reset()
            progress_dialog.setLabelText("Downloading Results...")

            # 获取总共需要下载的文件数量
            total_files = len(remote_files)

            for index, file in enumerate(remote_files):
                local_file_path = os.path.join(local_result_folder, file)
                remote_file_path = os.path.join(remote_path, file)
                sftp.get(remote_file_path, local_file_path)

                # 更新进度条
                progress_dialog.setValue(index + 1)
                if progress_dialog.wasCanceled():
                    break

            progress_dialog.close()
            sftp.close()

            print("VASP calculation completed successfully.")
            ssh.close()

        except Exception as e:
            print("An error occurred:", e)


if __name__ == '__main__':

    app = QApplication(sys.argv)
    window = ASE_ui()
    window.show()
    sys.exit(app.exec_())    