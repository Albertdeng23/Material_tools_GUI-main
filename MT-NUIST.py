import paramiko
import sys
import os
import time
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from functools import partial
from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtGui import *
from PyQt5.QtCore import *

# ASE 相关库导入
from ase.visualize import view
from ase import Atoms
from ase.calculators.emt import EMT
from ase.calculators.vasp import Vasp
from ase.io import read

# 自定义模块导入 (假设这些文件在同一目录下)
# 注意：如果没有这些文件，运行时会报错，请确保文件存在或注释掉相关功能
try:
    from CIFDownload import CIFDownloader
    from DescriptorDesign import FileSelectionWidget
except ImportError:
    print("警告: 未找到 CIFDownload 或 DescriptorDesign 模块，相关功能将不可用。")
    CIFDownloader = None
    FileSelectionWidget = None


### TODO 待办事项列表
#   1. 使用 VASP 生成 POSCAR, INCAR 之类的输入文件逻辑完善

###----------------------------------------------------------------------------------------------------------####
##  工具函数部分
###----------------------------------------------------------------------------------------------------------####

def plot_band_structure(eigenval_file, figure):
    """
    读取 VASP 的 EIGENVAL 文件并在给定的 matplotlib figure 上绘制能带图。
    
    参数:
        eigenval_file (str): EIGENVAL 文件路径
        figure (matplotlib.figure.Figure): 用于绘图的画布对象
    """
    # 读取 EIGENVAL 文件内容
    with open(eigenval_file, 'r') as f:
        data = f.readlines()

    # 解析数据：获取费米能级 (第6行第1个数据)
    # 注意：VASP EIGENVAL 格式通常第6行包含 (电子数, k点数, 能带数) 或类似信息，具体取决于版本
    # 这里假设第6行包含费米能级信息，或者代码逻辑是针对特定格式编写的
    # 通常 EIGENVAL 第6行是: <num_electrons> <num_kpoints> <num_bands>
    # 真正的费米能级通常在 DOSCAR 或 OUTCAR 中，或者 EIGENVAL 某些特定格式下。
    # 此处代码逻辑假设第6行包含相关信息，需根据实际文件确认。
    try:
        efermi = float(data[5].split()[0]) 
    except (IndexError, ValueError):
        efermi = 0.0 # 如果读取失败，默认为0

    # 读取能带数据，跳过前8行头信息，只读取第2列(能量值)
    # 减去费米能级以将 0 eV 对齐到费米面
    eigenval_raw = np.loadtxt(eigenval_file, skiprows=8, usecols=(1,)) - efermi

    # 绘制能带图
    nkpoints = len(eigenval_raw)
    kpoints = np.linspace(0, 1, nkpoints) # 生成归一化的 k 点路径

    figure.clear() # 清除旧图
    ax = figure.add_subplot(111)
    ax.plot(kpoints, eigenval_raw, color='b')

    # 设置坐标轴和标签
    ax.axhline(y=0, color='k', ls='--') # y=0 参考线
    ax.set_xlabel('k-points (Normalized)')
    ax.set_ylabel('Energy (eV)')

    # 显示高对称点竖线 (假设每 0.25 为一个高对称点区间)
    for i in range(1, nkpoints - 1):
        if kpoints[i] % 0.25 == 0:
            ax.axvline(x=kpoints[i], color='k', ls='--')

    # 用红线强调费米能级
    ax.axhline(y=0, color='r', linewidth=0.5, ls='-')

    # 添加标题
    ax.set_title('Band Structure')

###----------------------------------------------------------------------------------------------------------####
##  UI 组件类定义
###----------------------------------------------------------------------------------------------------------####

class Band_Structure(QMainWindow):
    """
    能带结构可视化窗口类。
    包含文件选择、绘图显示和图片保存功能。
    """
    def __init__(self):
        super().__init__()
        self.setWindowTitle("能带图绘制")
        self.setGeometry(100, 100, 600, 500) # 调整了初始大小

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.create_main_layout()

    def create_main_layout(self):
        # 设置窗口样式表 (CSS 风格)
        self.setStyleSheet("""
            QWidget { font-size: 14px; }
            QPushButton {
                background-color: #5AA469; color: white;
                border-style: outset; border-width: 2px; border-radius: 10px;
                border-color: beige; font: bold 14px; padding: 6px; margin: 6px;
            }
            QPushButton:hover { background-color: #66c280; }
            QLabel { min-height: 3em; align: center; font: 14px; }
            QScrollArea { border: 0; }
        """)

        main_widget = QWidget(self)
        self.setCentralWidget(main_widget)
        main_layout = QVBoxLayout(main_widget)
        main_layout.setContentsMargins(20, 20, 20, 20)
        main_layout.setSpacing(10)

        # 创建滚动区域以容纳绘图
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        main_layout.addWidget(scroll_area)

        scroll_widget = QWidget()
        scroll_area.setWidget(scroll_widget)

        scroll_layout = QVBoxLayout(scroll_widget)
        self.image_label = QLabel()
        self.image_label.setAlignment(Qt.AlignCenter)
        self.image_label.setMinimumSize(400, 300)
        scroll_layout.addWidget(self.image_label)

        # 按钮布局
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

        # 信息提示标签
        self.label = QLabel()
        self.label.setAlignment(Qt.AlignCenter)
        scroll_layout.addWidget(self.label)

    def select_file(self, file_type):
        """打开文件选择对话框"""
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(self, '选择文件')
        if file_path:
            setattr(self, file_type, file_path) # 动态设置属性，如 self.EIGENVAL
            self.label.setText(f"{file_type} 文件已选择：{file_path}")

    def plot_band_structure(self):
        """调用绘图函数并更新界面"""
        if not hasattr(self, 'EIGENVAL'):
            self.label.setText("请先选择 EIGENVAL 文件！")
            return
            
        eigenval_file = self.EIGENVAL
        try:
            plot_band_structure(eigenval_file, self.figure)
            self.update_image_label()
            self.label.setText("绘制成功")
        except Exception as e:
            self.label.setText(f"绘制失败: {str(e)}")

    def update_image_label(self):
        """将 matplotlib 的 canvas 转换为 QPixmap 显示在 Label 上"""
        self.canvas.draw()
        w, h = self.canvas.get_width_height()
        image = QImage(self.canvas.buffer_rgba(), w, h, QImage.Format_RGBA8888)
        self.image_label.setPixmap(QPixmap.fromImage(image))

    def save_figure(self):
        """保存当前图像"""
        save_dialog = QFileDialog()
        save_path, _ = save_dialog.getSaveFileName(self, '保存图片', '', 'PNG (*.png);;JPEG (*.jpg *.jpeg)')
        if save_path:
            self.figure.savefig(save_path)
            self.label.setText(f"图片已保存：{save_path}")


class Ui_AseAtomInput(object):
    """
    主界面的 UI 布局定义类。
    通常由 Qt Designer 生成，这里包含了手写的布局代码。
    """
    def setupUi(self, AseAtomInput):
        AseAtomInput.setObjectName("AseAtomInput")
        AseAtomInput.resize(686, 418)
        # 定义全局样式
        styleSheet = """
            QWidget { font-size: 14px; }
            QPushButton {
                background-color: #5AA469; color: white;
                border-style: outset; border-width: 2px; border-radius: 10px;
                border-color: beige; font: bold 14px; padding: 6px; margin: 1px;
            }
            QPushButton:hover { background-color: #66c280; }
            QLineEdit { font: 14px Arial; border: 1px solid #ccc; padding: 3px; margin: 4px; }
            QLabel { font: bold italic; color: #555; }
            QGraphicsView { border-style: outset; border-width: 2px; border-radius: 5px; border-color: #ccc; }
            QMenuBar { background-color: #5AA469; color: white; }
            QMenuBar::item { background-color: #5AA469; color: white; }
            QMenuBar::item::selected { background-color: #66c280; }
            QMenu { background-color: white; color: #5AA469; }
            QMenu::item::selected { background-color: #66c280; color: white; }
        """
        AseAtomInput.setStyleSheet(styleSheet)
        
        self.centralwidget = QtWidgets.QWidget(AseAtomInput)
        self.centralwidget.setObjectName("centralwidget")
        
        # 右侧图形视图区域 (预留)
        self.graphicsView = QtWidgets.QGraphicsView(self.centralwidget)
        self.graphicsView.setGeometry(QtCore.QRect(380, 40, 261, 241))
        self.graphicsView.setObjectName("graphicsView")
        
        # 底部按钮区域
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
        
        # 左侧输入区域
        self.layoutWidget1 = QtWidgets.QWidget(self.centralwidget)
        self.layoutWidget1.setGeometry(QtCore.QRect(10, 10, 351, 287))
        self.layoutWidget1.setObjectName("layoutWidget1")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.layoutWidget1)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        
        # 原子符号输入
        self.label = QtWidgets.QLabel(self.layoutWidget1)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.AtomInPut = QtWidgets.QLineEdit(self.layoutWidget1)
        self.AtomInPut.setObjectName("AtomInPut")
        self.verticalLayout.addWidget(self.AtomInPut)
        
        # 坐标输入 X
        self.label_2 = QtWidgets.QLabel(self.layoutWidget1)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.X_Input = QtWidgets.QLineEdit(self.layoutWidget1)
        self.X_Input.setText("")
        self.X_Input.setObjectName("X_Input")
        self.verticalLayout.addWidget(self.X_Input)
        
        # 坐标输入 Y
        self.label_3 = QtWidgets.QLabel(self.layoutWidget1)
        self.label_3.setObjectName("label_3")
        self.verticalLayout.addWidget(self.label_3)
        self.Y_Input = QtWidgets.QLineEdit(self.layoutWidget1)
        self.Y_Input.setText("")
        self.Y_Input.setObjectName("Y_Input")
        self.verticalLayout.addWidget(self.Y_Input)
        
        # 坐标输入 Z
        self.label_4 = QtWidgets.QLabel(self.layoutWidget1)
        self.label_4.setObjectName("label_4")
        self.verticalLayout.addWidget(self.label_4)
        self.Z_Input = QtWidgets.QLineEdit(self.layoutWidget1)
        self.Z_Input.setText("")
        self.Z_Input.setObjectName("Z_Input")
        self.verticalLayout.addWidget(self.Z_Input)
        
        # 计算器名称输入
        self.label_5 = QtWidgets.QLabel(self.layoutWidget1)
        self.label_5.setObjectName("label_5")
        self.verticalLayout.addWidget(self.label_5)
        self.Caluator_Input = QtWidgets.QLineEdit(self.layoutWidget1)
        self.Caluator_Input.setText("")
        self.Caluator_Input.setObjectName("Caluator_Input")
        self.verticalLayout.addWidget(self.Caluator_Input)
        
        AseAtomInput.setCentralWidget(self.centralwidget)
        
        # 菜单栏设置
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
        
        # 菜单动作
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
        self.AddAtom_Buttom.setText(_translate("AseAtomInput", "添加原子"))
        self.DelAtom.setText(_translate("AseAtomInput", "删除原子"))
        self.ViewStructure_Button.setText(_translate("AseAtomInput", "查看结构"))
        self.SetCaluator_Button.setText(_translate("AseAtomInput", "设置计算器"))
        self.Process_Button.setText(_translate("AseAtomInput", "开始计算"))
        self.Cancel_Buttom.setText(_translate("AseAtomInput", "取消"))
        self.label.setText(_translate("AseAtomInput", "原子符号"))
        self.label_2.setText(_translate("AseAtomInput", "X 坐标"))
        self.label_3.setText(_translate("AseAtomInput", "Y 坐标"))
        self.label_4.setText(_translate("AseAtomInput", "Z 坐标"))
        self.label_5.setText(_translate("AseAtomInput", "计算器类型"))
        self.menu.setTitle(_translate("AseAtomInput", "文件"))
        self.menuConnect_Server.setTitle(_translate("AseAtomInput", "VASP"))
        self.menuAbout.setTitle(_translate("AseAtomInput", "关于"))
        self.menuPlot.setTitle(_translate("AseAtomInput", "绘图"))
        self.menuMaterials_Project.setTitle(_translate("AseAtomInput", "Materials Project"))
        self.actionFrom.setText(_translate("AseAtomInput", "导入 CIF"))
        self.actionServer.setText(_translate("AseAtomInput", "连接服务器"))
        self.actionVision.setText(_translate("AseAtomInput", "版本信息"))
        self.actionAbout_Author.setText(_translate("AseAtomInput", "关于作者"))
        self.actionBand_Structure.setText(_translate("AseAtomInput", "能带结构"))
        self.actionCIF_Download.setText(_translate("AseAtomInput", "CIF 下载"))
        self.actionDescriptor_Design.setText(_translate("AseAtomInput", "描述符设计"))

class MyApp(QWidget):
    """版本信息显示窗口"""
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle('软件信息')
        self.setGeometry(300, 300, 500, 350)
        self.setMinimumSize(500, 350)
        self.setStyleSheet("background-color: #F5DEB3;")

        layout = QVBoxLayout()

        version_label = QLabel('Version 1.2')
        version_label.setFont(QFont('Arial', 20, QFont.Bold))
        version_label.setStyleSheet("color: #FFFFFF; margin-bottom: 20px;")
        version_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(version_label)

        github_label = QLabel('<a href="https://github.com/Albertdeng23/Material_tools_GUI-main">GitHub: https://github.com/Albertdeng23/Material_tools_GUI-main</a>')
        github_label.setFont(QFont('Arial', 16))
        github_label.setStyleSheet("color: #FFFFFF;")
        github_label.setAutoFillBackground(True)
        p = github_label.palette()
        p.setColor(github_label.backgroundRole(), Qt.transparent)
        github_label.setPalette(p)
        github_label.setOpenExternalLinks(True)
        github_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(github_label)

        release_label = QLabel('<a href="https://github.com/Albertdeng23/Material_tools_GUI-main/releases">Release: https://github.com/Albertdeng23/Material_tools_GUI-main/releases</a>')
        release_label.setFont(QFont('Arial', 16))
        release_label.setStyleSheet("color: #FFFFFF; margin-top: 10px;")
        release_label.setAutoFillBackground(True)
        p = release_label.palette()
        p.setColor(release_label.backgroundRole(), Qt.transparent)
        release_label.setPalette(p)
        release_label.setOpenExternalLinks(True)
        release_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(release_label)

        self.setLayout(layout)

class AboutAuthorDialog(QDialog):
    """关于作者对话框"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('关于作者')
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

###----------------------------------------------------------------------------------------------------------####
##  多线程与进度条工具
###----------------------------------------------------------------------------------------------------------####

class ProgressUpdater(QObject):
    """用于在后台线程中发送进度信号的类"""
    progress_signal = pyqtSignal(int)

    def __init__(self, total_files):
        super().__init__()
        self.total_files = total_files

    def update_progress(self):
        """模拟文件传输进度"""
        for index in range(1, self.total_files + 1):
            self.progress_signal.emit(index)
            time.sleep(0.01)  # 模拟传输延迟

class ThreadedProgressUpdater(QRunnable):
    """QRunnable 包装器，用于在 QThreadPool 中运行"""
    def __init__(self, progress_updater):
        super().__init__()
        self.progress_updater = progress_updater

    def run(self):
        self.progress_updater.update_progress()

###----------------------------------------------------------------------------------------------------------####
##  主程序逻辑类 (ASE_ui)
###----------------------------------------------------------------------------------------------------------####

class ASE_ui(Ui_AseAtomInput, QMainWindow): 
    """
    主应用程序窗口类，继承自 UI 类和 QMainWindow。
    负责连接 UI 信号与后端逻辑。
    """
    def __init__(self) -> None:
        super().__init__()
        self.atom_calculator = AtomCalculator()
        self.myAppInstance = MyApp()
        self.about_author_dialog = AboutAuthorDialog(self)
        self.initUI()

    def initUI(self):
        self.setupUi(self)
        self.retranslateUi(self)

        # 信号槽连接
        self.DelAtom.clicked.connect(self.deleteAtom)
        self.AddAtom_Buttom.clicked.connect(self.addAtom)
        self.ViewStructure_Button.clicked.connect(self.showStructure)
        self.SetCaluator_Button.clicked.connect(self.setCalculator)
        self.Process_Button.clicked.connect(self.startCalculation)
        
        # 菜单动作连接
        self.actionFrom.triggered.connect(self.importCIF)
        self.actionServer.triggered.connect(self.connectRemoteServer)
        self.actionVision.triggered.connect(self.versionButton)
        self.actionAbout_Author.triggered.connect(self.showAboutAuthorDialog)
        self.actionBand_Structure.triggered.connect(self.plot_Band)
        self.actionCIF_Download.triggered.connect(self.showCFIDownloader)
        self.actionDescriptor_Design.triggered.connect(self.showDescriptorDesign)

    def showAboutAuthorDialog(self):
        self.about_author_dialog.exec_()

    def addAtom(self):
        """从输入框读取数据并添加原子"""
        element = self.AtomInPut.text()
        try:
            x = float(self.X_Input.text())
            y = float(self.Y_Input.text())
            z = float(self.Z_Input.text())
            self.atom_calculator.addAtom(element, x, y, z)
            # 清空输入框
            self.AtomInPut.clear()
            self.X_Input.clear()
            self.Y_Input.clear()
            self.Z_Input.clear()
        except ValueError:
            QMessageBox.warning(self, "输入错误", "坐标必须是数字！")

    def versionButton(self):
        self.myAppInstance.show()

    def deleteAtom(self):
        """删除上一个添加的原子"""
        self.atom_calculator.deleteAtom()

    def finishInput(self):
        """结束输入 (关闭窗口)"""
        self.close()

    def showStructure(self):
        """调用 ASE 的 view 函数显示原子结构"""
        if self.atom_calculator.atoms is not None:
            view(self.atom_calculator.atoms)
        else:
            QMessageBox.information(self, "提示", "当前没有原子结构可显示。")

    def setCalculator(self):
        """设置计算器类型"""
        calculator_name = self.Caluator_Input.text()
        self.atom_calculator.setCalculator(calculator_name)
        if calculator_name.lower() == "vasp":
            self.showVaspParamDialog()
        else:
            # 简单的反馈
            show_Calculator = QErrorMessage(self)
            show_Calculator.showMessage(f"已设置计算器: {calculator_name}")
    
    def plot_Band(self):
        """打开能带图绘制窗口"""
        self.Band_structure = Band_Structure()
        self.Band_structure.show()

    def showVaspParamDialog(self):
        """显示 VASP 参数设置对话框"""
        self.dialog = VaspParamDialog(self)
        if self.dialog.exec_() == QDialog.Accepted:
            self.vasp_params = self.dialog.getParams()
            show_Calculator = QErrorMessage(self)
            show_Calculator.showMessage(f"计算器: VASP, 参数: {self.vasp_params}")

    def setOptimizer(self):
        """设置优化器 (预留功能)"""
        optimizer_name = self.optimizer_input.text()
        self.atom_calculator.setOptimizer(optimizer_name)
        show_Optimizer = QErrorMessage(self)
        show_Optimizer.showMessage(f"优化器为: {optimizer_name}")

    def startCalculation(self):
        """开始计算"""
        self.atom_calculator.startCalculation()
        energy = self.atom_calculator.energy
        forces = self.atom_calculator.forces
        message = QErrorMessage(self)
        message.showMessage(f"能量: {energy}, 受力: {forces}")

    def importCIF(self):
        """导入 CIF 文件"""
        dialog = CIFImportDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            atoms = dialog.getAtoms()
            if atoms is not None:
                self.atom_calculator.atoms = atoms
                show_Structure = QErrorMessage(self)
                show_Structure.showMessage('成功导入原子信息')
    
    def connectRemoteServer(self):
        """连接远程服务器逻辑"""
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
                    print('连接成功')

                    # 运行 VASP 计算流程
                    use_times = self.atom_calculator.run_vasp_calculation(self.ssh , self.file_path)
                    QMessageBox.information(self, '成功', f'成功连接到 {self.hostname} \n 结果文件路径 : .\\Results \n 耗时: {use_times /60 : .2f} 分钟。')

                except paramiko.AuthenticationException:
                    QMessageBox.warning(self, '认证错误', '请检查您的账号和密码')
                    self.connectRemoteServer() # 重新尝试

                except Exception as e:
                    error_dialog = QErrorMessage(self)
                    error_dialog.showMessage(f'连接失败：{str(e)}')

    def showCFIDownloader(self):
        if CIFDownloader:
            self.CIFDownload = CIFDownloader()
            self.CIFDownload.show()
        else:
            QMessageBox.warning(self, "错误", "CIFDownloader 模块未加载")

    def showDescriptorDesign(self):
        if FileSelectionWidget:
            self.DescriptorDesigner = FileSelectionWidget()
            self.DescriptorDesigner.show()
        else:
            QMessageBox.warning(self, "错误", "DescriptorDesign 模块未加载")

###----------------------------------------------------------------------------------------------------------####
##  辅助对话框类
###----------------------------------------------------------------------------------------------------------####

class RemoteServerDialog(QDialog):
    """远程服务器连接信息输入对话框"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('连接远程服务器')
        self.setModal(True)
        self.server_info = None

        self.hostname_label = QLabel('主机名 (IP):')
        self.username_label = QLabel('账号:')
        self.password_label = QLabel('密码:')

        self.hostname_input = QLineEdit()
        self.username_input = QLineEdit()
        self.password_input = QLineEdit()
        self.password_input.setEchoMode(QLineEdit.Password)

        self.folder_button = QPushButton('选择本地 VASP 文件夹')
        self.folder_button.clicked.connect(self.selectFolder)

        # 显示选择的路径
        self.selected_folder_input = QLineEdit()
        self.selected_folder_input.setReadOnly(True)

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
        layout.addWidget(self.selected_folder_input)
        layout.addWidget(button_box)
        self.setLayout(layout)

        self.folder_path = None

    def selectFolder(self):
        folder = QFileDialog.getExistingDirectory(self, '选择 VASP 文件所在目录', '/')
        if folder:
            self.folder_path = folder
            self.selected_folder_input.setText(folder)

    def getServerInfo(self):
        self.server_info = {
            'hostname': self.hostname_input.text(),
            'username': self.username_input.text(),
            'password': self.password_input.text(),
            'folder': self.folder_path,
        }
        return self.server_info

class CIFImportDialog(QDialog):
    """CIF 文件导入对话框"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('从 CIF 导入')
        self.setModal(True)
        self.filepath = None

        self.filepath_label = QLabel('CIF 路径:')
        self.filepath_input = QLineEdit()
        self.btn_browse = QPushButton('浏览', self)
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
            error_dialog.showMessage(f'无法导入 CIF 文件：{str(e)}')
            return None

###----------------------------------------------------------------------------------------------------------####
##  VASP 输出显示与高亮
###----------------------------------------------------------------------------------------------------------####

class SSHReader(QThread):
    """后台线程：实时读取 SSH 通道的输出"""
    output_signal = pyqtSignal(str)

    def __init__(self, channel):
        super().__init__()
        self.channel = channel

    def run(self):
        # 当通道未关闭时循环读取
        while not self.channel.exit_status_ready():
            if self.channel.recv_ready():
                chunk = self.channel.recv(4096).decode(errors='replace')
                self.output_signal.emit(chunk)
                self.msleep(1) # 让出时间片

         # 进程退出后，读取剩余数据
        while self.channel.recv_ready():
            chunk = self.channel.recv(4096).decode(errors='replace')
            self.output_signal.emit(chunk)
            self.msleep(1)

class VASPSyntaxHighlighter(QSyntaxHighlighter):
    """VASP 输出日志的语法高亮器"""
    def highlightBlock(self, text):
        yellow_format = QTextCharFormat()
        yellow_format.setForeground(QColor(255, 255, 0))

        light_blue_format = QTextCharFormat()
        light_blue_format.setForeground(QColor(173, 216, 230))

        green_format = QTextCharFormat()
        green_format.setForeground(QColor(0, 255, 0))

        # 高亮 DAV 迭代行
        if text.strip().startswith('DAV'):
            self.setFormat(0, len(text), yellow_format)
            try:
                dav_index = text.index('DAV')
                self.setFormat(dav_index, 3, light_blue_format)
            except ValueError:
                pass
        # 高亮表头
        elif "N       E                     dE             d eps       ncg     rms          rms(c)" in text:
            self.setFormat(0, len(text), green_format)

class VASPOutputWidget(QWidget):
    """显示 VASP 实时输出的窗口"""
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)

        self.text_area = QTextEdit()
        self.highlighter = VASPSyntaxHighlighter(self.text_area.document())
        
        self.text_area.setReadOnly(True)
        self.text_area.setLineWrapMode(QTextEdit.NoWrap) # 禁止自动换行

        font = QFont("Courier New", 12)
        self.text_area.setFont(font)

        # 设置深色背景
        palette = self.text_area.palette()
        palette.setColor(QPalette.Base, QColor(50, 50, 50))
        palette.setColor(QPalette.Text, QColor(200, 200, 200))
        self.text_area.setPalette(palette)

        scroll_area.setWidget(self.text_area)

        layout.addWidget(scroll_area)
        self.setWindowTitle('VASP Output')
        self.setLayout(layout)
        self.resize(1280, 720)

    def append_text(self, text):
        self.text_area.moveCursor(QtGui.QTextCursor.End)
        self.text_area.insertPlainText(text)
        self.text_area.moveCursor(QtGui.QTextCursor.End)

class VaspParamDialog(QDialog):
    """VASP 参数设置对话框"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('设置 VASP 参数')
        self.setModal(True)

        self.xc_label = QLabel('交换关联泛函 (xc):')
        self.xc_input = QLineEdit()

        self.kpts_label = QLabel('K点设置 (kpts):')
        self.kpts_input = QLineEdit()

        self.encut_label = QLabel('截断能 (encut):')
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
        try:
            kpts = eval(self.kpts_input.text()) # 注意：eval 有安全风险，生产环境建议用 ast.literal_eval
            encut = int(self.encut_input.text())
        except:
            kpts = (1,1,1)
            encut = 400
        return {'xc': xc, 'kpts': kpts, 'encut': encut , 'command':f'""'} 

###----------------------------------------------------------------------------------------------------------####
##  核心逻辑类 (AtomCalculator)
###----------------------------------------------------------------------------------------------------------####

class AtomCalculator:
    """
    后端逻辑核心类。
    负责管理原子结构、设置计算器、以及执行远程 VASP 计算。
    """
    def __init__(self):
        self.atoms_list = []
        self.atoms = None
        self.calculator_name = None
        self.optimizer_name = None
        self.energy = None
        self.forces = None

    def addAtom(self, element, x, y, z):
        """添加原子到列表并更新结构"""
        atom_info = {'symbol': element, 'position': (x, y, z)}
        self.atoms_list.append(atom_info)
        self.updateStructure()

    def deleteAtom(self):
        """删除最后一个原子"""
        if self.atoms_list:
            self.atoms_list.pop()
            self.updateStructure()

    def updateStructure(self):
        """根据原子列表重建 ASE Atoms 对象"""
        symbols = [atom['symbol'] for atom in self.atoms_list]
        positions = [atom['position'] for atom in self.atoms_list]
        if symbols:
            self.atoms = Atoms(symbols=symbols, positions=positions)
        else:
            self.atoms = None

    def setCalculator(self, calculator_name):
        self.calculator_name = calculator_name

    def setOptimizer(self, optimizer_name):
        self.optimizer_name = optimizer_name

    def startCalculation(self):
        """本地计算逻辑 (EMT 或 VASP 接口)"""
        if self.atoms is not None:
            if self.calculator_name == "emt":
                emt = EMT()
                self.atoms.set_calculator(emt)
                self.energy = self.atoms.get_potential_energy()
                self.forces = self.atoms.get_forces()

            elif self.calculator_name == "vasp":
                # 这里的 VASP 是 ASE 的本地接口，需要本地安装 VASP
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
                pass # 待实现

            else:
                self.energy = None
                self.forces = None
        else:
            self.energy = None
            self.forces = None

    def run_vasp_calculation(self, ssh, file_path):
        """
        在远程服务器上运行 VASP 计算。
        
        步骤:
        1. 在服务器创建工作目录。
        2. 上传本地输入文件 (INCAR, POSCAR, POTCAR, KPOINTS)。
        3. 执行 vasp_std 命令。
        4. 实时显示输出。
        5. 计算完成后下载结果。
        """
        start_time = time.time()

        # 1. 创建远程工作目录
        print("正在创建工作目录...")
        _, stdout, _ = ssh.exec_command('mkdir -p VASP_calculation')
        stdout.channel.recv_exit_status() # 等待命令执行完毕

        # 清空目录 (慎用 rm -rf)
        _, stdout, _ = ssh.exec_command('rm -rf ./VASP_calculation/*')
        stdout.channel.recv_exit_status()

        # 2. 传输文件到远程服务器
        sftp = ssh.open_sftp()
        local_files = ['INCAR', 'POSCAR', 'POTCAR', 'KPOINTS']
        local_files_with_path = [os.path.join(file_path, file) for file in local_files]
        remote_path = 'VASP_calculation/'

        total_files = len(local_files_with_path)

        # 创建进度对话框
        progress_dialog = QProgressDialog("正在传输文件...", "取消", 0, total_files)
        progress_dialog.setWindowTitle("文件传输")
        progress_dialog.setWindowModality(Qt.WindowModal)
        progress_dialog.setAutoReset(False)
        progress_dialog.setAutoClose(False)

        # 启动进度更新线程
        progress_updater = ProgressUpdater(total_files)
        progress_updater.progress_signal.connect(progress_dialog.setValue)
        progress_updater_thread = QThread()
        progress_updater.moveToThread(progress_updater_thread)
        progress_updater_thread.started.connect(progress_updater.update_progress)
        progress_updater_thread.start()

        # 执行上传
        for index, file in enumerate(local_files_with_path):
            try:
                sftp.put(file, os.path.join(remote_path, os.path.basename(file)))
            except FileNotFoundError:
                print(f"错误: 本地文件未找到: {file}")

        # 清理进度线程
        progress_updater_thread.quit()
        progress_updater_thread.wait()
        progress_dialog.close()
        sftp.close()

        # 3. 执行 VASP 计算
        transport = ssh.get_transport()
        channel = transport.open_session()
        # 注意: 这里假设远程服务器环境变量中已有 vasp_std 命令
        channel.exec_command(f'cd {remote_path} && vasp_std')

        # 4. 实时显示输出
        widget = VASPOutputWidget()
        thread = SSHReader(channel)
        thread.output_signal.connect(widget.append_text)

        thread.start()
        widget.show()

        # 等待计算线程结束 (阻塞 UI 直到计算完成，实际应用中可能需要优化以免界面卡死)
        while thread.isRunning():
            QApplication.processEvents()

        # 5. 将结果传回本地 Windows 文件夹
        local_result_folder = 'Results'
        os.makedirs(local_result_folder, exist_ok=True)
        sftp = ssh.open_sftp()
        try:
            remote_files = sftp.listdir(remote_path)
        except IOError:
            remote_files = []

        # 重置进度对话框用于下载
        progress_dialog.reset()
        progress_dialog.setLabelText("正在下载结果...")
        progress_dialog.setRange(0, len(remote_files))
        progress_dialog.show()

        for index, file in enumerate(remote_files):
            local_file_path = os.path.join(local_result_folder, file)
            remote_file_path = os.path.join(remote_path, file)
            try:
                sftp.get(remote_file_path, local_file_path)
            except Exception as e:
                print(f"下载 {file} 失败: {e}")

            progress_dialog.setValue(index + 1)
            if progress_dialog.wasCanceled():
                break

        progress_dialog.close()
        sftp.close()
        ssh.close()
        
        print("VASP 计算成功完成。")
        end_time = time.time()
        return end_time - start_time


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = ASE_ui()
    window.show()
    sys.exit(app.exec_())
