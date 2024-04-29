import paramiko
import sys
# from Visualizer.demo import Ui_AseAtomInput
from PyQt5.QtWidgets import *
from PyQt5 import *
from ase.visualize import view
from ase import Atoms
from ase.calculators.emt import EMT
from ase.calculators.vasp import Vasp
from ase.io import read
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtGui import QFont, QColor, QPalette
from PyQt5.QtCore import Qt
import os

### TODO 
#   1.使用VASP生成POSCAR , INCAR 之类的输入文件

### mainUI类
class Ui_AseAtomInput(object):
    def setupUi(self, AseAtomInput):
        AseAtomInput.setObjectName("AseAtomInput")
        AseAtomInput.resize(686, 424)
        self.centralwidget = QtWidgets.QWidget(AseAtomInput)
        self.centralwidget.setObjectName("centralwidget")
        self.graphicsView = QtWidgets.QGraphicsView(self.centralwidget)
        self.graphicsView.setGeometry(QtCore.QRect(380, 40, 261, 241))
        self.graphicsView.setObjectName("graphicsView")
        self.widget = QtWidgets.QWidget(self.centralwidget)
        self.widget.setGeometry(QtCore.QRect(10, 330, 651, 30))
        self.widget.setObjectName("widget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.widget)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.AddAtom_Buttom = QtWidgets.QPushButton(self.widget)
        self.AddAtom_Buttom.setObjectName("AddAtom_Buttom")
        self.horizontalLayout.addWidget(self.AddAtom_Buttom)
        self.DelAtom = QtWidgets.QPushButton(self.widget)
        self.DelAtom.setObjectName("DelAtom")
        self.horizontalLayout.addWidget(self.DelAtom)
        self.ViewStructure_Button = QtWidgets.QPushButton(self.widget)
        self.ViewStructure_Button.setObjectName("ViewStructure_Button")
        self.horizontalLayout.addWidget(self.ViewStructure_Button)
        self.SetCaluator_Button = QtWidgets.QPushButton(self.widget)
        self.SetCaluator_Button.setObjectName("SetCaluator_Button")
        self.horizontalLayout.addWidget(self.SetCaluator_Button)
        self.Process_Button = QtWidgets.QPushButton(self.widget)
        self.Process_Button.setObjectName("Process_Button")
        self.horizontalLayout.addWidget(self.Process_Button)
        self.Cancel_Buttom = QtWidgets.QPushButton(self.widget)
        self.Cancel_Buttom.setObjectName("Cancel_Buttom")
        self.horizontalLayout.addWidget(self.Cancel_Buttom)
        self.widget1 = QtWidgets.QWidget(self.centralwidget)
        self.widget1.setGeometry(QtCore.QRect(10, 10, 351, 287))
        self.widget1.setObjectName("widget1")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.widget1)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtWidgets.QLabel(self.widget1)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.AtomInPut = QtWidgets.QLineEdit(self.widget1)
        self.AtomInPut.setObjectName("AtomInPut")
        self.verticalLayout.addWidget(self.AtomInPut)
        self.label_2 = QtWidgets.QLabel(self.widget1)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.X_Input = QtWidgets.QLineEdit(self.widget1)
        self.X_Input.setText("")
        self.X_Input.setObjectName("X_Input")
        self.verticalLayout.addWidget(self.X_Input)
        self.label_3 = QtWidgets.QLabel(self.widget1)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.verticalLayout.addWidget(self.label_3)
        self.Y_Input = QtWidgets.QLineEdit(self.widget1)
        self.Y_Input.setText("")
        self.Y_Input.setObjectName("Y_Input")
        self.verticalLayout.addWidget(self.Y_Input)
        self.label_4 = QtWidgets.QLabel(self.widget1)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_4.setFont(font)
        self.label_4.setObjectName("label_4")
        self.verticalLayout.addWidget(self.label_4)
        self.Z_Input = QtWidgets.QLineEdit(self.widget1)
        self.Z_Input.setText("")
        self.Z_Input.setObjectName("Z_Input")
        self.verticalLayout.addWidget(self.Z_Input)
        self.label_5 = QtWidgets.QLabel(self.widget1)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.verticalLayout.addWidget(self.label_5)
        self.Caluator_Input = QtWidgets.QLineEdit(self.widget1)
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
        self.menu.addAction(self.actionFrom)
        self.menuConnect_Server.addAction(self.actionServer)
        self.menuAbout.addAction(self.actionVision)
        self.menuAbout.addAction(self.actionAbout_Author)
        self.menubar.addAction(self.menu.menuAction())
        self.menubar.addAction(self.menuConnect_Server.menuAction())
        self.menubar.addAction(self.menuAbout.menuAction())

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
        self.menuConnect_Server.setTitle(_translate("AseAtomInput", "Connect Server"))
        self.menuAbout.setTitle(_translate("AseAtomInput", "About"))
        self.actionFrom.setText(_translate("AseAtomInput", "Import cif"))
        self.actionServer.setText(_translate("AseAtomInput", "Server"))
        self.actionVision.setText(_translate("AseAtomInput", "Version"))
        self.actionAbout_Author.setText(_translate("AseAtomInput", "About Author"))

### verison类
class MyApp(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
    
    def initUI(self):
        # 设置窗口基本属性
        self.setWindowTitle('Software Information')
        self.setGeometry(300, 300, 400, 300)
        self.setMinimumSize(400, 300)

        # 设置窗口背景色
        palette = QPalette()
        palette.setColor(QPalette.Window, QColor(240, 240, 240))
        self.setPalette(palette)

        # 设置布局
        layout = QVBoxLayout()

        # 设置并添加版本标签
        version_label = QLabel('Version 1.2')
        version_label.setFont(QFont('Arial', 18, QFont.Bold))
        version_label.setStyleSheet("color: blue; margin-bottom: 20px;")
        version_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(version_label)

        # 设置并添加GitHub标签
        github_label = QLabel('GitHub: https://example.com')
        github_label.setFont(QFont('Arial', 16))
        github_label.setStyleSheet("background-color: white; color: black; padding: 10px; border: 2px solid #0066CC;")
        github_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(github_label)

        # 设置并添加Release标签
        release_label = QLabel('Release: https://example.com/release')
        release_label.setFont(QFont('Arial', 16))
        release_label.setStyleSheet("background-color: white; color: black; padding: 10px; border: 2px solid #0066CC;")
        release_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(release_label)

        # 应用布局
        self.setLayout(layout)

## about author类


### 创建MainWindow
class ASE_ui(Ui_AseAtomInput,QMainWindow): 
    def __init__(self,) -> None:
        super().__init__()
        self.atom_calculator = AtomCalculator()
        self.myAppInstance = MyApp()
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
                hostname = server_info['hostname']
                username = server_info['username']
                password = server_info['password']


                try:
                    
                    self.atom_calculator.run_vasp_calculation(hostname , username , password)

                    QMessageBox.information(self, f'successful', f'successfully connect to {hostname}')

                except Exception as e:
                    error_dialog = QErrorMessage(self)
                    error_dialog.showMessage(f'connect faild：{str(e)}')

### DONE 链接远程服务器类，用于获取远程服务器的连接信息
class RemoteServerDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('connect Remote Server')
        self.setModal(True)
        self.server_info = None

        self.hostname_label = QLabel('hosthame:')
        self.username_label = QLabel('acount:')
        self.password_label = QLabel('passwd:')

        self.hostname_input = QLineEdit()
        self.username_input = QLineEdit()
        self.password_input = QLineEdit()
        self.password_input.setEchoMode(QLineEdit.Password)

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
        layout.addWidget(button_box)
        self.setLayout(layout)

    def getServerInfo(self):
        self.server_info = {
            'hostname': self.hostname_input.text(),
            'username': self.username_input.text(),
            'password': self.password_input.text(),
        }
        return self.server_info
### 从CIF文件导入类
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

### VASP 设置参数类
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


    def run_vasp_calculation(self , hostname , username , passwd):
        try:
            # 连接到linux服务器，基于paramiko
            ssh = paramiko.SSHClient()
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh.connect(hostname, username=username, password=passwd)

            # Create a folder for VASP calculation
            ssh.exec_command('mkdir -p VASP_calculation')
            
            # Transfer files to the remote server
            sftp = ssh.open_sftp()
            local_files = ['INCAR', 'POSCAR', 'POTCAR', 'KPOINTS']
            remote_path = 'VASP_calculation/'
            for file in local_files:
                sftp.put(file, os.path.join(remote_path, file))
            sftp.close()


            # Execute VASP calculation
            stdin, stdout, stderr = ssh.exec_command('cd {} && vasp_std'.format(remote_path))
            # Wait for the command to finish
            stdout.channel.recv_exit_status()

            # Transfer results back to local Windows folder
            local_result_folder = 'Results'
            os.makedirs(local_result_folder, exist_ok=True)
            for file in local_files:
                sftp.get(os.path.join(remote_path, file), os.path.join(local_result_folder, file))
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

    