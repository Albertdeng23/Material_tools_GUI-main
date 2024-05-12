import sys
from PyQt5.QtWidgets import *
from PyQt5.QtWidgets import QMessageBox
from chgnet.model.model import CHGNet
from pymatgen.core import Structure
import traceback
from chgnet.model.model import CHGNet
from chgnet.model.dynamics import MolecularDynamics
from pymatgen.core import Structure
import warnings
from PyQt5.QtCore import QThread, pyqtSignal
import time
# 启动PyQt应用
app = QApplication(sys.argv)

# AI推断
# 用于执行Direct Inference的线程
class DirectInferenceThread(QThread):
    result_signal = pyqtSignal(dict, float)

    def __init__(self, cif_path):
        QThread.__init__(self)
        self.cif_path = cif_path

    def run(self):
        start_time = time.time()
        chgnet = CHGNet.load()
        structure = Structure.from_file(self.cif_path)
        prediction = chgnet.predict_structure(structure)
        end_time = time.time()
        elapsed_time = end_time - start_time
        self.result_signal.emit(prediction, elapsed_time)
# 直接推断结果的窗口
class windowDirectInference(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Direct Inference Results")
        self.layout = QVBoxLayout(self)
        
        # 使用QTextEdit展示结果
        self.result_text = QTextEdit(self)
        self.result_text.setReadOnly(True)
        self.layout.addWidget(self.result_text)

    def displayResult(self, prediction, elapsed_time):
        result_str = f"Calculation Time: {elapsed_time:.2f} seconds\n\n"
        for key, unit in [
            ("energy", "eV/atom"),
            ("forces", "eV/A"),
            ("stress", "GPa"),
            ("magmom", "mu_B"),
        ]:
            result_str += f"CHGNet-predicted {key} ({unit}):\n{prediction[key[0]]}\n\n"
        self.result_text.setText(result_str)

# 分子动力学参数设置窗口
class WindowSetMolecularDynamics(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        
        self.setWindowTitle("Set Molecular Dynamics Parameters")
        
        self.form_layout = QFormLayout()
        
        # 创建输入字段
        self.ensemble_input = QLineEdit(self)
        self.ensemble_input.setText("nvt")
        self.temperature_input = QLineEdit(self)
        self.temperature_input.setText("1000")
        self.timestep_input = QLineEdit(self)
        self.timestep_input.setText("2")
        self.traj_input = QLineEdit(self)
        self.traj_input.setText("md_out.traj")
        self.logfile_input = QLineEdit(self)
        self.logfile_input.setText("md_out.log")
        self.loginterval_input = QLineEdit(self)
        self.loginterval_input.setText("100")
        self.steps_input = QLineEdit(self)
        self.steps_input.setText("50")
        self.device_combo = QComboBox(self)
        self.device_combo.addItem("cpu", "cpu")
        self.device_combo.addItem("cuda", "cuda")

        # 将输入字段加入布局
        self.form_layout.addRow("Ensemble", self.ensemble_input)
        self.form_layout.addRow("Temperature (K)", self.temperature_input)
        self.form_layout.addRow("Time Step (fs)", self.timestep_input)
        self.form_layout.addRow("Trajectory File", self.traj_input)
        self.form_layout.addRow("Log File", self.logfile_input)
        self.form_layout.addRow("Log Interval", self.loginterval_input)
        self.form_layout.addRow("Number of Steps", self.steps_input)
        self.form_layout.addRow("Device", self.device_combo)
        self.setLayout(self.form_layout)

        # 添加OK 和 Cancel按钮
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel, 
            parent=self
        )
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self.form_layout.addWidget(self.buttons)

    def getParameters(self):
        return {
            "ensemble": self.ensemble_input.text(),
            "temperature": int(self.temperature_input.text()),
            "timestep": float(self.timestep_input.text()),
            "trajectory": self.traj_input.text(),
            "logfile": self.logfile_input.text(),
            "loginterval": int(self.loginterval_input.text()),
            "steps": int(self.steps_input.text()),
            "device": self.device_combo.currentText()
        }



# 主窗口
class CHGnetApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = "CHGnet GUI"
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle(self.title)
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)
        
        # 文件选择按钮和标签
        self.file_button = QPushButton("Load CIF file")
        self.file_button.clicked.connect(self.openFileNameDialog)
        self.file_label = QLabel("Please load a CIF file to proceed")
        layout.addWidget(self.file_button)
        layout.addWidget(self.file_label)

        # 结果显示区域
        self.result_text = QTextEdit(central_widget)
        self.result_text.setReadOnly(True)
        layout.addWidget(self.result_text)
        
        # 预测类型选择下拉菜单和执行按钮
        self.prediction_type_combo = QComboBox(central_widget)
        self.prediction_type_combo.addItem("Direct Inference")
        self.prediction_type_combo.addItem("Molecular Dynamics")
        self.prediction_type_combo.addItem("Structure Optimization")
        layout.addWidget(self.prediction_type_combo)

        self.run_button = QPushButton("Run CHGnet")
        self.run_button.clicked.connect(self.runCHGnet)
        layout.addWidget(self.run_button)

    def openFileNameDialog(self):
        # 打开文件对话框
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getOpenFileName(self, "Load CIF File", "", 
                                                  "CIF Files (*.cif);;All Files (*)",
                                                  options=options)
        if fileName:
            self.file_label.setText(f"Loaded file: {fileName}")

    ### 主函数
    def runCHGnet(self):
        loaded_file = self.file_label.text().replace('Loaded file: ', '').strip()
        
        # 检查用户是否选择了预测类型和加载了文件
        if 'Loaded file' not in self.file_label.text() or not loaded_file:
            QMessageBox.information(self, 'Load File', 'Please load a CIF file before running the model.')
            return

        prediction_type = self.prediction_type_combo.currentText()
        if prediction_type == 'Direct Inference':
            self.performDirectInference(loaded_file)

        elif prediction_type == 'Molecular Dynamics':
            md_params_dialog = WindowSetMolecularDynamics(self)
            if md_params_dialog.exec_():
                md_params = md_params_dialog.getParameters()
                self.performMolecularDynamics(loaded_file, md_params)

        elif prediction_type == 'Structure Optimization':
            # 将来实现结构优化功能
            pass
        else:
            QMessageBox.information(self, 'Invalid Operation', 'Please select a valid operation to perform.')

    ## 直接推理
    def performDirectInference(self, cif_path):
        # 创建一个新窗口
        self.direct_inference_window = windowDirectInference()
        self.direct_inference_window.show()

        # 启动线程进行计算
        self.thread = DirectInferenceThread(cif_path)
        self.thread.result_signal.connect(self.direct_inference_window.displayResult)
        self.thread.start()

    ### 分子动力学模拟
    def performMolecularDynamics(self, file_path, parameters):
        try:
            # 设置警告过滤
            warnings.filterwarnings("ignore", module="pymatgen")
            warnings.filterwarnings("ignore", module="ase")

            # 加载CHGNet模型和结构
            chgnet = CHGNet.load()
            structure = Structure.from_file(file_path)

            # 创建分子动力学模拟对象
            md = MolecularDynamics(
                atoms=structure,
                model=chgnet,
                ensemble=parameters["ensemble"],
                temperature=parameters["temperature"],
                timestep=parameters["timestep"],
                trajectory=parameters["trajectory"],
                logfile=parameters["logfile"],
                loginterval=parameters["loginterval"],
                use_device=parameters["device"],
            )

            # 运行模拟
            md.run(parameters["steps"])
            
            # 显示提示消息
            self.result_text.setText("Molecular Dynamics simulation completed successfully!")

        except Exception as e:
            QMessageBox.critical(self, 'Simulation Error', 'An error occurred during the simulation:\n' + str(e) + '\n\n' + traceback.format_exc())
            self.result_text.setText('')


    
window = CHGnetApp()
window.show()
sys.exit(app.exec_())