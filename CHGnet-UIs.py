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
from PyQt5.QtCore import *

# 启动PyQt应用
app = QApplication(sys.argv)

## 分子动力学参数设置窗口
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
# 工作线程类
class WorkerThread(QThread):
    # 定义工作完成和错误信号
    done = pyqtSignal()
    error = pyqtSignal(Exception)

    def __init__(self, function, *args, **kwargs):
        super().__init__()
        self.function = function
        self.args = args
        self.kwargs = kwargs
      
    def run(self):
        try:
            self.function(*self.args, **self.kwargs)
            self.done.emit()
        except Exception as e:
            self.error.emit(e)


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

                # 显示进度条对话框
                self.progress_dialog = QProgressDialog("Running Molecular Dynamics simulation...", None, 0, 100, self)
                self.progress_dialog.setWindowModality(Qt.WindowModal)
                self.progress_dialog.setValue(0)
                self.progress_dialog.show()
                
                # 创建工作线程并连接信号
                self.worker_thread = WorkerThread(self.performMolecularDynamics, loaded_file, md_params)
                self.worker_thread.done.connect(self.onSimulationDone)
                self.worker_thread.error.connect(self.onSimulationError)
                # 启动工作线程
                self.worker_thread.start()
        elif prediction_type == 'Structure Optimization':
            # 将来实现结构优化功能
            pass
        else:
            QMessageBox.information(self, 'Invalid Operation', 'Please select a valid operation to perform.')

    ### 分子动力学模拟
    def performMolecularDynamics(self, file_path, parameters):
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
    def onSimulationDone(self):
        # 计算完成，关闭进度条对话框
        self.progress_dialog.setValue(100)
        self.progress_dialog.close()
        QMessageBox.information(self, 'Simulation Completed', 'Molecular Dynamics simulation completed successfully!')

    def onSimulationError(self, e):
        # 发生错误，关闭进度条并显示错误信息
        self.progress_dialog.close()
        QMessageBox.critical(self, 'Simulation Error', 'An error occurred during the simulation:\n' + str(e))

    ## 直接推理
    def performDirectInference(self, file_path):
        try:
            # 加载预训练的CHGNet模型
            chgnet = CHGNet.load()
            # 从CIF文件创建结构
            structure = Structure.from_file(file_path)
            # 进行预测
            prediction = chgnet.predict_structure(structure)
            
            # 显示结果
            result = ''
            for key, unit in [
                ("energy", "eV/atom"),
                ("forces", "eV/A"),
                ("stress", "GPa"),
                ("magmom", "mu_B"),
            ]:
                result += f"CHGNet-predicted {key} ({unit}):\n{prediction[key]}\n\n"
            self.result_text.setText(result)
        except Exception as e:
            QMessageBox.critical(self, 'Prediction Error', 'An error occurred during the prediction:\n' + str(e) + '\n\n' + traceback.format_exc())
            self.result_text.setText('')

window = CHGnetApp()
window.show()
sys.exit(app.exec_())