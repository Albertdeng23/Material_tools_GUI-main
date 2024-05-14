import sys
from PyQt5.QtWidgets import *
from PyQt5.QtWidgets import QMessageBox
from chgnet.model.model import CHGNet
from pymatgen.core import Structure
import traceback
from chgnet.model.model import CHGNet
from chgnet.model.dynamics import MolecularDynamics
from chgnet.model import StructOptimizer
from pymatgen.core import Structure
from PyQt5.QtCore import QThread, pyqtSignal , Qt
from ase.io import read
from ase.visualize import view


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
class MolecularDynamicsThread(QThread):
    finished = pyqtSignal(str)  # Signal to emit output string

    def __init__(self, structure_file, md_params):
        QThread.__init__(self)
        self.structure_file = structure_file
        self.md_params = md_params

    def run(self):
        try:
            # Load structure from CIF file
            structure = Structure.from_file(self.structure_file)

            # Load pretrained CHGNet model
            chgnet = CHGNet.load()

            # Determine which device to use: 'cpu' or 'cuda'
            if self.md_params['device'].lower() == 'cpu':
                use_device = 'cpu'
            else:
                use_device = 'cuda'

            # Initialize MolecularDynamics with parameters, including the device
            md = MolecularDynamics(
                atoms=structure,
                model=chgnet,
                ensemble=self.md_params['ensemble'],
                temperature=self.md_params['temperature'],
                timestep=self.md_params['timestep'],
                trajectory=self.md_params['trajectory'],
                logfile=self.md_params['logfile'],
                loginterval=self.md_params['loginterval'],
                use_device=use_device  # Adding the device choice here
            )
            
            # Run Molecular Dynamics simulation with specified steps
            md.run(self.md_params['steps'])

            self.finished.emit("Molecular Dynamics simulation completed successfully.\n"
                               f"Output trajectory saved to: {self.md_params['trajectory']}\n"
                               f"Log saved to: {self.md_params['logfile']}")

        except Exception as e:
            # Send the error to the main thread
            self.finished.emit(f"An error occurred during simulation: {e}")

class StructureOptimizationThread(QThread):
    finished = pyqtSignal(str)  # Signal to emit output string
    # Signal to emit output with the relaxed structure
    relaxed_structure_signal = pyqtSignal(object)  

    def __init__(self, structure_file):
        QThread.__init__(self)
        self.structure_file = structure_file

    def run(self):
        try:
            # Load structure from CIF file
            structure = Structure.from_file(self.structure_file)
            
            # Initialize and run the StructOptimizer
            relaxer = StructOptimizer()
            result = relaxer.relax(structure)
            
            # Emit the relaxed structure through the signal
            self.relaxed_structure_signal.emit(result["final_structure"])
            
            # Prepare a message with the relaxation results
            message = "Structure optimization completed successfully.\n"
            message += f"Relaxed total energy in eV: {result['trajectory'].energies[-1]}"
            self.finished.emit(message)
        
        except Exception as e:
            # Send the error to the main thread
            self.finished.emit(f"An error occurred during structure optimization: {e}")

class CHGnetApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = "CHGnet GUI"
        self.initUI()
        self.output_text = ''
        
    def initUI(self):
        self.setWindowTitle(self.title)

    # 使用更大的字体和美观的界面
        self.setStyleSheet("""
            QMainWindow {
                font-size: 16px;
            }
            QPushButton {
                background-color: #5AA469;
                color: white;
                border-style: outset;
                border-width: 2px;
                border-radius: 10px;
                border-color: beige;
                font: bold 14px;
                min-width: 10em;
                padding: 6px;
            }
            QPushButton:hover {
                background-color: #66c280;
            }
            QLabel {
                font: 14px;
            }
            QTextEdit {
                border: 1px solid #ccc;
                padding: 6px;
            }
            QComboBox {
                min-height: 2em;
                border-radius: 5px;
            }
        """)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)
        
        # Adjust spacing and alignment
        layout.setSpacing(10) 
        layout.setAlignment(Qt.AlignTop)
        
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

        # 对于layout里的组件设置外边距
        layout.setContentsMargins(20, 20, 20, 20)

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
            ### 实现DirectInference
            self.performDirectInference(loaded_file)

        elif prediction_type == 'Molecular Dynamics':
            ### 实现Molecular Dynamics
            md_params_dialog = WindowSetMolecularDynamics(self)
            if md_params_dialog.exec_():
                md_params = md_params_dialog.getParameters()
                self.performMolecularDynamics(loaded_file, md_params)

        elif prediction_type == 'Structure Optimization':
           ## 实现Structure Optimization
           self.performStructureOptimization(loaded_file)
        else:
            QMessageBox.information(self, 'Invalid Operation', 'Please select a valid operation to perform.')
            
    def performDirectInference(self, loaded_file):
        try:
            # 加载CHGnet模型
            chgnet = CHGNet.load()
            
            # 从文件加载结构
            structure = Structure.from_file(loaded_file)
            
            # 使用模型对结构进行预测
            prediction = chgnet.predict_structure(structure)
            
            # 准备结果输出
            result_output = ""
            for key, unit in [
                ("energy", "eV/atom"),
                ("forces", "eV/A"),
                ("stress", "GPa"),
                ("magmom", "mu_B"),
            ]:
                result_output += f"CHGNet-predicted {key} ({unit}):\n{prediction[key[0]]}\n\n"
            
            # 显示结果
            self.result_text.setText(result_output)
        except Exception as e:
            # 如果出现任何错误，将其显示为警告框，并在QTextEdit中显示
            QMessageBox.warning(self, 'Error', f'Failed to perform direct inference: {e}')
            self.result_text.setText(f'Error during direct inference: {e}')
    
    ### MolecularDynamics
    def performMolecularDynamics(self, loaded_file, md_params):
        self.result_text.append("Starting Molecular Dynamics simulation...")

        # Creates a worker thread to perform computation
        self.thread = MolecularDynamicsThread(loaded_file, md_params)
        
        # Connecting signals
        self.thread.finished.connect(self.onMolecularDynamicsFinished)
        
        # Starting thread
        self.thread.start()
    def onMolecularDynamicsFinished(self, result):
        # Once calculation is done, update the UI with results
        self.result_text.append(result)

    ### Structure Optimization
    def performStructureOptimization(self, loaded_file):
        try:
            # 更新界面以显示正在进行的结构优化
            self.result_text.setText("Performing Structure Optimization...\n")
            QApplication.processEvents()

            # 从加载的CIF文件创建Structure对象
            structure = Structure.from_file(loaded_file)

            # 创建StructOptimizer实例并执行结构优化
            relaxer = StructOptimizer()
            result = relaxer.relax(structure)

            # 从结果中取得最终结构和能量信息
            final_structure = result["final_structure"]
            relaxed_total_energy = result['trajectory'].energies[-1]

            # 将结果更新到界面，显示最终结构和优化后的总能量
            result_text = "CHGNet relaxed structure:\n{}\n\n".format(final_structure)
            result_text += "relaxed total energy in eV: {}\n".format(relaxed_total_energy)
            self.result_text.setText(result_text)

        except Exception as e:
            # 出错时向用户信息面板输出错误信息
            self.result_text.setText("An error occurred during structure optimization:\n" + str(e))
            traceback.print_exc()
            
    def performStructureOptimization(self, loaded_file):
        self.result_text.append("Starting Structure Optimization...")
        
        # Creates a worker thread to perform structure optimization
        self.optimization_thread = StructureOptimizationThread(loaded_file)
        
        # Connecting signals
        self.optimization_thread.finished.connect(self.onStructureOptimizationFinished)
        self.optimization_thread.relaxed_structure_signal.connect(self.displayRelaxedStructure)
        
        # Starting thread
        self.optimization_thread.start()

    def onStructureOptimizationFinished(self, result):
        # Once calculation is done, update the UI with results
        self.result_text.append(result)

    def displayRelaxedStructure(self, relaxed_structure):
        # This method could be used to display the relaxed structure in any form you prefer.
        # Here is an example of adding the relaxed structure to a QTextEdit.
        self.result_text.append("CHGNet relaxed structure:\n" + str(relaxed_structure))
        
if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = CHGnetApp()
    window.show()
    sys.exit(app.exec_())