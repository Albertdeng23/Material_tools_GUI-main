import sys
import traceback
from PyQt5.QtWidgets import *
from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtCore import QThread, pyqtSignal, Qt

# 科学计算相关库导入
from pymatgen.core import Structure
from chgnet.model.model import CHGNet
from chgnet.model.dynamics import MolecularDynamics
from chgnet.model import StructOptimizer

###----------------------------------------------------------------------------------------------------------####
##  辅助窗口类
###----------------------------------------------------------------------------------------------------------####

class WindowSetMolecularDynamics(QDialog):
    """
    分子动力学 (MD) 参数设置对话框。
    允许用户输入系综、温度、步长等参数。
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        
        self.setWindowTitle("设置分子动力学参数")
        
        # 使用表单布局，方便排列 "标签-输入框" 对
        self.form_layout = QFormLayout()
        
        # --- 创建输入字段并设置默认值 ---
        self.ensemble_input = QLineEdit(self)
        self.ensemble_input.setText("nvt")  # 默认系综：NVT
        
        self.temperature_input = QLineEdit(self)
        self.temperature_input.setText("1000")  # 默认温度：1000 K
        
        self.timestep_input = QLineEdit(self)
        self.timestep_input.setText("2")  # 默认时间步长：2 fs
        
        self.traj_input = QLineEdit(self)
        self.traj_input.setText("md_out.traj")  # 轨迹文件输出路径
        
        self.logfile_input = QLineEdit(self)
        self.logfile_input.setText("md_out.log")  # 日志文件输出路径
        
        self.loginterval_input = QLineEdit(self)
        self.loginterval_input.setText("100")  # 日志记录间隔步数
        
        self.steps_input = QLineEdit(self)
        self.steps_input.setText("5000")  # 总模拟步数
        
        # 设备选择 (CPU 或 CUDA)
        self.device_combo = QComboBox(self)
        self.device_combo.addItem("cpu", "cpu")
        self.device_combo.addItem("cuda", "cuda")

        # --- 将输入字段加入布局 ---
        self.form_layout.addRow("系综 (Ensemble)", self.ensemble_input)
        self.form_layout.addRow("温度 (Temperature K)", self.temperature_input)
        self.form_layout.addRow("时间步长 (Time Step fs)", self.timestep_input)
        self.form_layout.addRow("轨迹文件 (Trajectory File)", self.traj_input)
        self.form_layout.addRow("日志文件 (Log File)", self.logfile_input)
        self.form_layout.addRow("记录间隔 (Log Interval)", self.loginterval_input)
        self.form_layout.addRow("模拟步数 (Number of Steps)", self.steps_input)
        self.form_layout.addRow("计算设备 (Device)", self.device_combo)
        
        self.setLayout(self.form_layout)

        # --- 添加 确定(OK) 和 取消(Cancel) 按钮 ---
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel, 
            parent=self
        )
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self.form_layout.addWidget(self.buttons)

    def getParameters(self):
        """
        获取用户在对话框中输入的参数，并转换为相应的数据类型。
        返回:
            dict: 包含所有MD参数的字典。
        """
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

###----------------------------------------------------------------------------------------------------------####
##  后台工作线程类
###----------------------------------------------------------------------------------------------------------####

class MolecularDynamicsThread(QThread):
    """
    分子动力学模拟的工作线程。
    使用 QThread 将耗时的 MD 模拟放在后台运行，避免阻塞主界面。
    """
    finished = pyqtSignal(str)  # 信号：任务完成或出错时发送文本消息

    def __init__(self, structure_file, md_params):
        QThread.__init__(self)
        self.structure_file = structure_file
        self.md_params = md_params

    def run(self):
        """线程的主执行函数"""
        try:
            # 1. 从 CIF 文件加载结构
            structure = Structure.from_file(self.structure_file)

            # 2. 加载预训练的 CHGNet 模型
            # 注意：这可能需要下载模型或加载本地缓存，比较耗时
            chgnet = CHGNet.load()

            # 3. 确定使用的计算设备 ('cpu' 或 'cuda')
            if self.md_params['device'].lower() == 'cpu':
                use_device = 'cpu'
            else:
                use_device = 'cuda'

            # 4. 初始化分子动力学模拟器
            md = MolecularDynamics(
                atoms=structure,
                model=chgnet,
                ensemble=self.md_params['ensemble'],
                temperature=self.md_params['temperature'],
                timestep=self.md_params['timestep'],
                trajectory=self.md_params['trajectory'],
                logfile=self.md_params['logfile'],
                loginterval=self.md_params['loginterval'],
                use_device=use_device  # 传入设备参数
            )
            
            # 5. 运行指定步数的模拟
            md.run(self.md_params['steps'])

            # 6. 发送成功信号
            self.finished.emit("分子动力学模拟成功完成。\n"
                               f"轨迹文件已保存至: {self.md_params['trajectory']}\n"
                               f"日志文件已保存至: {self.md_params['logfile']}")

        except Exception as e:
            # 捕获异常并将错误信息发送回主线程
            self.finished.emit(f"模拟过程中发生错误: {e}")


class StructureOptimizationThread(QThread):
    """
    结构优化的工作线程。
    用于在后台执行晶体结构弛豫。
    """
    finished = pyqtSignal(str)  # 信号：发送文本消息（如能量信息）
    relaxed_structure_signal = pyqtSignal(object)  # 信号：发送优化后的 Structure 对象

    def __init__(self, structure_file):
        QThread.__init__(self)
        self.structure_file = structure_file

    def run(self):
        """线程的主执行函数"""
        try:
            # 1. 从 CIF 文件加载结构
            structure = Structure.from_file(self.structure_file)
            
            # 2. 初始化并运行结构优化器
            relaxer = StructOptimizer()
            result = relaxer.relax(structure)
            
            # 3. 通过信号发送优化后的结构对象
            self.relaxed_structure_signal.emit(result["final_structure"])
            
            # 4. 准备并发送包含能量信息的文本消息
            message = "结构优化成功完成。\n"
            # 获取轨迹中最后一步的能量
            message += f"优化后的总能量 (eV): {result['trajectory'].energies[-1]}"
            self.finished.emit(message)
        
        except Exception as e:
            # 捕获异常并发送错误信息
            self.finished.emit(f"结构优化过程中发生错误: {e}")

###----------------------------------------------------------------------------------------------------------####
##  主应用程序窗口类
###----------------------------------------------------------------------------------------------------------####

class CHGnetApp(QMainWindow):
    """
    CHGNet GUI 主窗口。
    提供文件加载、功能选择（直接推理、MD、结构优化）和结果显示功能。
    """
    def __init__(self):
        super().__init__()
        self.title = "CHGnet GUI 工具"
        self.initUI()
        self.output_text = ''
        
    def initUI(self):
        self.setWindowTitle(self.title)

        # 设置样式表：调整字体大小、按钮颜色和布局间距
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
        
        # 调整布局的间距和对齐方式
        layout.setSpacing(10) 
        layout.setAlignment(Qt.AlignTop)
        
        # --- 1. 文件选择区域 ---
        self.file_button = QPushButton("加载 CIF 文件")
        self.file_button.clicked.connect(self.openFileNameDialog)
        self.file_label = QLabel("请加载一个 CIF 文件以继续")
        layout.addWidget(self.file_button)
        layout.addWidget(self.file_label)

        # --- 2. 结果显示区域 ---
        self.result_text = QTextEdit(central_widget)
        self.result_text.setReadOnly(True) # 设置为只读
        layout.addWidget(self.result_text)
        
        # --- 3. 功能选择与执行区域 ---
        self.prediction_type_combo = QComboBox(central_widget)
        self.prediction_type_combo.addItem("直接推理 (Direct Inference)")
        self.prediction_type_combo.addItem("分子动力学 (Molecular Dynamics)")
        self.prediction_type_combo.addItem("结构优化 (Structure Optimization)")
        layout.addWidget(self.prediction_type_combo)

        self.run_button = QPushButton("运行 CHGnet")
        self.run_button.clicked.connect(self.runCHGnet)
        layout.addWidget(self.run_button)

        # 设置整体布局的外边距
        layout.setContentsMargins(20, 20, 20, 20)

    def openFileNameDialog(self):
        """打开文件选择对话框"""
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getOpenFileName(self, "加载 CIF 文件", "", 
                                                  "CIF Files (*.cif);;All Files (*)",
                                                  options=options)
        if fileName:
            self.file_label.setText(f"Loaded file: {fileName}")

    def runCHGnet(self):
        """
        主执行函数。
        根据下拉菜单的选择，分发任务到相应的处理函数。
        """
        # 解析文件路径
        loaded_file = self.file_label.text().replace('Loaded file: ', '').strip()
        
        # 检查用户是否已加载文件
        if 'Loaded file' not in self.file_label.text() or not loaded_file:
            QMessageBox.information(self, '加载文件', '请在运行模型前先加载 CIF 文件。')
            return

        prediction_type = self.prediction_type_combo.currentText()
        
        # 根据选择执行不同逻辑
        if 'Direct Inference' in prediction_type:
            # 执行直接推理 (静态计算)
            self.performDirectInference(loaded_file)

        elif 'Molecular Dynamics' in prediction_type:
            # 执行分子动力学
            md_params_dialog = WindowSetMolecularDynamics(self)
            if md_params_dialog.exec_(): # 如果用户点击了 OK
                md_params = md_params_dialog.getParameters()
                self.performMolecularDynamics(loaded_file, md_params)

        elif 'Structure Optimization' in prediction_type:
           # 执行结构优化
           self.performStructureOptimization(loaded_file)
        else:
            QMessageBox.information(self, '无效操作', '请选择一个有效的操作。')
            
    def performDirectInference(self, loaded_file):
        """
        执行直接推理（静态计算）。
        注意：此处在主线程运行，如果模型加载很慢可能会短暂卡顿界面。
        """
        try:
            self.result_text.setText("正在加载模型并进行推理，请稍候...")
            QApplication.processEvents() # 刷新界面显示

            # 加载 CHGnet 模型
            chgnet = CHGNet.load()
            
            # 从文件加载结构
            structure = Structure.from_file(loaded_file)
            
            # 使用模型对结构进行预测
            prediction = chgnet.predict_structure(structure)
            
            # 格式化输出结果
            result_output = ""
            for key, unit in [
                ("energy", "eV/atom"),
                ("forces", "eV/A"),
                ("stress", "GPa"),
                ("magmom", "mu_B"),
            ]:
                # 提取预测字典中的对应值
                val = prediction.get(key[0] if key == "energy" else key) # energy 在字典中通常是 'e' 或 'energy'，需根据实际库调整，这里假设库返回标准键
                # CHGNet predict_structure 返回的键通常是 'e', 'f', 's', 'm'
                # 为了通用性，这里直接打印整个 prediction 字典可能更安全，或者根据 CHGNet 版本调整
                # 下面代码假设 prediction 字典包含完整的键名
                result_output += f"CHGNet-predicted {key} ({unit}):\n{prediction.get(key[0])}\n\n"
            
            # 显示结果
            self.result_text.setText(str(prediction)) # 简单起见，直接显示完整字典，或者使用上面的格式化
        except Exception as e:
            QMessageBox.warning(self, '错误', f'直接推理失败: {e}')
            self.result_text.setText(f'推理过程中出错: {e}')
            traceback.print_exc()
    
    def performMolecularDynamics(self, loaded_file, md_params):
        """启动分子动力学模拟线程"""
        self.result_text.append("正在启动分子动力学模拟...")

        # 创建工作线程
        self.thread = MolecularDynamicsThread(loaded_file, md_params)
        
        # 连接信号：当线程完成时调用 onMolecularDynamicsFinished
        self.thread.finished.connect(self.onMolecularDynamicsFinished)
        
        # 启动线程
        self.thread.start()

    def onMolecularDynamicsFinished(self, result):
        """MD 模拟完成后的回调函数"""
        self.result_text.append(result)

    def performStructureOptimization(self, loaded_file):
        """启动结构优化线程"""
        self.result_text.append("正在启动结构优化...")
        
        # 创建结构优化工作线程
        self.optimization_thread = StructureOptimizationThread(loaded_file)
        
        # 连接信号
        self.optimization_thread.finished.connect(self.onStructureOptimizationFinished)
        self.optimization_thread.relaxed_structure_signal.connect(self.displayRelaxedStructure)
        
        # 启动线程
        self.optimization_thread.start()

    def onStructureOptimizationFinished(self, result):
        """结构优化完成后的回调函数（处理文本消息）"""
        self.result_text.append(result)

    def displayRelaxedStructure(self, relaxed_structure):
        """
        显示优化后的结构。
        参数:
            relaxed_structure: pymatgen Structure 对象
        """
        self.result_text.append("CHGNet 优化后的结构:\n" + str(relaxed_structure))
        
if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = CHGnetApp()
    window.show()
    sys.exit(app.exec_())
