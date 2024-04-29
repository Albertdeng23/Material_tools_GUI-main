from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QLineEdit, QPushButton, QMenuBar, QAction, QMessageBox
from pymatgen.ext.matproj import MPRester
from pymatgen.io.cif import CifWriter
import sys
from PyQt5.QtGui import QFont


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("生成CIF文件")
        self.resize(400, 250)

        # 创建API密钥和结构类型的标签
        api_key_label = QLabel("API密钥：", self)
        api_key_label.move(50, 50)
        api_key_label.setStyleSheet("font-size: 25px;")  # 设置字体大小
        self.api_key_entry = QLineEdit(self)
        self.api_key_entry.setGeometry(200, 50, 200, 30)  # 设置输入框的位置和大小（left, top, width, height）
        self.api_key_entry.setStyleSheet("font-size: 20px;")  # 设置字体大小

        structure_type_label = QLabel("结构类型：", self)
        structure_type_label.move(50, 100)
        structure_type_label.setStyleSheet("font-size: 25px;")
        self.structure_type_entry = QLineEdit(self)
        self.structure_type_entry.setGeometry(200, 100, 200, 30)
        self.structure_type_entry.setStyleSheet("font-size: 20px;")

        # 创建生成按钮
        generate_button = QPushButton("生成CIF文件", self)
        generate_button.setGeometry(150, 150, 200, 50)  # 设置按钮的位置和大小
        generate_button.setStyleSheet("font-size: 16px;")
        generate_button.clicked.connect(self.generate_cif_files)

    def generate_cif_files(self):
        api_key = self.api_key_entry.text()
        structure_type = self.structure_type_entry.text()

        mpr = MPRester(api_key)
        entries = mpr.get_entries(structure_type, inc_structure=True)

        for entry in entries:
            material_id = entry.entry_id
            structure = entry.structure
            cif_writer = CifWriter(structure, write_magmoms=False)
            cif_writer.write_file(f"{material_id}.cif")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())