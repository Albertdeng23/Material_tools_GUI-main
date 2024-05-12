from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QLineEdit, QPushButton, QFileDialog, QMessageBox, QProgressDialog
from pymatgen.ext.matproj import MPRester
from pymatgen.io.cif import CifWriter
from PyQt5.QtCore import Qt
import sys
from PyQt5.QtGui import QPixmap

class CIFDownloader(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("生成CIF文件")
        self.resize(900, 450)

        self.selected_folder_label = QLabel("已选择的路径：", self)
        self.selected_folder_label.setGeometry(50, 50, 300, 30)
        self.selected_folder_label.setStyleSheet("font-size: 25px;")

        self.selected_folder_path_textbox = QLineEdit(self)
        self.selected_folder_path_textbox.setGeometry(220, 50, 500, 30)
        self.selected_folder_path_textbox.setReadOnly(True)
        self.selected_folder_path_textbox.setStyleSheet("font-size: 20px;")

        browse_button = QPushButton("浏览文件目录", self)
        browse_button.setGeometry(50, 220, 200, 50)
        browse_button.setStyleSheet("font-size: 20px;")
        browse_button.clicked.connect(self.browse_folder)

        generate_button = QPushButton("生成CIF文件", self)
        generate_button.setGeometry(250, 220, 200, 50)
        generate_button.setStyleSheet("font-size: 20px;")
        generate_button.clicked.connect(self.generate_cif_files)

        api_key_label = QLabel("API密钥：", self)
        api_key_label.move(50, 115)
        api_key_label.setStyleSheet("font-size: 25px;")
        self.api_key_entry = QLineEdit(self)
        self.api_key_entry.setGeometry(180, 115, 400, 30)
        self.api_key_entry.setStyleSheet("font-size: 20px;")
        default_value = "SwXzTEMOZZf2Pr7B36xvIeJTOTc2IBjK"
        self.api_key_entry.setText(default_value)

        structure_type_label = QLabel("结构类型：", self)
        structure_type_label.move(50, 180)
        structure_type_label.setStyleSheet("font-size: 25px;")
        self.structure_type_entry = QLineEdit(self)
        self.structure_type_entry.setGeometry(180, 180, 400, 30)
        self.structure_type_entry.setStyleSheet("font-size: 20px;")

        self.selected_folder_path = ""

        self.tips_lable = QLabel("tips.可填写的结构类型:化学式（比如：'Li2O'）化学系统（比如：'oxide'） \n     材料ID（比如：'mp-1234'）物质名称（比如：'NaCl'）",self)
        self.tips_lable .setGeometry(10,350,1000,60)
        self.tips_lable.setStyleSheet("font-size:25px;")

        CDlogo = QPixmap(r"D:\conda\CDlogo.png")
        scaled_CDlogo = CDlogo.scaled(200, 200)
        self.logo_label = QLabel(self)
        self.logo_label.setPixmap(scaled_CDlogo)
        self.logo_label.resize(scaled_CDlogo .width(), scaled_CDlogo .height())
        self.logo_label.setGeometry(600,100,200,200)

    def generate_cif_files(self):
        api_key = self.api_key_entry.text()
        structure_type = self.structure_type_entry.text()

        if not self.selected_folder_path:
            QMessageBox.warning(self, "警告", "请先选择文件路径")
            return

        mpr = MPRester(api_key)
        entries = mpr.get_entries(structure_type, inc_structure=True)
        total_entries = len(entries)

        # 弹出窗口显示 "请稍后"
        QMessageBox.information(self, "请稍候", "数据正在准备中，请稍候...")

        # 使用 QProgressBar 显示下载进度
        progress_dialog = QProgressDialog("下载中...", "取消", 0, total_entries, self)
        progress_dialog.setWindowModality(Qt.WindowModal)

        for idx, entry in enumerate(entries, 1):
            material_id = entry.entry_id
            structure = entry.structure
            cif_writer = CifWriter(structure, write_magmoms=False)
            cif_writer.write_file(f"{self.selected_folder_path}/{material_id}.cif")
            progress_dialog.setValue(idx)

            if progress_dialog.wasCanceled():
                break

        if progress_dialog.wasCanceled():
            QMessageBox.warning(self, "下载失败", "未完成下载")
        else:
            QMessageBox.information(self, "完成", f"CIF文件已生成到目录：{self.selected_folder_path}")

    def browse_folder(self):
        folder_path = QFileDialog.getExistingDirectory(self, "选择文件目录", ".")
        if folder_path:
            self.selected_folder_path = folder_path
            self.selected_folder_path_textbox.setText(folder_path)
            self.selected_folder_label.setText(f"已选择的路径：{folder_path}")

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = CIFDownloader()
    window.show()
    sys.exit(app.exec_())