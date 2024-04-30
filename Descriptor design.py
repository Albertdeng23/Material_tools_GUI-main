import sys
import os
from pymatgen.io.cif import CifParser
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QPushButton, QFileDialog, QTextEdit, QMessageBox, QLabel, QProgressDialog
from collections import Counter
from mendeleev import element
import csv
from PyQt5.QtCore import Qt


class FileSelectionWidget(QWidget):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        layout = QVBoxLayout()

        self.setWindowTitle("请选择CIF文件夹")
        self.setGeometry(100, 100, 600, 200)

        self.folder_path_label = QLabel("CIF文件夹路径:")
        layout.addWidget(self.folder_path_label)

        self.folder_path_text = QTextEdit()
        self.folder_path_text.setFixedHeight(30)
        layout.addWidget(self.folder_path_text)

        self.folder_button = QPushButton("浏览文件夹")
        self.folder_button.clicked.connect(self.folderDialog)
        layout.addWidget(self.folder_button)


        self.save_path_label = QLabel("保存路径:")
        layout.addWidget(self.save_path_label)

        self.save_path_text = QTextEdit()
        self.save_path_text.setFixedHeight(30)
        layout.addWidget(self.save_path_text)

        self.save_path_button = QPushButton("选择保存路径")
        self.save_path_button.clicked.connect(self.selectSavePath)
        layout.addWidget(self.save_path_button)

        self.confirm_button = QPushButton("确认")
        self.confirm_button.clicked.connect(self.executeCode)
        layout.addWidget(self.confirm_button)

        self.setLayout(layout)

    def folderDialog(self):
        folderPath = QFileDialog.getExistingDirectory(self, '选择文件夹')
        if folderPath:
            self.folder_path = folderPath
            self.folder_path_text.setText(self.folder_path)

    def selectSavePath(self):
        savePath, _ = QFileDialog.getSaveFileName(self, '选择保存路径', '', 'CSV Files (*.csv)')
        if savePath:
            self.save_path = savePath
            self.save_path_text.setText(self.save_path)

    def executeCode(self):
        if hasattr(self, 'folder_path') and hasattr(self, 'save_path'):
            print("正在处理，请稍候.....")

            descriptor_list = ["abundance_crust", 'abundance_sea', 'atomic_number', 'atomic_radius', 'atomic_radius_rahm', 'atomic_volume',
                                           'atomic_weight', 'atomic_weight_uncertainty', 'block','c6', 'c6_gb', 'cas', 'covalent_radius_bragg',
                                           'covalent_radius_cordero','element1_covalent_radius_pyykko','covalent_radius_pyykko_double', 'covalent_radius_pyykko_triple',
                                           'cpk_color','density','description', 'dipole_polarizability','dipole_polarizability_unc','discoverers',
                                           'discovery_location', 'discovery_year','ec','econf', 'electron_affinity', 'en_allen', 'en_ghosh','en_pauling',
                                           'evaporation_heat', 'fusion_heat','gas_basicity','geochemical_class', 'glawe_number', 'goldschmidt_class','group',
                                           'group_id', 'heat_of_formation','is_monoisotopic', 'is_radioactive', 'isotopes','jmol_color','lattice_constant',
                                           'lattice_structure', 'mendeleev_number','metallic_radius', 'metallic_radius_c12','molar_heat_capacity', 'molcas_gv_color',
                                           'name','name_origin','period', 'pettifor_number', 'phase_transitions','proton_affinity','screening_constants', 'sources',
                                           'specific_heat_capacity', 'symbol', 'thermal_conductivity','uses','vdw_radius', 'vdw_radius_alvarez', 'vdw_radius_batsanov',
                                           'vdw_radius_bondi', 'vdw_radius_dreiding','vdw_radius_mm3', 'vdw_radius_rt', 'vdw_radius_truhlar','vdw_radius_uff']

            fieldnames = [f"element{number + 1}_{descriptor}" for number in range(count_element) for descriptor in descriptor_list]
            element_list = []

            progress_dialog = QProgressDialog("正在处理，请稍候.....", "取消", 0, 0, self)
            progress_dialog.setWindowTitle("处理中")
            progress_dialog.setWindowModality(Qt.WindowModal)
            progress_dialog.setMinimumDuration(0)

            file_count = len([filename for filename in os.listdir(self.folder_path) if filename.endswith(".cif")])
            processed_count = 0

            progress_dialog.setMaximum(file_count)

            with open(self.save_path, mode='a+', newline='', encoding='utf-8') as file:
                writer = csv.DictWriter(file, fieldnames=fieldnames)
                writer.writeheader()

                for filename in os.listdir(self.folder_path):
                    if filename.endswith(".cif"):
                        cif_parser = CifParser(f'{self.folder_path}/' + filename)
                        structure = cif_parser.get_structures()[0]
                        species = structure.species
                        lattice = structure.lattice
                        symbol_list = [element_.symbol for element_ in species]
                        count_element = len(symbol_list)
                        writer_dict = {}

                        for number in range(0, count_element-1):
                            writer_dict[f"element{count_element + 1}_x"] = structure.cart_coords[number][0]
                            writer_dict[f"element{count_element + 1}_y"] = structure.cart_coords[number][1]
                            writer_dict[f"element{count_element + 1}_z"] = structure.cart_coords[number][2]

                        

                        for number in range(0,count_element):
                            for each_descriptor in descriptor_list:
                                writer.writerow({f"element{number + 1}_{each_descriptor}":element(symbol_list[number]).__dict__[each_descriptor]})

                        processed_count += 1
                        progress_dialog.setValue(processed_count)

                        if progress_dialog.wasCanceled():
                            break

                progress_dialog.close()

            QMessageBox.information(self, "提示", "处理完成")
        else:
            QMessageBox.warning(self, "警告", "请先选择文件夹和保存路径")

if __name__ == '__main__':
    app = QApplication(sys.argv)

    window = FileSelectionWidget()
    window.show()

    sys.exit(app.exec_())