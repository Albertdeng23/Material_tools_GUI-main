import paramiko
import os

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Remote VASP Calculation")

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        self.layout = QVBoxLayout(self.central_widget)

        self.button = QPushButton("Run VASP Calculation")
        self.button.clicked.connect(self.run_vasp_calculation)

        self.layout.addWidget(self.button)

    def run_vasp_calculation(self):
        try:
            # Connect to the remote Linux server
            ssh = paramiko.SSHClient()
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh.connect('10.14.135.64', username='customer', password='234324')

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

if __name__ == "__main__":
    app = QApplication([])
    window = MainWindow()
    window.show()
    app.exec_()