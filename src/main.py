from PyQt5.QtWidgets import QApplication
import sys

from ui.main_window import MoleculeAppWindow

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MoleculeAppWindow()
    window.show()
    sys.exit(app.exec_())
