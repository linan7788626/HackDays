import os
import sys
from PyQt5.QtGui import *

# Create window
app = QApplication(sys.argv)
w = QWidget()
w.setWindowTitle("PyQT4 Pixmap @ pythonspot.com ")

# Create widget
label = QLabel(w)
pixmap = QPixmap(os.getcwd() + '/logo.png')
label.setPixmap(pixmap)
w.resize(pixmap.width(),pixmap.height())

# Draw window
w.show()
app.exec_()
