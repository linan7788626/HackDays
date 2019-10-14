import sys
import seaborn as sns
import matplotlib.pyplot as plt

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

tips = sns.load_dataset("tips")


def seabornplot():
    g = sns.FacetGrid(tips, col="sex", hue="time", palette="Set1",
                                hue_order=["Dinner", "Lunch"])
    g.map(plt.scatter, "total_bill", "tip", edgecolor="w")
    return g.fig


class MainWindow(QMainWindow):
    send_fig = pyqtSignal(str)

    def __init__(self):
        super(MainWindow, self).__init__()

        self.main_widget = QWidget(self)

        self.fig = seabornplot()
        self.canvas = FigureCanvas(self.fig)

        self.canvas.setSizePolicy(QSizePolicy.Expanding,
                      QSizePolicy.Expanding)
        self.canvas.updateGeometry()
        self.button = QPushButton("Button")
        self.label = QLabel("A plot:")

        self.layout = QGridLayout(self.main_widget)
        self.layout.addWidget(self.button)
        self.layout.addWidget(self.label)
        self.layout.addWidget(self.canvas)

        self.setCentralWidget(self.main_widget)
        self.show()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    win = MainWindow()
    sys.exit(app.exec_())
