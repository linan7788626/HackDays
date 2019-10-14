from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

from matplotlib.figure import Figure
import sys
import seaborn as sns

tips = sns.load_dataset("tips")

class MainWindow(QMainWindow):
    send_fig = pyqtSignal(str)

    def __init__(self):
        super(MainWindow, self).__init__()

        self.main_widget = QWidget(self)

        self.fig = Figure()
        self.ax1 = self.fig.add_subplot(121)
        self.ax2 = self.fig.add_subplot(122, sharex=self.ax1, sharey=self.ax1)
        self.axes=[self.ax1, self.ax2]
        self.canvas = FigureCanvas(self.fig)

        self.canvas.setSizePolicy(QSizePolicy.Expanding,
                                  QSizePolicy.Expanding)
        self.canvas.updateGeometry()

        self.dropdown1 = QComboBox()
        self.dropdown1.addItems(["sex", "time", "smoker"])
        self.dropdown2 = QComboBox()
        self.dropdown2.addItems(["sex", "time", "smoker", "day"])
        self.dropdown2.setCurrentIndex(2)

        self.dropdown1.currentIndexChanged.connect(self.update)
        self.dropdown2.currentIndexChanged.connect(self.update)
        self.label = QLabel("A plot:")

        self.layout = QGridLayout(self.main_widget)
        self.layout.addWidget(QLabel("Select category for subplots"))
        self.layout.addWidget(self.dropdown1)
        self.layout.addWidget(QLabel("Select category for markers"))
        self.layout.addWidget(self.dropdown2)

        self.layout.addWidget(self.canvas)

        self.setCentralWidget(self.main_widget)
        self.show()
        self.update()

    def update(self):

        colors=["b", "r", "g", "y", "k", "c"]
        self.ax1.clear()
        self.ax2.clear()
        cat1 = self.dropdown1.currentText()
        cat2 = self.dropdown2.currentText()
        print cat1, cat2

        for i, value in enumerate(tips[cat1].unique().get_values()):
            print "value ", value
            df = tips.loc[tips[cat1] == value]
            self.axes[i].set_title(cat1 + ": " + value)
            for j, value2 in enumerate(df[cat2].unique().get_values()):
                print "value2 ", value2
                df.loc[ tips[cat2] == value2 ].plot(kind="scatter", x="total_bill", y="tip",
                                                ax=self.axes[i], c=colors[j], label=value2)
        self.axes[i].legend()
        self.fig.canvas.draw_idle()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    win = MainWindow()
    sys.exit(app.exec_())
