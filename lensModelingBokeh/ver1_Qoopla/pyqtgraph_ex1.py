# import python standard modules
import sys
# import 3rd party libraries
from PyQt5 import QtCore, QtGui
import numpy as np
import pyqtgraph as pg
#import local python

class Window(QtGui.QMainWindow):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.setWindowTitle('Simple Sinc Function UI')
        self.window = pg.GraphicsWindow()
        self.window.setBackground('w')

        self.contour_plot = self.window.addViewBox()
        self.setCentralWidget(self.window)
        # add a slider to change the coda of the sinc function
        self.slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        # initialize the time axis (this will not change)
        self.t = np.linspace(-0.500, 0.500, num=1000, endpoint=True)
        # a will change the coda of the sinc function
        self.a = 10
        #
        # add the menubar with the method createMenuBar()
        self.createMenuBar()
        # add the dock widget with the method createDockWidget()
        self.createDockWidget()
        #
        # first set the default value to a
        self.slider.setValue(self.a)
        # when the slider is changed, it emits a signal and sends an integer value
        # we send that value to a method called slider value changed that updates the value a
        self.slider.valueChanged.connect(self.sliderValueChanged)
        # finally draw the curve
    def createMenuBar(self):
        # file menu actions
        exit_action = QtGui.QAction('&Exit', self)
        exit_action.triggered.connect(self.close)
        # create an instance of menu bar
        menubar = self.menuBar()
        # add file menu and file menu actions
        file_menu = menubar.addMenu('&File')
        file_menu.addAction(exit_action)

    def createDockWidget(self):
        my_dock_widget = QtGui.QDockWidget()
        my_dock_widget.setObjectName('Control Panel')
        my_dock_widget.setAllowedAreas(QtCore.Qt.TopDockWidgetArea | QtCore.Qt.BottomDockWidgetArea)
        my_dock_widget.setTitleBarWidget(QtGui.QWidget(my_dock_widget))
        # create a widget to house user control widgets like sliders
        my_house_widget = QtGui.QWidget()
        # every widget should have a layout, right?
        my_house_layout = QtGui.QVBoxLayout()
        # add the slider initialized in __init__() to the layout
        my_house_layout.addWidget(self.slider)
        # apply the 'house' layout to the 'house' widget
        my_house_widget.setLayout(my_house_layout)
        # set the house widget 'inside' the dock widget
        my_dock_widget.setWidget(my_house_widget)
        # now add the dock widget to the main window
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, my_dock_widget)

    def sliderValueChanged(self, int_value):
        self.a = int_value
        self.drawCurve()

    def drawCurve(self):
        x = np.linspace(0,6.28,30)
        y = x[:]
        xx,yy = np.meshgrid(x,y)
        img = pg.ImageItem(np.sin(xx)+np.cos(yy))
        # self.contour_plot.addItem(img)
        # self.contour_plot.setAspectLocked()

        curves = []
        c = pg.IsocurveItem(level=1,pen='r')
        c.setParentItem(img)
        c.setZValue(10)
        curves.append(c)
        curves[0].setData(np.sin(xx)+np.cos(yy))

        # self.contour_plot.addItem(c)

        # self.my_sinc = np.sin(self.a * np.pi * self.t) / (self.a * np.pi * self.t)

def main():
    app = QtGui.QApplication(sys.argv)
    app.setApplicationName('Simple Sinc Function UI')
    window = Window()
    window.show()
    app.exec_()

if __name__ == '__main__':
    main()
