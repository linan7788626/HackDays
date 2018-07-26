import sys
import os
import utils


from PyQt4 import QtGui,QtCore

class SlideShowPics(QtGui.QMainWindow):

    """ SlideShowPics class defines the methods for UI and
        working logic
    """
    def __init__(self, imgLst, parent=None):
        super(SlideShowPics, self).__init__(parent)
        # self._path = path
        self._imageCache = []
        self._imagesInList = imgLst
        self._pause = False
        self._count = 0
        self.animFlag = True
        self.updateTimer = QtCore.QTimer()
        self.connect(self.updateTimer, QtCore.SIGNAL("timeout()"), self.nextImage)
        self.prepairWindow()
        self.nextImage()

    def prepairWindow(self):
        # Centre UI
        screen = QtGui.QDesktopWidget().screenGeometry(self)
        size = self.geometry()
        self.move((screen.width()-size.width())/2, (screen.height()-size.height())/2)
        self.setStyleSheet("QWidget{background-color: #000000;}")
        self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.buildUi()
        #self.showFullScreen()
        self.playPause()

    def buildUi(self):
        self.label = QtGui.QLabel()
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.setCentralWidget(self.label)

    def nextImage(self):
        """ switch to next image or previous image
        """
        if self._imagesInList:
            if self._count == len(self._imagesInList):
                self._count = 0

            self.showImageByPath(
                    self._imagesInList[self._count])

            if self.animFlag:
                self._count += 1
            else:
                self._count -= 1


    def showImageByPath(self, path):
        if path:
            image = QtGui.QImage(path)
            pp = QtGui.QPixmap.fromImage(image)
            self.label.setPixmap(pp.scaled(
                    self.label.size(),
                    QtCore.Qt.KeepAspectRatio,
                    QtCore.Qt.SmoothTransformation))

    def playPause(self):
        if not self._pause:
            self._pause = True
            self.updateTimer.start(25)
            return self._pause
        else:
            self._pause = False
            self.updateTimer.stop()

    def keyPressEvent(self, keyevent):
        """ Capture key to exit, next image, previous image,
            on Escape , Key Right and key left respectively.
        """
        event = keyevent.key()
        if event == QtCore.Qt.Key_Escape:
            self.close()
        if event == QtCore.Qt.Key_Left:
            self.animFlag = False
            self.nextImage()
        if event == QtCore.Qt.Key_Right:
            self.animFlag = True
            self.nextImage()
        if event == 32:
            self._pause = self.playPause()

def main(paths):
    if isinstance(paths, list):
        imgLst = utils.imageFilePaths(paths)
    elif isinstance(paths, str):
        imgLst =  utils.imageFilePaths([paths])
    else:
        print " You can either enter a list of paths or single path"
    app = QtGui.QApplication(sys.argv)

    if imgLst:
        window =  SlideShowPics(imgLst)
        window.show()
        window.raise_()
        app.exec_()
    else:
        msgBox = QtGui.QMessageBox()
        msgBox.setText("No Image found in any of the paths below\n\n%s" % paths)
        msgBox.setStandardButtons(msgBox.Cancel | msgBox.Open);
        if msgBox.exec_() == msgBox.Open:
            main(str(QtGui.QFileDialog.getExistingDirectory(None,
                "Select Directory to SlideShow",
                os.getcwd())))

class SlidersGroup(QtGui.QGroupBox):

    valueChanged = QtCore.pyqtSignal(int)

    def __init__(self, orientation, title, parent=None):
        super(SlidersGroup, self).__init__(title, parent)

        self.slider = QtGui.QSlider(orientation)
        self.slider.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.slider.setTickPosition(QtGui.QSlider.TicksBothSides)
        self.slider.setTickInterval(10)
        self.slider.setSingleStep(1)

        self.scrollBar = QtGui.QScrollBar(orientation)
        self.scrollBar.setFocusPolicy(QtCore.Qt.StrongFocus)

        self.dial = QtGui.QDial()
        self.dial.setFocusPolicy(QtCore.Qt.StrongFocus)

        self.slider.valueChanged.connect(self.scrollBar.setValue)
        self.scrollBar.valueChanged.connect(self.dial.setValue)
        self.dial.valueChanged.connect(self.slider.setValue)
        self.dial.valueChanged.connect(self.valueChanged)

        if orientation == QtCore.Qt.Horizontal:
            direction = QtGui.QBoxLayout.TopToBottom
        else:
            direction = QtGui.QBoxLayout.LeftToRight

        slidersLayout = QtGui.QBoxLayout(direction)
        slidersLayout.addWidget(self.slider)
        slidersLayout.addWidget(self.scrollBar)
        slidersLayout.addWidget(self.dial)
        self.setLayout(slidersLayout)

    def setValue(self, value):
        self.slider.setValue(value)

    def setMinimum(self, value):
        self.slider.setMinimum(value)
        self.scrollBar.setMinimum(value)
        self.dial.setMinimum(value)

    def setMaximum(self, value):
        self.slider.setMaximum(value)
        self.scrollBar.setMaximum(value)
        self.dial.setMaximum(value)

    def invertAppearance(self, invert):
        self.slider.setInvertedAppearance(invert)
        self.scrollBar.setInvertedAppearance(invert)
        self.dial.setInvertedAppearance(invert)

    def invertKeyBindings(self, invert):
        self.slider.setInvertedControls(invert)
        self.scrollBar.setInvertedControls(invert)
        self.dial.setInvertedControls(invert)


class Window(QtGui.QWidget):
    def __init__(self):
        super(Window, self).__init__()

        self.horizontalSliders = SlidersGroup(QtCore.Qt.Horizontal,
                "Horizontal")
        self.verticalSliders = SlidersGroup(QtCore.Qt.Vertical, "Vertical")

        self.stackedWidget = QtGui.QStackedWidget()
        self.stackedWidget.addWidget(self.horizontalSliders)
        self.stackedWidget.addWidget(self.verticalSliders)

        self.createControls("Controls")

        self.horizontalSliders.valueChanged.connect(self.verticalSliders.setValue)
        self.verticalSliders.valueChanged.connect(self.valueSpinBox.setValue)
        self.valueSpinBox.valueChanged.connect(self.horizontalSliders.setValue)

        layout = QtGui.QHBoxLayout()
        layout.addWidget(self.controlsGroup)
        layout.addWidget(self.stackedWidget)
        self.setLayout(layout)

        self.minimumSpinBox.setValue(0)
        self.maximumSpinBox.setValue(20)
        self.valueSpinBox.setValue(5)

        self.setWindowTitle("Sliders")

    def createControls(self, title):
        self.controlsGroup = QtGui.QGroupBox(title)

        minimumLabel = QtGui.QLabel("Minimum value:")
        maximumLabel = QtGui.QLabel("Maximum value:")
        valueLabel = QtGui.QLabel("Current value:")

        invertedAppearance = QtGui.QCheckBox("Inverted appearance")
        invertedKeyBindings = QtGui.QCheckBox("Inverted key bindings")

        self.minimumSpinBox = QtGui.QSpinBox()
        self.minimumSpinBox.setRange(-100, 100)
        self.minimumSpinBox.setSingleStep(1)

        self.maximumSpinBox = QtGui.QSpinBox()
        self.maximumSpinBox.setRange(-100, 100)
        self.maximumSpinBox.setSingleStep(1)

        self.valueSpinBox = QtGui.QSpinBox()
        self.valueSpinBox.setRange(-100, 100)
        self.valueSpinBox.setSingleStep(1)

        orientationCombo = QtGui.QComboBox()
        orientationCombo.addItem("Horizontal slider-like widgets")
        orientationCombo.addItem("Vertical slider-like widgets")

        orientationCombo.activated.connect(self.stackedWidget.setCurrentIndex)
        self.minimumSpinBox.valueChanged.connect(self.horizontalSliders.setMinimum)
        self.minimumSpinBox.valueChanged.connect(self.verticalSliders.setMinimum)
        self.maximumSpinBox.valueChanged.connect(self.horizontalSliders.setMaximum)
        self.maximumSpinBox.valueChanged.connect(self.verticalSliders.setMaximum)
        invertedAppearance.toggled.connect(self.horizontalSliders.invertAppearance)
        invertedAppearance.toggled.connect(self.verticalSliders.invertAppearance)
        invertedKeyBindings.toggled.connect(self.horizontalSliders.invertKeyBindings)
        invertedKeyBindings.toggled.connect(self.verticalSliders.invertKeyBindings)

        controlsLayout = QtGui.QGridLayout()
        controlsLayout.addWidget(minimumLabel, 0, 0)
        controlsLayout.addWidget(maximumLabel, 1, 0)
        controlsLayout.addWidget(valueLabel, 2, 0)
        controlsLayout.addWidget(self.minimumSpinBox, 0, 1)
        controlsLayout.addWidget(self.maximumSpinBox, 1, 1)
        controlsLayout.addWidget(self.valueSpinBox, 2, 1)
        controlsLayout.addWidget(invertedAppearance, 0, 2)
        controlsLayout.addWidget(invertedKeyBindings, 1, 2)
        controlsLayout.addWidget(orientationCombo, 3, 0, 1, 3)
        self.controlsGroup.setLayout(controlsLayout)


if __name__ == '__main__':

    curntPaths = os.getcwd()
    if len(sys.argv) > 1:
        curntPaths = sys.argv[1:]
    main(curntPaths)
