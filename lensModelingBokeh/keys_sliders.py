import sys
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QWidget, QSlider, QLabel, QVBoxLayout


class Widget(QWidget):
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        self.v_layout = QVBoxLayout()

        self.slider = QSlider()
        self.slider.setOrientation(Qt.Horizontal)
        self.label = QLabel('Slider at position 0')

        self.v_layout.addWidget(self.label)
        self.v_layout.addWidget(self.slider)

        self.setLayout(self.v_layout)

        self.slider.valueChanged.connect(self.slider_moved)

    def keyPressEvent(self, event):
        if event.key()==Qt.Key_Right:
            self.slider.setValue(self.slider.value() + 1)
        elif event.key()==Qt.Key_Left:
            self.slider.setValue(self.slider.value() - 1)
        else:
            QWidget.keyPressEvent(self, event)

    def slider_moved(self, position):
        self.label.setText('Slider at position %d' % position)


if __name__ == '__main__':
    app = QApplication(sys.argv)

    widget = Widget()
    widget.show()

    sys.exit(app.exec_())
