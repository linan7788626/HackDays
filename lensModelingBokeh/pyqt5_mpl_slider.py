import sys

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar

import lensing_libs as ll
import numpy as np
import pylab as pl

class AppForm(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('Interactive Lens Modeling')

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()

        self.textbox.setText('1 2 3 4')
        self.on_draw()

    def save_plot(self):
        file_choices = "PNG (*.png)|*.png"

        path = QFileDialog.getSaveFileName(self,
                        'Save file', '',
                        file_choices)
        if path:
            self.canvas.print_figure(path, dpi=self.dpi)
            self.statusBar().showMessage('Saved to %s' % path, 2000)

    def on_about(self):
        msg = """ A demo of using PyQt with matplotlib:

         * Use the matplotlib navigation bar
         * Add values to the text box and press Enter (or click "Draw")
         * Show or hide the grid
         * Drag the slider to modify the width of the bars
         * Save the plot to a file using the File menu
         * Click on a bar to receive an informative message
        """
        QMessageBox.about(self, "About the demo", msg.strip())

    def on_pick(self, event):
        # The event received here is of the type
        # matplotlib.backend_bases.PickEvent
        #
        # It carries lots of information, of which we're using
        # only a small amount here.
        #
        box_points = event.artist.get_bbox().get_points()
        msg = "You've clicked on a bar with coords:\n %s" % box_points

        QMessageBox.information(self, "Click!", msg)

    def on_draw(self):
        """ Redraws the figure
        """

        bsz = 8.0 # (arcsec)
        nnn = 128
        dsx = bsz/nnn

        xi1 = np.linspace(-bsz/2.0,bsz/2.0,nnn)+0.5*dsx
        xi2 = np.linspace(-bsz/2.0,bsz/2.0,nnn)+0.5*dsx
        xi1,xi2 = np.meshgrid(xi1,xi2)
        #----------------------------------------------------------------------
        l_xcen = self.slider_x1.value()/100.0*bsz-bsz/2.0   	# x position of center (also try (0.0,0.14)
        l_ycen = self.slider_x2.value()/100.0*bsz-bsz/2.0  	# y position of center
        l_re = self.slider_re.value()/100.0*4.0   # Einstein radius of lens.
        l_rc = 0.0   # Core size of lens (in units of Einstein radius).
        l_axrat = self.slider_ql.value()/100.0   # Axis ratio of lens.
        l_pa = self.slider_la.value()/100.0*360.0  # Orintation of lens.

        lpar = np.asarray([l_xcen,l_ycen,l_re,l_rc,l_axrat,l_pa])
        #----------------------------------------------------------------------
        g_amp = 1.0   	# peak brightness value
        g_sig = self.slider_rh.value()/200.0  	# Gaussian "sigma" (i.e., size)
        g_xcen =  self.slider_y1.value()/100.0*bsz-bsz/2.0  	# x position of center (also try (0.0,0.14)
        g_ycen =  self.slider_y2.value()/100.0*bsz-bsz/2.0  	# y position of center
        g_axrat =  self.slider_qs.value()/100.0 	# minor-to-major axis ratio
        g_pa =  self.slider_sa.value()/100.0*360.0      	# major-axis position angle (degrees) c.c.w. from x axis

        gpar = np.asarray([g_amp,g_sig,g_xcen,g_ycen,g_axrat,g_pa])
        #----------------------------------------------------------------------
        g_limage,g_source,mua,yi1,yi2 = ll.lensed_images(xi1,xi2,lpar,gpar)

        #--------------------------lens images contour------------------------
        levels = [0.5,]

        self.axesb.clear()
        self.axesb.set_xlim(-bsz/2.0,bsz/2.0)
        self.axesb.set_ylim(-bsz/2.0,bsz/2.0)
        self.axesb.contour(xi1,xi2,g_limage,levels,colors=('b'))
        self.axesb.contour(xi1,xi2,mua,colors=('r'),linewidths = 2.0)

        #str = unicode(self.textbox.text())
        #self.data = list(map(int, self.textbox.text().split()))

        #x = range(len(self.data))

        # clear the axes and redraw the plot anew
        #
        #self.axesa.grid(self.grid_cb.isChecked())

        self.axesa.clear()
        self.axesa.set_xlim(-bsz/2.0,bsz/2.0)
        self.axesa.set_ylim(-bsz/2.0,bsz/2.0)
        self.axesa.contour(xi1,xi2,g_source,levels,colors=('b'))
        self.axesa.contour(yi1,yi2,mua,0,colors=('g'),linewidths = 2.0)

        #self.axes.bar(
            #left=x,
            #height=[tmp*self.slider_y1.value()/100.0 for tmp in self.data],#+list([self.slider_y1.value() / 100.0])*len(self.data),
            #width=self.slider_re.value() / 100.0,
            #align='center',
            #alpha=0.44,
            #picker=5)

        self.canvas.draw()

    def create_main_frame(self):
        self.main_frame = QWidget()
        self.main_frame.resize(800, 400)

        # Create the mpl Figure and FigCanvas objects.
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 100
        #self.fig = Figure((5.0, 5.0), dpi=self.dpi)
        self.fig = pl.figure(num=None,figsize=(10,5),dpi=80, facecolor='w', edgecolor='k')
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)

        # Since we have only one plot, we can use add_axes
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #

        #self.axes = self.fig.add_subplot(111)
        self.axesa = pl.axes([0.05,0.1,0.4,0.8])
        self.axesb = pl.axes([0.55,0.1,0.4,0.8])

        # Bind the 'pick' event for clicking on one of the bars
        #
        self.canvas.mpl_connect('pick_event', self.on_pick)

        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)

        # Other GUI controls
        #
        self.textbox = QLineEdit()
        self.textbox.setMinimumWidth(200)
        self.textbox.editingFinished.connect(self.on_draw)

        self.draw_button = QPushButton("&Optmize")
        self.draw_button.clicked.connect(self.on_draw)

        self.grid_cb = QCheckBox("Show &Grid")
        self.grid_cb.setChecked(False)
        self.grid_cb.stateChanged.connect(self.on_draw) #int

        slider_label_re = QLabel(r'$R_e$')
        self.slider_re = QSlider(Qt.Horizontal)
        self.slider_re.setRange(1, 100)
        self.slider_re.setValue(40)
        self.slider_re.setTracking(True)
        self.slider_re.setTickPosition(QSlider.TicksBothSides)
        self.slider_re.valueChanged.connect(self.on_draw)#int

        slider_label_x1 = QLabel(r'$x_1$')
        self.slider_x1 = QSlider(Qt.Horizontal)
        self.slider_x1.setRange(1, 100)
        self.slider_x1.setValue(50)
        self.slider_x1.setTracking(True)
        self.slider_x1.setTickPosition(QSlider.TicksBothSides)
        self.slider_x1.valueChanged.connect(self.on_draw)#int

        slider_label_x2 = QLabel(r'$x_2$')
        self.slider_x2 = QSlider(Qt.Horizontal)
        self.slider_x2.setRange(1, 100)
        self.slider_x2.setValue(50)
        self.slider_x2.setTracking(True)
        self.slider_x2.setTickPosition(QSlider.TicksBothSides)
        self.slider_x2.valueChanged.connect(self.on_draw)#int

        slider_label_ql = QLabel(r'$q_l$')
        self.slider_ql = QSlider(Qt.Horizontal)
        self.slider_ql.setRange(1, 100)
        self.slider_ql.setValue(50)
        self.slider_ql.setTracking(True)
        self.slider_ql.setTickPosition(QSlider.TicksBothSides)
        self.slider_ql.valueChanged.connect(self.on_draw)#int

        slider_label_la = QLabel(r'$\theta_l$')
        self.slider_la = QSlider(Qt.Horizontal)
        self.slider_la.setRange(1, 100)
        self.slider_la.setValue(30)
        self.slider_la.setTracking(True)
        self.slider_la.setTickPosition(QSlider.TicksBothSides)
        self.slider_la.valueChanged.connect(self.on_draw)#int

        slider_label_rh = QLabel(r'$R_{\rm half}$')
        self.slider_rh = QSlider(Qt.Horizontal)
        self.slider_rh.setRange(1, 100)
        self.slider_rh.setValue(50)
        self.slider_rh.setTracking(True)
        self.slider_rh.setTickPosition(QSlider.TicksBothSides)
        self.slider_rh.valueChanged.connect(self.on_draw)#int

        slider_label_y1 = QLabel(r'$y_1$')
        self.slider_y1 = QSlider(Qt.Horizontal)
        self.slider_y1.setRange(1, 100)
        self.slider_y1.setValue(50)
        self.slider_y1.setTracking(True)
        self.slider_y1.setTickPosition(QSlider.TicksBothSides)
        self.slider_y1.valueChanged.connect(self.on_draw)#int

        slider_label_y2 = QLabel(r'$y_2$')
        self.slider_y2 = QSlider(Qt.Horizontal)
        self.slider_y2.setRange(1, 100)
        self.slider_y2.setValue(50)
        self.slider_y2.setTracking(True)
        self.slider_y2.setTickPosition(QSlider.TicksBothSides)
        self.slider_y2.valueChanged.connect(self.on_draw)#int

        slider_label_qs = QLabel(r'$q_s$')
        self.slider_qs = QSlider(Qt.Horizontal)
        self.slider_qs.setRange(1, 100)
        self.slider_qs.setValue(50)
        self.slider_qs.setTracking(True)
        self.slider_qs.setTickPosition(QSlider.TicksBothSides)
        self.slider_qs.valueChanged.connect(self.on_draw)#int

        slider_label_sa = QLabel(r'\theta_s')
        self.slider_sa = QSlider(Qt.Horizontal)
        self.slider_sa.setRange(1, 100)
        self.slider_sa.setValue(70)
        self.slider_sa.setTracking(True)
        self.slider_sa.setTickPosition(QSlider.TicksBothSides)
        self.slider_sa.valueChanged.connect(self.on_draw)#int
        #
        # Layout with box sizers
        #
        hbox = QHBoxLayout()

        for w in [self.textbox, self.draw_button, self.grid_cb]:
            hbox.addWidget(w)
            hbox.setAlignment(w, Qt.AlignVCenter)

        vbox = QVBoxLayout()
        vbox.addWidget(self.mpl_toolbar)
        vbox.addWidget(self.canvas)

        hbox_s1 = QHBoxLayout()
        for w in [slider_label_re, self.slider_re, slider_label_rh, self.slider_rh]:
            hbox_s1.addWidget(w)
            hbox_s1.setAlignment(w, Qt.AlignVCenter)

        hbox_s2 = QHBoxLayout()
        for w in [slider_label_x1, self.slider_x1, slider_label_y1, self.slider_y1]:
            hbox_s2.addWidget(w)
            hbox_s2.setAlignment(w, Qt.AlignVCenter)

        hbox_s3 = QHBoxLayout()
        for w in [slider_label_x2, self.slider_x2, slider_label_y2, self.slider_y2]:
            hbox_s3.addWidget(w)
            hbox_s3.setAlignment(w, Qt.AlignVCenter)

        hbox_s4 = QHBoxLayout()
        for w in [slider_label_ql, self.slider_ql, slider_label_qs, self.slider_qs]:
            hbox_s4.addWidget(w)
            hbox_s4.setAlignment(w, Qt.AlignVCenter)

        hbox_s5 = QHBoxLayout()
        for w in [slider_label_la, self.slider_la, slider_label_sa, self.slider_sa]:
            hbox_s5.addWidget(w)
            hbox_s5.setAlignment(w, Qt.AlignVCenter)

        vbox.addLayout(hbox)
        vbox.addLayout(hbox_s1)
        vbox.addLayout(hbox_s2)
        vbox.addLayout(hbox_s3)
        vbox.addLayout(hbox_s4)
        vbox.addLayout(hbox_s5)

        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)

    def create_status_bar(self):
        self.status_text = QLabel("This is a statusbar")
        self.statusBar().addWidget(self.status_text, 1)

    def create_menu(self):
        self.file_menu = self.menuBar().addMenu("&File")

        load_file_action = self.create_action("&Save plot",
            shortcut="Ctrl+S", slot=self.save_plot,
            tip="Save the plot")
        quit_action = self.create_action("&Quit", slot=self.close,
            shortcut="Ctrl+Q", tip="Close the application")

        self.add_actions(self.file_menu,
            (load_file_action, None, quit_action))

        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About",
            shortcut='F1', slot=self.on_about,
            tip='About the demo')

        self.add_actions(self.help_menu, (about_action,))

    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(  self, text, slot=None, shortcut=None,
                        icon=None, tip=None, checkable=False,
                        signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            action.triggered.connect(slot)
        if checkable:
            action.setCheckable(True)
        return action


def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.resize(896, 768)
    form.setWindowTitle('Hoopla')
    form.show()
    app.exec_()

if __name__ == "__main__":
    main()
