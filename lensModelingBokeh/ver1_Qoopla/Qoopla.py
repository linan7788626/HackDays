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
from PIL import Image
import ofits as ff

nbins = 1000

class AppForm(QMainWindow):
    def __init__(self, parent=None):
        self.bsz = 5.0 # (arcsec)
        self.nnn = 500
        self.dsx = self.bsz/self.nnn

        self.xi1 = np.linspace(-self.bsz/2.0,self.bsz/2.0,self.nnn)+0.5*self.dsx
        self.xi2 = np.linspace(-self.bsz/2.0,self.bsz/2.0,self.nnn)+0.5*self.dsx
        self.xi1,self.xi2 = np.meshgrid(self.xi1,self.xi2)

        self.file_name = "./ds9_images.png"
        a_tmp = Image.open(self.file_name)
        a_tmp = np.array(a_tmp)
        self.lgals = np.flipud(a_tmp[:,:,0]/256.0)

        self.srcs_name = "./ds9_srcs.png"
        s_tmp = Image.open(self.srcs_name)
        s_tmp = np.array(s_tmp)
        # # sgals_tmp = (s_tmp[:,:,0]+s_tmp[:,:,1]+s_tmp[:,:,2])/256.0/2.0
        # # sgals_tmp = (s_tmp[:,:,0])/256.0
        # self.sgals = self.lgals*0.0
        # xoff = -15
        # yoff = -10
        # self.sgals[200+xoff:233+xoff,200+yoff:233+yoff] = sgals_tmp
        self.sgals = np.flipud(s_tmp[:,:,0]/256.0)

        st_re = nbins/2
        st_x1 = nbins/2
        st_x2 = nbins/2
        st_rc = 0*nbins
        st_ql = nbins/2
        st_pa = nbins/3

        # self.lpar = np.asarray([(st_x1*self.bsz/nbins-self.bsz/2.0),(st_x2*self.bsz/nbins-self.bsz/2.0),st_re*4.0/nbins,st_rc*0.0,st_ql*1.0/nbins,st_pa*360.0/nbins])
        self.lpar = np.asarray([(st_x1*self.bsz/nbins-self.bsz/2.0 + self.bsz/nbins/2.0),
                                (st_x2*self.bsz/nbins-self.bsz/2.0 + self.bsz/nbins/2.0),
                                st_re*2.0/nbins, st_rc*0.0, st_ql*1.0/nbins, st_pa*360.0/nbins])
        #----------------------------------------------------------------------
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('Interactive Lens Modeling')

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()

        self.textbox.setText('67.26 0.32 0.68')
        self.on_draw()

        #= QLabel("Re = %.2f" % self.lpar[2])
        #label_nombre.setText("hola")

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

        l_xcen = self.slider_x1.value()*self.bsz/nbins-self.bsz/2.0 + self.bsz/nbins/2.0  	# x position of center (also try (0.0,0.14)
        l_ycen = self.slider_x2.value()*self.bsz/nbins-self.bsz/2.0 + self.bsz/nbins/2.0 	# y position of center
        l_re = self.slider_re.value()*2.0/nbins   # Einstein radius of lens.
        l_rc = 0.0   # Core size of lens (in units of Einstein radius).
        l_axrat = self.slider_ql.value()*1.0/nbins   # Axis ratio of lens.
        l_pa = self.slider_la.value()*360.0/nbins  # Orintation of lens.

        self.lpar = np.asarray([l_xcen,l_ycen,l_re,l_rc,l_axrat,l_pa])
        #----------------------------------------------------------------------
        g_amp = 1.0   	# peak brightness value
        g_sig = self.slider_rh.value()*0.5/nbins  	# Gaussian "sigma" (i.e., size)
        g_xcen =  self.slider_y1.value()*self.bsz/4.0/nbins-self.bsz/8.0 + self.bsz/nbins/8.0 	# x position of center (also try (0.0,0.14)
        g_ycen =  self.slider_y2.value()*self.bsz/4.0/nbins-self.bsz/8.0 + self.bsz/nbins/8.0 	# y position of center
        g_axrat =  self.slider_qs.value()*1.0/nbins 	# minor-to-major axis ratio
        g_pa =  self.slider_sa.value()*360.0/nbins      	# major-axis position angle (degrees) c.c.w. from x axis

        self.spar = np.asarray([g_xcen,g_ycen,g_sig,g_amp,g_axrat,g_pa])

        print "plot***", self.spar
        #----------------------------------------------------------------------
        g_limage,g_source,mua,yi1,yi2 = ll.lensed_images(self.xi1,self.xi2,self.lpar,self.spar)
        #--------------------------lens images contour------------------------

        levels = [0.5,]
        levels_imgs = [0.0,0.08,0.1,0.2,0.3,0.4,0.5]

        self.axesa.clear()
        self.axesa.set_xlim(-self.bsz/8.0,self.bsz/8.0)
        self.axesa.set_ylim(-self.bsz/8.0,self.bsz/8.0)
        self.axesa.set_xticks([-0.4, -0.2, 0.0, 0.2, 0.4])
        self.axesa.set_yticks([-0.4, -0.2, 0.0, 0.2, 0.4])
        self.axesa.contourf(self.xi1,self.xi2,self.sgals,levels_imgs,cmap=pl.cm.gray)
        self.axesa.contour(self.xi1,self.xi2,g_source,levels,colors=('deepskyblue'))
        self.axesa.contour(yi1,yi2,mua,0,colors=('r'),linewidths = 2.0)

        self.axesb.clear()
        # pl.xticks([-2.0, -1.0, 0.0, 1.0, 2.0])
        # pl.yticks([-2.0, -1.0, 0.0, 1.0, 2.0])
        self.axesb.set_xlim(-self.bsz/2.0,self.bsz/2.0)
        self.axesb.set_ylim(-self.bsz/2.0,self.bsz/2.0)
        self.axesb.contourf(self.xi1,self.xi2,self.lgals,levels_imgs,cmap=pl.cm.gray)
        self.axesb.contour(self.xi1,self.xi2,g_limage,levels,colors=('deepskyblue'))
        self.axesb.contour(self.xi1,self.xi2,mua,colors=('r'),linewidths = 2.0)

        self.canvas.draw()

        self.slider_label_re.setText("Re: %.2f" % self.lpar[2])

        # print "Drag", self.lpar

    # @profile
    def create_main_frame(self):
        self.main_frame = QWidget()
        #self.main_frame.resize(800, 400)
        self.main_frame.setFixedSize(820, 700)

        self.dpi = 80
        self.fig = pl.figure(num=None,figsize=(10,5), dpi=self.dpi, facecolor='w', edgecolor='k')
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)

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

        self.cb_cosmo = QComboBox()
        self.cb_cosmo.addItems(["Planck15", "WMP9", "Concordance", "Einstein-deSitter"])
        #self.cb_cosmo.currentIndexChanged.connect(self.selectionchange)

        self.cb_models = QComboBox()
        self.cb_models.addItems(["SIE", "NFW", "PIEMD"])
        #self.cb_models.currentIndexChanged.connect(self.selectionchange)

        self.cb_codes = QComboBox()
        self.cb_codes.addItems(["Scipy", "Q-Lens", "Lensed", "LensTool", "Gravlens", "Glafic"])
        #self.cb.currentIndexChanged.connect(self.selectionchange)


        self.draw_button = QPushButton("&Optimize")
        #self.draw_button.clicked.connect(lambda: self.optmization_actions("lensed"))
        self.draw_button.clicked.connect(self.optimization_actions)

        #self.grid_cb = QCheckBox("Show &Grid")
        #self.grid_cb.setChecked(False)
        #self.grid_cb.stateChanged.connect(self.on_draw) #int


        self.slider_label_re = QLabel("Re = %.2f" % self.lpar[2])
        self.slider_re = QSlider(Qt.Horizontal)
        self.slider_re.setMinimumWidth(270)
        self.slider_re.setMaximumWidth(270)
        self.slider_re.setRange(1, nbins)
        self.slider_re.setValue(nbins/2)
        self.slider_re.setTracking(True)
        #self.slider_re.setTickPosition(QSlider.TicksBothSides)
        self.slider_re.valueChanged.connect(self.on_draw)#int

        slider_label_x1 = QLabel("Lens Position x1")
        self.slider_x1 = QSlider(Qt.Horizontal)
        self.slider_x1.setMinimumWidth(270)
        self.slider_x1.setMaximumWidth(270)
        self.slider_x1.setRange(1, nbins)
        self.slider_x1.setValue(nbins/2)
        self.slider_x1.setTracking(True)
        #self.slider_x1.setTickPosition(QSlider.TicksBothSides)
        self.slider_x1.valueChanged.connect(self.on_draw)#int

        slider_label_x2 = QLabel("Lens Position x2")
        self.slider_x2 = QSlider(Qt.Horizontal)
        self.slider_x2.setMinimumWidth(270)
        self.slider_x2.setMaximumWidth(270)
        self.slider_x2.setRange(1, nbins)
        self.slider_x2.setValue(nbins/2)
        self.slider_x2.setTracking(True)
        #self.slider_x2.setTickPosition(QSlider.TicksBothSides)
        self.slider_x2.valueChanged.connect(self.on_draw)#int

        slider_label_ql = QLabel("Ellipticity L")
        self.slider_ql = QSlider(Qt.Horizontal)
        self.slider_ql.setMinimumWidth(270)
        self.slider_ql.setMaximumWidth(270)
        self.slider_ql.setRange(1, nbins)
        self.slider_ql.setValue(nbins/2)
        self.slider_ql.setTracking(True)
        #self.slider_ql.setTickPosition(QSlider.TicksBothSides)
        self.slider_ql.valueChanged.connect(self.on_draw)#int

        slider_label_la = QLabel("Orintation L")
        self.slider_la = QSlider(Qt.Horizontal)
        self.slider_la.setMinimumWidth(270)
        self.slider_la.setMaximumWidth(270)
        self.slider_la.setRange(1, nbins)
        self.slider_la.setValue(nbins/3)
        self.slider_la.setTracking(True)
        #self.slider_la.setTickPosition(QSlider.TicksBothSides)
        self.slider_la.valueChanged.connect(self.on_draw)#int

        slider_label_rh = QLabel("Half-light Radius")
        self.slider_rh = QSlider(Qt.Horizontal)
        self.slider_rh.setMinimumWidth(270)
        self.slider_rh.setMaximumWidth(270)
        self.slider_rh.setRange(1, nbins)
        self.slider_rh.setValue(nbins/2)
        self.slider_rh.setTracking(True)
        #self.slider_rh.setTickPosition(QSlider.TicksBothSides)
        self.slider_rh.valueChanged.connect(self.on_draw)#int

        slider_label_y1 = QLabel("Source Position y1")
        self.slider_y1 = QSlider(Qt.Horizontal)
        self.slider_y1.setMinimumWidth(270)
        self.slider_y1.setMaximumWidth(270)
        self.slider_y1.setRange(1, nbins)
        self.slider_y1.setValue(nbins/2)
        self.slider_y1.setTracking(True)
        #self.slider_y1.setTickPosition(QSlider.TicksBothSides)
        self.slider_y1.valueChanged.connect(self.on_draw)#int

        slider_label_y2 = QLabel("Source Position y2")
        self.slider_y2 = QSlider(Qt.Horizontal)
        self.slider_y2.setMinimumWidth(270)
        self.slider_y2.setMaximumWidth(270)
        self.slider_y2.setRange(1, nbins)
        self.slider_y2.setValue(nbins/2)
        self.slider_y2.setTracking(True)
        #self.slider_y2.setTickPosition(QSlider.TicksBothSides)
        self.slider_y2.valueChanged.connect(self.on_draw)#int

        slider_label_qs = QLabel("Ellipticity S")
        self.slider_qs = QSlider(Qt.Horizontal)
        self.slider_qs.setMinimumWidth(270)
        self.slider_qs.setMaximumWidth(270)
        self.slider_qs.setRange(1, nbins)
        self.slider_qs.setValue(nbins/2)
        self.slider_qs.setTracking(True)
        #self.slider_qs.setTickPosition(QSlider.TicksBothSides)
        self.slider_qs.valueChanged.connect(self.on_draw)#int

        slider_label_sa = QLabel("Orintation S")
        self.slider_sa = QSlider(Qt.Horizontal)
        self.slider_sa.setMinimumWidth(270)
        self.slider_sa.setMaximumWidth(270)
        self.slider_sa.setRange(1, nbins)
        self.slider_sa.setValue(nbins/3)
        self.slider_sa.setTracking(True)
        #self.slider_sa.setTickPosition(QSlider.TicksBothSides)
        self.slider_sa.valueChanged.connect(self.on_draw)#int

        #
        # Layout with box sizers
        #
        hbox = QHBoxLayout()

        for w in [self.textbox, self.cb_cosmo, self.cb_models, self.cb_codes, self.draw_button]:
            hbox.addWidget(w)
            hbox.setAlignment(w, Qt.AlignVCenter)

        vbox = QVBoxLayout()
        vbox.addWidget(self.mpl_toolbar)
        vbox.addWidget(self.canvas)

        hbox_s1 = QHBoxLayout()
        for w in [self.slider_label_re, self.slider_re, slider_label_rh, self.slider_rh]:
            hbox_s1.addWidget(w)
            hbox_s1.setAlignment(w, Qt.AlignLeft)

        hbox_s2 = QHBoxLayout()
        for w in [slider_label_x1, self.slider_x1, slider_label_y1, self.slider_y1]:
            hbox_s2.addWidget(w)
            hbox_s2.setAlignment(w, Qt.AlignLeft)

        hbox_s3 = QHBoxLayout()
        for w in [slider_label_x2, self.slider_x2, slider_label_y2, self.slider_y2]:
            hbox_s3.addWidget(w)
            hbox_s3.setAlignment(w, Qt.AlignLeft)

        hbox_s4 = QHBoxLayout()
        for w in [slider_label_ql, self.slider_ql, slider_label_qs, self.slider_qs]:
            hbox_s4.addWidget(w)
            hbox_s4.setAlignment(w, Qt.AlignLeft)

        hbox_s5 = QHBoxLayout()
        for w in [slider_label_la, self.slider_la, slider_label_sa, self.slider_sa]:
            hbox_s5.addWidget(w)
            hbox_s5.setAlignment(w, Qt.AlignLeft)

        vbox.addLayout(hbox)
        vbox.addLayout(hbox_s1)
        vbox.addLayout(hbox_s2)
        vbox.addLayout(hbox_s3)
        vbox.addLayout(hbox_s4)
        vbox.addLayout(hbox_s5)

        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)

    def optimization_actions(self):
        print "Optimizing..."
        lpar_new, spar_new = ff.optimize_pars(self.xi1, self.xi2,
                                              self.lgals.ravel(),
                                              self.lpar, self.spar)

        g_limage,g_source,mua,yi1,yi2 = ll.lensed_images(self.xi1,self.xi2,lpar_new,spar_new)
        print "Done."

        self.slider_re.setValue(round(lpar_new[2]/2.0*nbins))
        self.slider_x1.setValue(round((lpar_new[0]+self.bsz/2.0-self.bsz/nbins/2.0)/self.bsz*nbins))
        self.slider_x2.setValue(round((lpar_new[1]+self.bsz/2.0-self.bsz/nbins/2.0)/self.bsz*nbins))
        self.slider_ql.setValue(round(lpar_new[4]*nbins))
        self.slider_la.setValue(round(lpar_new[5]/360.0*nbins))

        self.slider_rh.setValue(round(spar_new[2]/0.5*nbins))
        self.slider_y1.setValue(round((spar_new[0]+self.bsz/8.0-self.bsz/nbins/8.0)/(self.bsz/4.0)*nbins))
        self.slider_y2.setValue(round((spar_new[1]+self.bsz/8.0-self.bsz/nbins/8.0)/(self.bsz/4.0)*nbins))
        self.slider_qs.setValue(round(spar_new[4]*nbins))
        self.slider_sa.setValue(round(spar_new[5]/360.0*nbins))

        levels = [0.5,]
        levels_imgs = [0.0,0.08,0.1,0.2,0.3,0.4,0.5]
        self.axesa.clear()
        # pl.xticks([-2.0, 1.0, 0.0, 1.0, 2.0])
        # pl.yticks([-2.0, 1.0, 0.0, 1.0, 2.0])
        self.axesa.set_xlim(-self.bsz/8.0,self.bsz/8.0)
        self.axesa.set_ylim(-self.bsz/8.0,self.bsz/8.0)
        self.axesa.contourf(self.xi1,self.xi2,self.sgals,levels_imgs,cmap=pl.cm.gray)
        self.axesa.contour(self.xi1,self.xi2,g_source,levels,colors=('deepskyblue'))
        self.axesa.contour(yi1,yi2,mua,0,colors=('g'),linewidths = 2.0)

        self.axesb.clear()
        # pl.xticks([-0.2,-0.1,0.0,0.1,0.2])
        # pl.yticks([-0.2,-0.1,0.0,0.1,0.2])
        self.axesb.set_xlim(-self.bsz/2.0,self.bsz/2.0)
        self.axesb.set_ylim(-self.bsz/2.0,self.bsz/2.0)
        self.axesb.contourf(self.xi1,self.xi2,self.lgals,levels_imgs,cmap=pl.cm.gray)
        self.axesb.contour(self.xi1,self.xi2,g_limage,levels,colors=('deepskyblue'))
        self.axesb.contour(self.xi1,self.xi2,mua,colors=('r'),linewidths = 2.0)

        self.canvas.draw()

        print "opt----", spar_new

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
    app.setWindowIcon(QIcon("./icon.png"))
    form = AppForm()
    form.resize(800, 700)
    form.setWindowTitle('Qoopla')
    form.show()
    app.exec_()

if __name__ == "__main__":
    main()
