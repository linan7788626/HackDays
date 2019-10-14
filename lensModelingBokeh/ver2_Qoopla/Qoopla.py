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

bsz = 5.0 # (arcsec)
nnn = 400
dsx = bsz/nnn

nbins = 200
nsbox = 100

re_max = 3.0
bpl = bsz
dpl = bpl/nbins
xbc1 = 0.0
xbc2 = 0.0
x1_min = -bsz/2.0+dpl/2.0+xbc1/2.0
x1_max =  bsz/2.0-dpl/2.0-xbc1/2.0
x2_min = -bsz/2.0+dpl/2.0+xbc2/2.0
x2_max =  bsz/2.0-dpl/2.0-xbc2/2.0
ql_max = 1.0
la_max = 360.

rh_max = 0.5
bps = bsz/4.0
dps = bps/nbins
ybc1 = 0.0
ybc2 = 0.0
y1_min = -bps/2.0+dps/2.0+ybc1/2.0
y1_max =  bps/2.0-dps/2.0-ybc1/2.0
y2_min = -bps/2.0+dps/2.0+ybc2/2.0
y2_max =  bps/2.0-dps/2.0-ybc2/2.0

qs_max = 1.0
sa_max = 360.

st_re = nbins/2
st_x1 = nbins/2
st_x2 = nbins/2
st_rc = 0*nbins
st_ql = nbins*3./4
st_la = nbins/3
st_rh = nbins/2
st_y1 = nbins/2
st_y2 = nbins/2
st_am = 0*nbins
st_qs = nbins/2
st_sa = nbins/3

slider_len = 270
sbox_len = 70

class AppForm(QMainWindow):
    def __init__(self, parent=None):
        self.xi1 = np.linspace(x1_min, x1_max, nnn)
        self.xi2 = np.linspace(x2_min, x2_max, nnn)
        self.xi1,self.xi2 = np.meshgrid(self.xi1,self.xi2)

        self.file_name = "./ering.jpg"
        a_tmp = Image.open(self.file_name)
        a_tmp = np.array(a_tmp)
        self.lgals = np.flipud(a_tmp[:,:,0]/256.0)

        # self.file_name = "./ds9_images_500.png"
        # a_tmp = Image.open(self.file_name)
        # a_tmp = np.array(a_tmp)
        # self.lgals = np.flipud(a_tmp[:,:,0]/256.0)

        # self.srcs_name = "./ds9_srcs_500.png"
        # s_tmp = Image.open(self.srcs_name)
        # s_tmp = np.array(s_tmp)
        # self.sgals = np.flipud(s_tmp[:,:,0]/256.0)

        self.lpar = np.asarray([(st_x1*dpl + x1_min),
                                (st_x2*dpl + x2_min),
                                st_re*re_max/nbins, st_rc*0.0, st_ql*ql_max/nbins, st_la*la_max/nbins])
        self.spar = np.asarray([(st_y1*dps + y1_min),
                                (st_y2*dps + y2_min),
                                st_rh*rh_max/nbins, st_sa+1.0, st_qs*qs_max/nbins, st_sa*sa_max/nbins])
        #----------------------------------------------------------------------
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('Interactive Lens Modeling')

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()

        self.textbox.setText('67.26 0.32 0.68')
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

        l_xcen = self.slider_x1.value()*dpl+x1_min  	# x position of center (also try (0.0,0.14)
        l_ycen = self.slider_x2.value()*dpl+x2_min 	# y position of center
        l_re = self.slider_re.value()*re_max/nbins   # Einstein radius of lens.
        l_rc = 0.0   # Core size of lens (in units of Einstein radius).
        l_axrat = self.slider_ql.value()*ql_max/nbins   # Axis ratio of lens.
        l_pa = self.slider_la.value()*la_max/nbins  # Orintation of lens.

        self.lpar = np.asarray([l_xcen,l_ycen,l_re,l_rc,l_axrat,l_pa])
        #----------------------------------------------------------------------
        g_amp = 1.0   	# peak brightness value
        g_sig = self.slider_rh.value()*rh_max/nbins  	# Gaussian "sigma" (i.e., size)
        g_xcen =  self.slider_y1.value()*dps+y1_min 	# x position of center (also try (0.0,0.14)
        g_ycen =  self.slider_y2.value()*dps+y2_min 	# y position of center
        g_axrat =  self.slider_qs.value()*qs_max/nbins 	# minor-to-major axis ratio
        g_pa =  self.slider_sa.value()*sa_max/nbins      	# major-axis position angle (degrees) c.c.w. from x axis

        self.spar = np.asarray([g_xcen,g_ycen,g_sig,g_amp,g_axrat,g_pa])

        #----------------------------------------------------------------------
        g_limage,g_source,mua,yi1,yi2 = ll.lensed_images(self.xi1,self.xi2,self.lpar,self.spar)
        #--------------------------lens images contour------------------------

        levels = [0.5,]
        levels_imgs = [0.0,0.08,0.1,0.2,0.3,0.4,0.5]

        self.axesa.clear()
        self.axesa.set_xlim(y1_min,y1_max)
        self.axesa.set_ylim(y2_min,y2_max)
        self.axesa.set_xticks([ybc1-0.4, ybc1-0.2, ybc1, ybc1+0.2, ybc1+0.4])
        self.axesa.set_yticks([ybc2-0.4, ybc2-0.2, ybc2, ybc2+0.2, ybc2+0.4])
        # self.axesa.contourf(self.xi1,self.xi2,self.sgals,levels_imgs,cmap=pl.cm.gray)
        self.axesa.contour(self.xi1,self.xi2,g_source,levels,colors=('deepskyblue'))
        self.axesa.contour(yi1,yi2,mua,0,colors=('r'),linewidths = 2.0)

        self.axesb.clear()
        self.axesb.set_xlim(x1_min,x1_max)
        self.axesb.set_ylim(x2_min,x2_max)
        self.axesb.set_xticks([xbc1-2, xbc1-1, xbc1, xbc1+1, xbc1+2])
        self.axesb.set_yticks([xbc2-2, xbc2-1, xbc2, xbc2+1, xbc2+2])
        self.axesb.contourf(self.xi1,self.xi2,self.lgals,levels_imgs,cmap=pl.cm.gray)
        self.axesb.contour(self.xi1,self.xi2,g_limage,levels,colors=('deepskyblue'))
        self.axesb.contour(self.xi1,self.xi2,mua,colors=('r'),linewidths = 2.0)

        self.canvas.draw()

#-------
    def update_slider_re(self):
        spinbox_value = round(self.sbox_re.value() / re_max * nbins)
        self.slider_re.setSliderPosition(spinbox_value)

    def update_sbox_re(self, value):
        self.sbox_re.setValue(float(value) * re_max / nbins)
#-------
    def update_slider_x1(self):
        spinbox_value = round((self.sbox_x1.value()-x1_min) / (x1_max-x1_min) * nbins)
        self.slider_x1.setSliderPosition(spinbox_value)

    def update_sbox_x1(self, value):
        self.sbox_x1.setValue(float(value) * (x1_max-x1_min) / nbins + x1_min)
#-------
    def update_slider_x2(self):
        spinbox_value = round((self.sbox_x2.value()-x2_min) / (x2_max-x2_min) * nbins)
        self.slider_x2.setSliderPosition(spinbox_value)

    def update_sbox_x2(self, value):
        self.sbox_x2.setValue(float(value) * (x2_max-x2_min) / nbins + x2_min)
#-------
    def update_slider_ql(self):
        spinbox_value = round(self.sbox_ql.value() / ql_max * nbins)
        self.slider_ql.setSliderPosition(spinbox_value)

    def update_sbox_ql(self, value):
        self.sbox_ql.setValue(float(value) * ql_max / nbins)
#-------
    def update_slider_la(self):
        spinbox_value = round(self.sbox_la.value() / la_max * nbins)
        self.slider_la.setSliderPosition(spinbox_value)

    def update_sbox_la(self, value):
        self.sbox_la.setValue(float(value) * la_max / nbins)
#-------
    def update_slider_rh(self):
        spinbox_value = round(self.sbox_rh.value() / rh_max * nbins)
        self.slider_rh.setSliderPosition(spinbox_value)

    def update_sbox_rh(self, value):
        self.sbox_rh.setValue(float(value) * rh_max / nbins)
#-------
    def update_slider_y1(self):
        spinbox_value = round((self.sbox_y1.value()-y1_min) / (y1_max-y1_min) * nbins)
        self.slider_y1.setSliderPosition(spinbox_value)

    def update_sbox_y1(self, value):
        self.sbox_y1.setValue(float(value) * (y1_max-y1_min) / nbins + y1_min)
#-------
    def update_slider_y2(self):
        spinbox_value = round((self.sbox_y2.value()-y2_min) / (y2_max-y2_min) * nbins)
        self.slider_y2.setSliderPosition(spinbox_value)

    def update_sbox_y2(self, value):
        self.sbox_y2.setValue(float(value) * (y2_max-y2_min) / nbins + y2_min)
#-------
    def update_slider_qs(self):
        spinbox_value = round(self.sbox_qs.value() / qs_max * nbins)
        self.slider_qs.setSliderPosition(spinbox_value)

    def update_sbox_qs(self, value):
        self.sbox_qs.setValue(float(value) * qs_max / nbins)
#-------
    def update_slider_sa(self):
        spinbox_value = round(self.sbox_sa.value() / sa_max * nbins)
        self.slider_sa.setSliderPosition(spinbox_value)

    def update_sbox_sa(self, value):
        self.sbox_sa.setValue(float(value) * sa_max / nbins)
#-------

    def create_main_frame(self):
        self.main_frame = QWidget()
        #self.main_frame.resize(800, 400)
        self.main_frame.setFixedSize(820, 700)

        self.dpi = 80
        self.fig = pl.figure(num=None,figsize=(10,5), dpi=self.dpi, facecolor='w', edgecolor='k')
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)

        self.axesa = pl.axes([0.07,0.1,0.4,0.88])
        self.axesb = pl.axes([0.57,0.1,0.4,0.88])

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


#-------
        self.slider_label_re = QLabel("Re")
        self.slider_re = QSlider(Qt.Horizontal)
        self.slider_re.setFocusPolicy(Qt.StrongFocus)
        self.slider_re.setMinimumWidth(slider_len)
        self.slider_re.setMaximumWidth(slider_len)
        self.slider_re.setRange(0, nbins)
        self.slider_re.setValue(st_re)
        self.slider_re.setTracking(True)
        self.slider_re.valueChanged.connect(self.on_draw)#int

        self.sbox_re = QDoubleSpinBox()
        self.sbox_re.setSingleStep(re_max/nsbox)
        self.sbox_re.setFocusPolicy(Qt.StrongFocus)
        self.sbox_re.setMinimumWidth(sbox_len)
        self.sbox_re.setMaximumWidth(sbox_len)
        self.sbox_re.setRange(0, re_max)
        self.sbox_re.setValue(re_max*st_re*1.0/nbins)

        self.slider_re.valueChanged.connect(self.update_sbox_re)
        self.sbox_re.editingFinished.connect(self.update_slider_re)
#-------
        self.slider_label_x1 = QLabel("x1")
        self.slider_x1 = QSlider(Qt.Horizontal)
        self.slider_x1.setFocusPolicy(Qt.StrongFocus)
        self.slider_x1.setMinimumWidth(slider_len)
        self.slider_x1.setMaximumWidth(slider_len)
        self.slider_x1.setRange(0, nbins)
        self.slider_x1.setValue(st_x1)
        self.slider_x1.setTracking(True)
        self.slider_x1.valueChanged.connect(self.on_draw)#int

        self.sbox_x1 = QDoubleSpinBox()
        self.sbox_x1.setSingleStep((x1_max-x1_min)/nsbox)
        self.sbox_x1.setFocusPolicy(Qt.StrongFocus)
        self.sbox_x1.setMinimumWidth(sbox_len)
        self.sbox_x1.setMaximumWidth(sbox_len)
        self.sbox_x1.setRange(x1_min, x1_max)
        self.sbox_x1.setValue(st_x1*dpl+x1_min)

        self.slider_x1.valueChanged.connect(self.update_sbox_x1)
        self.sbox_x1.editingFinished.connect(self.update_slider_x1)
#-------
        self.slider_label_x2 = QLabel("x2")
        self.slider_x2 = QSlider(Qt.Horizontal)
        self.slider_x2.setFocusPolicy(Qt.StrongFocus)
        self.slider_x2.setMinimumWidth(slider_len)
        self.slider_x2.setMaximumWidth(slider_len)
        self.slider_x2.setRange(0, nbins)
        self.slider_x2.setValue(st_x2)
        self.slider_x2.setTracking(True)
        self.slider_x2.valueChanged.connect(self.on_draw)#int

        self.sbox_x2 = QDoubleSpinBox()
        self.sbox_x2.setSingleStep((x2_max-x2_min)/nsbox)
        self.sbox_x2.setFocusPolicy(Qt.StrongFocus)
        self.sbox_x2.setMinimumWidth(sbox_len)
        self.sbox_x2.setMaximumWidth(sbox_len)
        self.sbox_x2.setRange(x2_min, x2_max)
        self.sbox_x2.setValue(st_x2*dpl+x2_min)

        self.slider_x2.valueChanged.connect(self.update_sbox_x2)
        self.sbox_x2.editingFinished.connect(self.update_slider_x2)
#-------
        self.slider_label_ql = QLabel("ql")
        self.slider_ql = QSlider(Qt.Horizontal)
        self.slider_ql.setFocusPolicy(Qt.StrongFocus)
        self.slider_ql.setMinimumWidth(slider_len)
        self.slider_ql.setMaximumWidth(slider_len)
        self.slider_ql.setRange(0, nbins)
        self.slider_ql.setValue(st_ql)
        self.slider_ql.setTracking(True)
        self.slider_ql.valueChanged.connect(self.on_draw)#int

        self.sbox_ql = QDoubleSpinBox()
        self.sbox_ql.setSingleStep(ql_max/nsbox)
        self.sbox_ql.setFocusPolicy(Qt.StrongFocus)
        self.sbox_ql.setMinimumWidth(sbox_len)
        self.sbox_ql.setMaximumWidth(sbox_len)
        self.sbox_ql.setRange(0, ql_max)
        self.sbox_ql.setValue(ql_max*st_ql*1.0/nbins)

        self.slider_ql.valueChanged.connect(self.update_sbox_ql)
        self.sbox_ql.editingFinished.connect(self.update_slider_ql)
#-------
        self.slider_label_la = QLabel("la")
        self.slider_la = QSlider(Qt.Horizontal)
        self.slider_la.setFocusPolicy(Qt.StrongFocus)
        self.slider_la.setMinimumWidth(slider_len)
        self.slider_la.setMaximumWidth(slider_len)
        self.slider_la.setRange(0, nbins)
        self.slider_la.setValue(st_la)
        self.slider_la.setTracking(True)
        self.slider_la.valueChanged.connect(self.on_draw)#int

        self.sbox_la = QDoubleSpinBox()
        self.sbox_la.setSingleStep(la_max/nsbox)
        self.sbox_la.setFocusPolicy(Qt.StrongFocus)
        self.sbox_la.setMinimumWidth(sbox_len)
        self.sbox_la.setMaximumWidth(sbox_len)
        self.sbox_la.setRange(0, la_max)
        self.sbox_la.setValue(la_max*st_la*1.0/nbins)

        self.slider_la.valueChanged.connect(self.update_sbox_la)
        self.sbox_la.editingFinished.connect(self.update_slider_la)
#-------
        self.slider_label_rh = QLabel("Rh")
        self.slider_rh = QSlider(Qt.Horizontal)
        self.slider_rh.setFocusPolicy(Qt.StrongFocus)
        self.slider_rh.setMinimumWidth(slider_len)
        self.slider_rh.setMaximumWidth(slider_len)
        self.slider_rh.setRange(0, nbins)
        self.slider_rh.setValue(st_rh)
        self.slider_rh.setTracking(True)
        self.slider_rh.valueChanged.connect(self.on_draw)#int

        self.sbox_rh = QDoubleSpinBox()
        self.sbox_rh.setSingleStep(rh_max/nsbox)
        self.sbox_rh.setFocusPolicy(Qt.StrongFocus)
        self.sbox_rh.setMinimumWidth(sbox_len)
        self.sbox_rh.setMaximumWidth(sbox_len)
        self.sbox_rh.setRange(0, rh_max)
        self.sbox_rh.setValue(rh_max*st_rh*1.0/nbins)

        self.slider_rh.valueChanged.connect(self.update_sbox_rh)
        self.sbox_rh.editingFinished.connect(self.update_slider_rh)
#-------
        self.slider_label_y1 = QLabel("y1")
        self.slider_y1 = QSlider(Qt.Horizontal)
        self.slider_y1.setFocusPolicy(Qt.StrongFocus)
        self.slider_y1.setMinimumWidth(slider_len)
        self.slider_y1.setMaximumWidth(slider_len)
        self.slider_y1.setRange(0, nbins)
        self.slider_y1.setValue(st_y1)
        self.slider_y1.setTracking(True)
        self.slider_y1.valueChanged.connect(self.on_draw)#int

        self.sbox_y1 = QDoubleSpinBox()
        self.sbox_y1.setSingleStep((y1_max-y1_min)/nsbox)
        self.sbox_y1.setFocusPolicy(Qt.StrongFocus)
        self.sbox_y1.setMinimumWidth(sbox_len)
        self.sbox_y1.setMaximumWidth(sbox_len)
        self.sbox_y1.setRange(y1_min, y1_max)
        self.sbox_y1.setValue(st_y1*dps+y1_min)

        self.slider_y1.valueChanged.connect(self.update_sbox_y1)
        self.sbox_y1.editingFinished.connect(self.update_slider_y1)
#-------
        self.slider_label_y2 = QLabel("y2")
        self.slider_y2 = QSlider(Qt.Horizontal)
        self.slider_y2.setFocusPolicy(Qt.StrongFocus)
        self.slider_y2.setMinimumWidth(slider_len)
        self.slider_y2.setMaximumWidth(slider_len)
        self.slider_y2.setRange(0, nbins)
        self.slider_y2.setValue(st_y2)
        self.slider_y2.setTracking(True)
        self.slider_y2.valueChanged.connect(self.on_draw)#int

        self.sbox_y2 = QDoubleSpinBox()
        self.sbox_y2.setSingleStep((y2_max-y2_min)/nsbox)
        self.sbox_y2.setFocusPolicy(Qt.StrongFocus)
        self.sbox_y2.setMinimumWidth(sbox_len)
        self.sbox_y2.setMaximumWidth(sbox_len)
        self.sbox_y2.setRange(y2_min, y2_max)
        self.sbox_y2.setValue(st_y2*dps+y2_min)

        self.slider_y2.valueChanged.connect(self.update_sbox_y2)
        self.sbox_y2.editingFinished.connect(self.update_slider_y2)
#-------
        self.slider_label_qs = QLabel("qs")
        self.slider_qs = QSlider(Qt.Horizontal)
        self.slider_qs.setFocusPolicy(Qt.StrongFocus)
        self.slider_qs.setMinimumWidth(slider_len)
        self.slider_qs.setMaximumWidth(slider_len)
        self.slider_qs.setRange(0, nbins)
        self.slider_qs.setValue(st_qs)
        self.slider_qs.setTracking(True)
        self.slider_qs.valueChanged.connect(self.on_draw)#int

        self.sbox_qs = QDoubleSpinBox()
        self.sbox_qs.setSingleStep(qs_max/nsbox)
        self.sbox_qs.setFocusPolicy(Qt.StrongFocus)
        self.sbox_qs.setMinimumWidth(sbox_len)
        self.sbox_qs.setMaximumWidth(sbox_len)
        self.sbox_qs.setRange(0.0, qs_max)
        self.sbox_qs.setValue(qs_max*st_qs*1.0/nbins)

        self.slider_qs.valueChanged.connect(self.update_sbox_qs)
        self.sbox_qs.editingFinished.connect(self.update_slider_qs)
#-------
        self.slider_label_sa = QLabel("sa")
        self.slider_sa = QSlider(Qt.Horizontal)
        self.slider_sa.setFocusPolicy(Qt.StrongFocus)
        self.slider_sa.setMinimumWidth(slider_len)
        self.slider_sa.setMaximumWidth(slider_len)
        self.slider_sa.setRange(0, nbins)
        self.slider_sa.setValue(st_sa)
        self.slider_sa.setTracking(True)
        self.slider_sa.valueChanged.connect(self.on_draw)#int

        self.sbox_sa = QDoubleSpinBox()
        self.sbox_sa.setSingleStep(sa_max/nsbox)
        self.sbox_sa.setFocusPolicy(Qt.StrongFocus)
        self.sbox_sa.setMinimumWidth(sbox_len)
        self.sbox_sa.setMaximumWidth(sbox_len)
        self.sbox_sa.setRange(0.0, sa_max)
        self.sbox_sa.setValue(sa_max*st_sa*1.0/nbins)

        self.slider_sa.valueChanged.connect(self.update_sbox_sa)
        self.sbox_sa.editingFinished.connect(self.update_slider_sa)
#-------

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
        for w in [self.slider_label_re, self.sbox_re, self.slider_re, self.slider_label_rh, self.sbox_rh, self.slider_rh]:
            hbox_s1.addWidget(w)
            hbox_s1.setAlignment(w, Qt.AlignLeft)

        hbox_s2 = QHBoxLayout()
        for w in [self.slider_label_x1, self.sbox_x1, self.slider_x1, self.slider_label_y1, self.sbox_y1, self.slider_y1]:
            hbox_s2.addWidget(w)
            hbox_s2.setAlignment(w, Qt.AlignLeft)

        hbox_s3 = QHBoxLayout()
        for w in [self.slider_label_x2, self.sbox_x2, self.slider_x2, self.slider_label_y2, self.sbox_y2, self.slider_y2]:
            hbox_s3.addWidget(w)
            hbox_s3.setAlignment(w, Qt.AlignLeft)

        hbox_s4 = QHBoxLayout()
        for w in [self.slider_label_ql, self.sbox_ql, self.slider_ql, self.slider_label_qs, self.sbox_qs, self.slider_qs]:
            hbox_s4.addWidget(w)
            hbox_s4.setAlignment(w, Qt.AlignLeft)

        hbox_s5 = QHBoxLayout()
        for w in [self.slider_label_la, self.sbox_la, self.slider_la, self.slider_label_sa, self.sbox_sa, self.slider_sa]:
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

        self.slider_re.setValue(round(lpar_new[2]/re_max*nbins))
        self.slider_x1.setValue(round((lpar_new[0]-x1_min)/dpl))
        self.slider_x2.setValue(round((lpar_new[1]-x2_min)/dpl))
        self.slider_ql.setValue(round(lpar_new[4]/ql_max*nbins))
        self.slider_la.setValue(round(lpar_new[5]/la_max*nbins))

        self.slider_rh.setValue(round(spar_new[2]/rh_max*nbins))
        self.slider_y1.setValue(round((spar_new[0]-y1_min)/dps))
        self.slider_y2.setValue(round((spar_new[1]-y2_min)/dps))
        self.slider_qs.setValue(round(spar_new[4]/qs_max*nbins))
        self.slider_sa.setValue(round(spar_new[5]/sa_max*nbins))

        levels = [0.5,]
        levels_imgs = [0.0,0.08,0.1,0.2,0.3,0.4,0.5]
        self.axesa.clear()
        self.axesa.set_xlim(y1_min,y1_max)
        self.axesa.set_ylim(y2_min,y2_max)
        self.axesa.set_xticks([ybc1-0.4, ybc1-0.2, ybc1, ybc1+0.2, ybc1+0.4])
        self.axesa.set_yticks([ybc2-0.4, ybc2-0.2, ybc2, ybc2+0.2, ybc2+0.4])
        # self.axesa.contourf(self.xi1,self.xi2,self.sgals,levels_imgs,cmap=pl.cm.gray)
        self.axesa.contour(self.xi1,self.xi2,g_source,levels,colors=('deepskyblue'))
        self.axesa.contour(yi1,yi2,mua,0,colors=('r'),linewidths = 2.0)

        self.axesb.clear()
        self.axesb.set_xlim(x1_min,x1_max)
        self.axesb.set_ylim(x2_min,x2_max)
        self.axesb.set_xticks([xbc1-2.0,xbc1-1.0,xbc1,xbc1+1.0,xbc1+2.0])
        self.axesb.set_yticks([xbc2-2.0,xbc2-1.0,xbc2,xbc2+1.0,xbc2+2.0])
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
