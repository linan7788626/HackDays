from pyqtgraph.Qt import QtGui,QtCore
import numpy as np
import pyqtgraph as pg
import sys

#data creation
app=QtGui.QApplication([])
x = np.linspace(0,6.28,30)
y = x[:]
xx,yy = np.meshgrid(x,y)
z = np.sin(xx)+np.cos(yy)

win = pg.GraphicsWindow()
win.setWindowTitle('esempio isocurve')
vb = win.addViewBox()
img = pg.ImageItem(z)
vb.addItem(img) #visualizes img in viewbox vb
vb.setAspectLocked() #set img proportions (?)

c = pg.IsocurveItem(data=z,level=2,pen='r')
c.setParentItem(img) #set img as parent of the isocurve item c (?)
c.setZValue(10) #I DO NOT KNOW THE MEANING OF THIS, JUST COPIED FROM THE EXAMPLE isocurve.py

vb.addItem(c)

win.show()
sys.exit(app.exec_())
