import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import numpy as np

pg.setConfigOptions(imageAxisOrder='row-major')

## Create image to display
arr = np.ones((100, 100), dtype=float)
arr[45:55, 45:55] = 0
arr[25, :] = 5
arr[:, 25] = 5
arr[75, :] = 5
arr[:, 75] = 5
arr[50, :] = 10
arr[:, 50] = 10
arr += np.sin(np.linspace(0, 20, 100)).reshape(1, 100)
arr += np.random.normal(size=(100,100))

## create GUI
app = QtGui.QApplication([])
w = pg.GraphicsWindow(size=(1000,800), border=True)
w.setWindowTitle('pyqtgraph example: ROI Examples')

text = """Data Selection From Image.<br>\n
Drag an ROI or its handles to update the selected image.<br>
Hold CTRL while dragging to snap to pixel boundaries<br>
and 15-degree rotation angles.
"""
w1 = w.addLayout(row=0, col=0)
label1 = w1.addLabel(text, row=0, col=0)
v1a = w1.addViewBox(row=1, col=0, lockAspect=True)
v1b = w1.addViewBox(row=1, col=1, lockAspect=True)
img1a = pg.ImageItem(arr)
v1a.addItem(img1a)
img1b = pg.ImageItem()
v1b.addItem(img1b)
v1a.disableAutoRange('xy')
v1b.disableAutoRange('xy')
v1a.autoRange()
v1b.autoRange()

rois = []
rois.append(pg.EllipseROI([30, 40], [30, 20], pen=(3,9)))
rois.append(pg.EllipseROI([60, 10], [30, 20], pen=(3,9)))

def update(roi):
    img1b.setImage(roi.getArrayRegion(arr, img1a), levels=(0, arr.max()))
    v1b.autoRange()

for roi in rois:
    roi.sigRegionChanged.connect(update)
    v1a.addItem(roi)

update(rois[-1])


if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
