import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

data=np.random.random([25,4])
data = data*100
len=np.sqrt(data[:,0].size)
x=np.reshape(data[:,0],(len,len))
y=np.reshape(data[:,1],(len,len))
z=np.reshape(data[:,3],(len,len))

l=plt.contourf(x,y,z,np.linspace(0,100,255))
contour_axis = plt.gca()

axmax = plt.axes([0.25, 0.1, 0.65, 0.03])  #slider location and size
axmin  = plt.axes([0.25, 0.15, 0.65, 0.03])
smax = Slider(axmax, 'Max',0, 100, 50)      #slider properties
smin = Slider(axmin, 'Min', 0, 100, 0)


def update(val):
    contour_axis.clear()
    contour_axis.contourf(x,y,z,np.linspace(smin.val,smax.val,255))
    plt.draw()
smax.on_changed(update)
smin.on_changed(update)

plt.show()
