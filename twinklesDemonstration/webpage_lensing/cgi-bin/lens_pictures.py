#!/usr/bin/env python 

import numpy as np
from PIL import Image
import scipy.ndimage as scn
import sys
import cgi 
import uuid
import urllib, cStringIO
import cgitb; cgitb.enable() 

print "Content-type: text/html" 
print 


form = cgi.FieldStorage()



def point_alphas(xi1, xi2, xc1=0.0, xc2=0.0, b=100.0):
    x1 = xi1 - xc1
    x2 = xi2 - xc2
    wx = np.sqrt(x1 * x1 + x2 * x2 + 0.00000001)
    al1 = b * x1 / wx
    al2 = b * x2 / wx
    return al1, al2


def ray_tracing(input_image, yi1, yi2, nnx1, nnx2):
    y1 = yi1 + nnx1 / 2.0
    y2 = yi2 + nnx2 / 2.0
    res = scn. map_coordinates(input_image, [y1.flatten(), y2.flatten()], order=1)
    return res.reshape((nnx1, nnx2))


URL = form.getvalue("url","empty.jpg"); 

xc1 = float(form.getvalue("xc1",0)); 
xc2 = float(form.getvalue("xc2",0)); 
b = float(form.getvalue("b",100.)); 

file_pic = cStringIO.StringIO(urllib.urlopen(URL).read())
img_in = Image.open(file_pic)
img_in = np.array(img_in, dtype=np.double)

nnx1, nnx2 = np.shape(img_in[:, :, 0])
xi1, xi2 = np.mgrid[0:nnx1 - 1:nnx1 * 1j, 0:nnx2 - 1:nnx2 * 1j]
xi1 = xi1 - nnx1 / 2.0
xi2 = xi2 - nnx2 / 2.0

ai1, ai2 = point_alphas(xi1, xi2,xc1,xc2,b)

yi1 = xi1 - ai1
yi2 = xi2 - ai2

img_out = img_in * 0.0

img_out[:, :, 0] = ray_tracing(img_in[:, :, 0], yi1, yi2, nnx1, nnx2)
img_out[:, :, 1] = ray_tracing(img_in[:, :, 1], yi1, yi2, nnx1, nnx2)
img_out[:, :, 2] = ray_tracing(img_in[:, :, 2], yi1, yi2, nnx1, nnx2)

img_out = np.uint8(img_out)
img = Image.fromarray(img_out)
imgname = "../lens/imgs/" + str(uuid.uuid4())  + ".jpg"
img.save(imgname, "JPEG")


print "<html><body>"
print "<img src=" + imgname +">" 
print "</body></html>"
print 
