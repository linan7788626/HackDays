import numpy as np
import astropy.io.fits as pyfits

name1="test.1"
name2="test.2"

a = np.loadtxt(name1)

b = np.loadtxt(name2)

oneone = a[:,0]
onetwo = a[:,1]

twoone = b[:,0]
twotwo = b[:,1]

oneone = oneone.astype('int')
onetwo = onetwo.astype('int')
twoone = twoone.astype('int')
twotwo = twotwo.astype('int')

idx1 = onetwo
idx2 = twotwo

tt = np.where(onetwo == 107) 
print oneone[tt]


#100 102
#140 107
#130 109
#133 115

#200 122
#250 188
#250 287
#983 443
#922 107
#488 432