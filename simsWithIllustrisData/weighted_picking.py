import numpy as np

a = np.random.random(1000)
b = np.random.normal(1.0,1.0,1000)
c = np.random.normal(0.7,0.2,1000)

print len(a)
print len(b)
print len(c)

a0 = 0.4
b0 = 1.2
c0 = 0.6

xsqa = (a-a0)**2.0
xsqb = (b-b0)**2.0
xsqc = (c-c0)**2.0

w1 = 1.0
w2 = 1.0
w3 = 1.0

xv = xsqa/w1+xsqb/w2+xsqc/w3


lkhd = np.exp(-xv)

idx = lkhd==lkhd.max()

print a0,b0,c0
print a[idx],b[idx],c[idx]
