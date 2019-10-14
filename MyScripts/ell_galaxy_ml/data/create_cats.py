import numpy as np


def creat_cats(nimages):

    sigma = np.random.random(nimages)*1.5+0.5
    ell = np.random.random(nimages)*0.5+0.5
    pha = np.random.random(nimages)*360.0

    for i in xrange(nimages):
        print sigma[i], ell[i], pha[i]

    return 0


if __name__ == '__main__':
    nimgs = 1000
    creat_cats(nimgs)
