"""
This script shows how to use our GalSim interface to create FITS images
that contain stars and galaxies
"""

import os
import galsim
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import StarObj, GalaxyBulgeObj, GalaxyDiskObj, GalaxyAgnObj, \
                                                 OpSim3_61DBObject
from lsst.sims.GalSimInterface import SNRdocumentPSF, GalSimStars, GalSimGalaxies, GalSimAgn

#if you want to use the actual LSST camera
#from lsst.obs.lsstSim import LsstSimMapper

def usage():
    print "--------------------------------------------------------"
    print "Usage: "
    print "python2.7 galSimGalaxiesOnly.py -l <lower redshift> -u <upper redshift>"
    print "--------------------------------------------------------"
    return 0;

import getopt
import sys

try:
    opts, args = getopt.getopt(sys.argv[1:], "h:l:u:")

except getopt.GetoptError, err:
    # Print debug info
    print getopt.GetoptError
    print str(err)
    sys.exit(2)

found_p1 = False
found_p2 = False
for opt, arg in opts:
    if opt in ("-h", "--help"):
        usage()
        sys.exit(2)
    elif opt in ("-l", ):
        zdn = arg
        found_p1 = True
    elif opt in ("-u", ):
        zup = arg
        found_p2 = True

if not (found_p1):
    print "*****************************"
    print "Lower redshift limit of the redshift bin was not given."
    usage()
    sys.exit(2)

if not (found_p2):
    print "*****************************"
    print "Upper redshift limit of the redshift bin was not given."
    usage()
    sys.exit(2)

class testGalSimStars(GalSimStars):
    #only draw images in the u and g band for speed
    bandpassNames = ['u','g','r','z']

    #convolve with a PSF; note that galaxies are not convolved with a PSF
    #PSF defined in galSimInterface/galSimUtilities.py
    PSF = SNRdocumentPSF()

    #If you want to use the LSST camera, uncomment the line below.
    #You can similarly assign any camera object you want here
    #camera = LsstSimMapper().camera

class testGalSimGalaxies(GalSimGalaxies):
    bandpassNames = ['u','g','r','z']

    PSF = SNRdocumentPSF()

    #If you want to use the LSST camera, uncomment the line below.
    #You can similarly assign any camera object you want here
    #camera = LsstSimMapper().camera




class testGalSimAgn(GalSimAgn):
    bandpassNames = ['u','g','r','z']

    #defined in galSimInterface/galSimUtilities.py
    PSF = SNRdocumentPSF()

    #If you want to use the LSST camera, uncomment the line below.
    #You can similarly assign any camera object you want here
    #camera = LsstSimMapper().camera


#select an OpSim pointing
obsMD = OpSim3_61DBObject()
obs_metadata = obsMD.getObservationMetaData(88625744, 0.05, makeCircBounds = True)

catName = './outputs_cats/galSim_galaxies_example.txt'


cstrn = "redshift>"+str(zdn)+" and redshift<="+str(zup)

bulges = CatalogDBObject.from_objid('galaxyBulge')
bulge_galSim = testGalSimGalaxies(bulges, obs_metadata=obs_metadata, column_outputs=['redshift'],constraint=cstrn)
bulge_galSim.noise_and_background = None
bulge_galSim.PSF = None
bulge_galSim.write_catalog(catName, write_header=True, write_mode='a')

print 'done with bulges'

disks = CatalogDBObject.from_objid('galaxyDisk')
disk_galSim = testGalSimGalaxies(disks, obs_metadata=obs_metadata, column_outputs=['redshift'],constraint=cstrn)
disk_galSim.copyGalSimInterpreter(bulge_galSim)
disk_galSim.noise_and_background = None
disk_galSim.PSF = None
disk_galSim.write_catalog(catName, write_header=True, write_mode='a')

print 'done with disks'

agn = CatalogDBObject.from_objid('galaxyAgn')
agn_galSim = testGalSimAgn(agn, obs_metadata=obs_metadata, column_outputs=['redshift'],constraint=cstrn)
agn_galSim.copyGalSimInterpreter(disk_galSim)
agn_galSim.noise_and_background = None
agn_galSim.PSF = None
agn_galSim.write_catalog(catName, write_header=True, write_mode='a')

agn_galSim.write_images(nameRoot='./outputs_fits/galaxy_'+str(zdn)+'_'+str(zup))

print 'done with galaxyAgn'
