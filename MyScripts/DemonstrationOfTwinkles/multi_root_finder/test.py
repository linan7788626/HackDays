"""
Test the GSL cython interface
"""

from nose.tools import assert_almost_equal

import numpy

import gsl_test
#import weave_test

def test_cython_gsl_bessel():
    x = 5.0
    assert_almost_equal(gsl_test.gsl_sf_bessel_jo(x), -1.775967713143382920e-01, 15)

#def test_weave_gsl_bessel():
    #x = 5.0
    #assert_almost_equal(weave_test.weave_gsl_bessel(x), -1.775967713143382920e-01, 15)
