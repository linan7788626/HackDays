from scipy import integrate
import numpy as np
from math import *

#zl = 0.183
#zs = 1.0
h = 0.7
Om0 = 0.27
Ol0 = 0.73
Otot = Om0+Ol0
Ok0 = 1.0-Otot
w = -1
rho_crit0 = 2.78e11 #M_sun Mpc^-3 *h*h
rho_bar0 = rho_crit0*Om0
#sigma8 = 0.801 #wmap 7th
sigma8 = 0.8
#----------------------------------------------------------------------------
apr =  206269.43		#1/1^{''}
vc = 2.9970e5			#km/s
#G = 6.67259e-11			#m^3/kg/s^2
G = 4.3e-9			#(Mpc/h)^1 (Msun/h)^-1 (km/s)^2
H0 = 100.0			#km/s/(Mpc/h)
pc = 3.085677e16		#m
kpc = 3.085677e19		#m
Mpc = 3.085677e22		#m
Msun = 1.98892e30       	#kg
#1Gly = 9.461e26cm = 9.461e24m = 9.461e21km
#1Mpc = 3.08568e24cm = 3.261566ly
dy = 31536000.0/365.0
yr = 31536000.0/365.0*365.25			#second
Gyr = yr*1e9			#second
#-----------------------------------------------------------------------------
#def flat(self):
#	omega_total = self.omega_m0 + self.omega_l0
#	if omega_total <= 1.0001 and omega_total >= 0.9999:
#		return True
#	else:
#		return False
#
#def open(self):
#	omega_total = self.omega_m0 + self.omega_l0
#	if omega_total <= 0.9999:
#		return True
#	else:
#		return False
#
#def closed(self):
#	omega_total = self.omega_m0 + self.omega_l0
#	if omega_total > 1.0001:
#		return True
#	else:
#		return False
#-----------------------------------------------------------------------------
def efunc(x):
	#res = 1.0/np.sqrt(Om0*(1.0+x)**3+Ok0*(1.0+x)**2+Ol0*(1.0+x)**(3*(1.0+w)))
	res = 1.0/np.sqrt(Om0*(1.0+x)**3+Ok0*(1.0+x)**2+Ol0)
	return res
#-----------------------------------------------------------------------------
def tfunc(x):
	res = efunc(x)/(1.0+x)
	return res
#-----------------------------------------------------------------------------
def Hz(x):
	res = H0/efunc(x)
	return res
#-----------------------------------------------------------------------------
def a(x):
	res = 1.0/(1.0+x)
	return res
#-----------------------------------------------------------------------------
def Dh():
	res = vc/H0
	return res
#-----------------------------------------------------------------------------
def Th():
	res = 1.0/H0*1e3*kpc/Gyr
	return res
#-----------------------------------------------------------------------------
def Tl(x):
	res = Th()*integrate.romberg(tfunc, 0, x)
	return res
#-----------------------------------------------------------------------------
def age_z(x):
	res = Th()*integrate.quad(tfunc, x, +np.inf)[0]
	return res
#-----------------------------------------------------------------------------
def Dc(x):
	#res = Dh()*integrate.romberg(efunc, 0, x)
	res = Dh()*integrate.quad(efunc, 0, x)[0]
	return res
#-----------------------------------------------------------------------------
def Dm(x):
	sOk0 = np.sqrt(np.abs(Ok0))
	if Ok0 > 0:
		res = Dh()/sOk0*np.sinh(sOk0*Dc(x)/Dh())
	elif Ok0 == 0.0 :
		res = Dc(x)
	else:
		res = Dh()/sOk0*np.sin(sOk0*Dc(x)/Dh())
	return res
#-----------------------------------------------------------------------------
def Da(x):
	res = Dc(x)/(1.0+x)
	return res
#-----------------------------------------------------------------------------
def Da2(x1, x2):
	dm1 = Dm(x1)
	dm2 = Dm(x2)
	tmp1 = dm2*np.sqrt(1.0+Ok0*dm1**2/Dh()**2)
	tmp2 = dm1*np.sqrt(1.0+Ok0*dm2**2/Dh()**2)

	res = 1.0/(1.0+x2)*(tmp1-tmp2)
	return res
#-----------------------------------------------------------------------------
def Dl(x):
	res = Dc(x)*(1.0+x)
	return res
#-----------------------------------------------------------------------------
def Dp(x1,x2):
	res = Dh()*integrate.romberg(tfunc, x1, x2)
	return res
#-----------------------------------------------------------------------------
def DistMod(x):
	res = 5.0*np.log10(Dl(x)*(Mpc/h)/pc/10.0)
	return res
#-----------------------------------------------------------------------------
def Omz(x):
	return Om0*(1.0+x)**3*efunc(x)**2

def Olz(x):
	return Ol0*(1.0+x)**0*efunc(x)**2

def Okz(x):
	return Ok0*(1.0+x)**2*efunc(x)**2
#-----------------------------------------------------------------------------
#def delta_c(self, redshift=None):
#	"""Over-density threshold for linear, spherical collapse."""
#	# Fitting function taken from NFW97
#	delta_c = 0.15*(12.0*numpy.pi)**(2.0/3.0)
#	if self.open():
#	    delta_c *= self.omega_m(redshift)**0.0185
#	if self.flat() and self.omega_m0 < 1.0001:
#	    delta_c *= self.omega_m(redshift)**0.0055
#
#	return delta_c/self.growth_factor(redshift)
#
#def delta_v(self, redshift=None):
#	"""Over-density for a collapsed, virialized halo."""
#	# Fitting function taken from NFW97
#	delta_v = 178.0
#	if self.open():
#	    delta_v /= self.omega_m(redshift)**0.7
#	if self.flat() and self.omega_m0 < 1.0001:
#	    delta_v /= self.omega_m(redshift)**0.55
#	return delta_v/self.growth_factor(redshift)
#-----------------------------------------------------------------------------
def rho_crit(x):
	"""Critical density in solar masses per cubic Mpc."""
	#rho_crit0 = 1.0e-29*1.0e-33*2.937999e+73
	rho_crit =  rho_crit0/efunc(x)/efunc(x)
	return rho_crit
def rho_bar(x):
	"""Matter density in solar masses per cubic Mpc."""
	return rho_crit(x)*Omz(x)
#-----------------------------------------------------------------------------
def sigma_crit(x1,x2):
	res = (vc*vc/4.0/np.pi/G*Da(x2)/Da(x1)/Da2(x1,x2))
	return res
#-----------------------------------------------------------------------------
#print '---------------------------------------------------------'
#print 'Cosmology : h			=', h
#print 'Cosmology : OmegaM0		=', Om0
#print 'Cosmology : OmegaL0		=', Ol0
#print 'Lensing   : zl			=', zl
#print 'Lensing   : zs			=', zs
#print 'Constant  : arcsec per rad	=', apr
#print 'Constant  : Hubble costant	= %.2e  km/s/(Mpc/h)' %H0
#print 'Constant  : Hubble distance	= %.2e  Mpc/h' %Dh()
#print 'Constant  : Hubble time		= %.2e  Gyr/h' %Th()
#print 'Constant  : G			= %.2e  Mpc/Msun/(km/s)^2' % G
#print 'Constant  : c			= %.2e  km/s' % vc
#print 'Constant  : Msun 		= %.2e  kg' % Msun
#print 'Constant  : pc			= %.2e  m' % pc
#print '---------------------------------------------------------'
#print 'For zl = %.3f:' % (zl)
#print 'Omega_m(z)			= %.2e' % Omz(zl)
#print 'Omega_Lambda(z)			= %.2e' % Olz(zl)
#print 'H(z)				= %.2e km/s/(Mpc/h)' % Hz(zl)
#print 'Comoving L.O.S. distance 	= %.2e Mpc/h' % Dc(zl)
#print 'Angular diameter distance	= %.2e Mpc/h' % Da(zl)
#print 'Luminosity distance		= %.2e Mpc/h' % Dl(zl)
#print 'arcsec per Mpc/h 		= %.2e ' % (1.0/(Dl(zl)/apr))
#print 'Critical Density 		= %.2e (Msun/h)/(Mpc/h)^3' % rho_crit(zl)
#print 'Average Density			= %.2e (Msun/h)/(Mpc/h)^3' % rho_bar(zl)
#print 'Lookback time			= %.2e Gyr/h' % Tl(zl)
#print 'Age of the Universe		= %.2e Gyr/h' % age_z(zl)
#print 'Scale Factor			= %.2e ' % a(zl)
#print 'Distance modulus 		= %.2e ' % DistMod(zl)
#print '---------------------------------------------------------'
#print 'For zs = %.3f:' % (zs)
#print 'Omega_m(z)			= %.2e' % Omz(zs)
#print 'Omega_Lambda(z)			= %.2e' % Olz(zs)
#print 'H(z)				= %.2e km/s/(Mpc/h)' % Hz(zs)
#print 'Comoving L.O.S. distance 	= %.2e Mpc/h' % Dc(zs)
#print 'Angular diameter distance	= %.2e Mpc/h' % Da(zs)
#print 'Luminosity distance		= %.2e Mpc/h' % Dl(zs)
#print 'arcsec per Mpc/h 		= %.2e ' % (1.0/(Dl(zs)/apr))
#print 'Critical Density 		= %.2e (Msun/h)/(Mpc/h)^3' % rho_crit(zs)
#print 'Average Density			= %.2e (Msun/h)/(Mpc/h)^3' % rho_bar(zs)
#print 'Lookback time			= %.2e Gyr/h' % Tl(zs)
#print 'Age of the Universe		= %.2e Gyr/h' % age_z(zs)
#print 'Scale Factor			= %.2e ' % a(zs)
#print 'Distance modulus 		= %.2e ' % DistMod(zs)
#print '---------------------------------------------------------'
#print 'Critical surface density 	= %.2e (Msun/h)/(Mpc/h)^2' % sigma_crit(zl,zs)
#print 'Dals				= %.2e Mpc/h' % (Da2(zl,zs))
#print '---------------------------------------------------------'
