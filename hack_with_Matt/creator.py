#!/usr/bin/python
import numpy as np
import sys
np.set_printoptions(threshold=np.nan)

#THIS PROGRAM GENERATES THE INSTANCE FILE FOR PHOSIM

#-----------------------------------------
#ENTER INPUT PARAMETERS HERE
filter_num="0"  #0=u, 1=g, 2=r, 3=i, 4=z, 5=y
name='input'  #name of the file with the lines of position, mag, semi-major and minor axes, angle, Sersic index and redshift
outname="instance.file"
number_of_galaxies=1  #<--How many lines are there in the input file?
SED_file="galaxySED/Exp.80E07.0005Z.spec.gz"  #<----write the SED file you want to use here
#-----------------------------------------

name='input'
m = np.loadtxt(name)

xc1=m[:,0]
xc2=m[:,1]
mag_tot=m[:,2]
Reff_a=m[:,3]
Reff_b=m[:,4]
Theta=m[:,5]
ndex=m[:,6]
zl=m[:,7]

line3="Unrefracted_Azimuth -0.021722"
line4="Unrefracted_Altitude 84.420088"
line5="Slalib_date 2000/1/2.389366"
line6="Opsim_rotskypos 0.000000"
line7="Opsim_obshistid 1099"
line8="Opsim_rottelpos 0.020649"
line9="Opsim_moondec -90"
line10="Opsim_rawseeing 0.7"
line11="Opsim_moonra 180"
line12="Opsim_expmjd 51545.389366"
line13="Opsim_moonalt -90"
line14="Opsim_sunalt -90"
line15="Opsim_filter "+filter_num
line16="Opsim_dist2moon 180.0"
line17="Opsim_moonphase 10.0"
line18="SIM_SEED 2222222"
line19="SIM_MINSOURCE 1"
line20="SIM_TELCONFIG 0"
line21="SIM_CAMCONFIG 1"
line22="SIM_NSNAP 1"
line23="SIM_VISTIME 16.5"
line24="SIM_DOMEINT 18.0"
line25="SIM_DOMEWAV 0"
line26="SIM_TEMPERATURE 20.0"
line27="SIM_TEMPVAR 0.0"
line28="SIM_PRESSURE 520.0"
line29="SIM_PRESSVAR 0.0"
line30="SIM_OVERDEPBIAS -45.0"
line31="SIM_CCDTEMP 173.0"
line32="SIM_ALTVAR 0.0"
line33="SIM_CONTROL 0"
line34="SIM_ACTUATOR 0.0"
line35="isDithered 0"
line36="ditherRaOffset 0.0"
line37="ditherDecOffset 0.0"

phile = open(outname,"w")

mean1=np.mean(xc1)
mean2=np.mean(xc2)
print >>phile, 'Unrefracted_RA_deg ', mean1
print >>phile, 'Unrefracted_Dec_deg ', mean2
print >>phile, line3
print >>phile, line4
print >>phile, line5
print >>phile, line6
print >>phile, line7
print >>phile, line8
print >>phile, line9
print >>phile, line10
print >>phile, line11
print >>phile, line12
print >>phile, line13
print >>phile, line14
print >>phile, line15
print >>phile, line16
print >>phile, line17
print >>phile, line18
print >>phile, line19
print >>phile, line20
print >>phile, line21
print >>phile, line22
print >>phile, line23
print >>phile, line24
print >>phile, line25
print >>phile, line26
print >>phile, line27
print >>phile, line28
print >>phile, line29
print >>phile, line30
print >>phile, line31
print >>phile, line32
print >>phile, line33
print >>phile, line34
print >>phile, line35
print >>phile, line36
print >>phile, line37

for h in range (0, number_of_galaxies):
	all_of_it="object "+str(h+1)+" "+str(xc1[h])+" "+str(xc2[h])+" "+\
	str(mag_tot[h])+" "+SED_file+" "+str(zl[h])+\
	" 0 0 0 0 0 sersic2D "+str(Reff_a[h])+" "+str(Reff_b[h])+" "+str(Theta[h])+" "+str(ndex[h])+" none none"
	print >>phile, all_of_it 

phile.close()