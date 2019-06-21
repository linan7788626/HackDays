#!/usr/bin/python
import numpy as np
import sys
np.set_printoptions(threshold=np.nan)
import csv
import math
import matplotlib.pyplot as plt
from numpy import zeros,abs,multiply,array,reshape

#---------------------------------------------
#This code finds the distances between stars in a set of multiple exposures
#This is the control panel, where we toggle things on and off

#You will need to make these directories if they do not exist:  'astro_out' and 'atmospheres_new'

start_image=9     							   #The starting number of images to run on
end_image=11 								   #The ending number of images to run on
start_iterate=0								   #The starting number of iterations of an image
end_iterate=10								   #The ending number of iterations of an image
the_directory=''							   #Declares the next value a string
the_directory='examples/'				           #The directory where the input catalogs are
prefix='star_'								   #The prefix of the image name, before the numbers
bin_increment=50							   #Size of steps to make in bins of distance
number_of_columns=(end_iterate-start_iterate)+1
#---------------------------------------------

#LOOP THROUGH ALL SETS OF INPUT IMAGES
for p in range (start_image, end_image):


	#LOOP THROUGH THE 10 EXPOSURES
	things=['far0','far1','far2','far3','far4','far5','far6','far7','far8','far9']

	for z in range (start_iterate, end_iterate):
		name=the_directory+"star_"+str(p)+"_"+str(z)+".cat"
		print "loading ", name
		m = np.loadtxt(name, skiprows=37)

		id1=m[:,0]
		x1=m[:,6]
		y1=m[:,7]

#	a=zip(id1,x1,y1,flag1)

		id3=np.array(id1)
		x3=np.array(x1)
		y3=np.array(y1)

		inds=x3.argsort()

		sorted_id=id3[inds]
		sorted_x=x3[inds]
		sorted_y=y3[inds]

		for d in range (0, len(x1)): #len(x1)): 
		
			x4=sorted_x[d]
			y4=sorted_y[d]
			num1=sorted_id[d]

#		var = raw_input("Please enter something: ")

			dx=(sorted_x-x4)
			dy=(sorted_y-y4)
			dist = (((dx**2)+(dy**2))**0.5)*0.2
			num2=(sorted_id)

			r=zip(dist,num2)

			c = zeros((len(r),2))
		
			for n in xrange(len(r)):
		#	print r[n][1]
				if r[n][1] != num1 and r[n][1] > num1: 
					c[n][1] = 1
					c[n][0] = 1	

			dist_fine=r*c

			if d == 0:
				phil = open(the_directory+"astro_"+str(p)+"_"+str(z)+".out","w")
			else:
				phil = open(the_directory+"astro_"+str(p)+"_"+str(z)+".out","a")
		
			for ab in range (0, len(dist)):
				if dist_fine[ab][1] != 0:
					zink = str(int(num1))+" "+str(int(dist_fine[ab][1]))+" "+str(dist_fine[ab][0])
					print >>phil, zink
		
			phil.close()
#-------------------------------------------------------------------------------------
		q = np.loadtxt(the_directory+"astro_"+str(p)+"_"+str(z)+".out")
		dista = q[:,2] 
		indy = dista.argsort()
		things[z] = dista[indy]

	for h in range (0, len(things[0])):

		if h == 0:
			phile = open(the_directory+"run"+str(p)+".out","w")
		else:
			phile = open(the_directory+"run"+str(p)+".out","a")

		everything = [things[0][h],things[1][h],things[2][h],things[3][h],things[4][h],things[5][h],things[6][h],things[7][h],things[8][h],things[9][h]]
		distancia = np.mean(everything)
		deviation = np.std(everything)*1000.

		print >>phile, distancia, deviation

	phile.close()

	finality = np.loadtxt(the_directory+"run"+str(p)+".out")

	distance=finality[:,0]
	dispersion=finality[:,1]

	bins=[0.,50.,100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,600.,650.,700.,750.,800.,850.,900.,950.,1000.,1050.]

	inds=np.digitize(distance,bins)

	#everything=zip(finality,inds)

	distance_means = [distance[inds == i].mean() for i in range (1, len(bins))]
	dispersion_means = [dispersion[inds == i].mean() for i in range (1, len(bins))]

	#var = raw_input("Please enter something: ")


	#np.histogram(distance)

	plt.plot(distance_means, dispersion_means, linewidth=1.0)
	plt.plot(distance_means, dispersion_means, 'bo')

	plt.xlabel('Distance (arcseconds)')
	plt.ylabel('Astrometric error (mas)')
	plt.axis([0,800,0,50])
	
plt.show()
