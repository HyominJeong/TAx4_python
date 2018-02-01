from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as plt
import os, sys

# Return position delta from Cartesian coordinate
# pos1, pos2 = [x,y,z] in numpy array
# output is float
# Hyomin Jeong
# 20180117
def posDelta(pos1,pos2):
	return(np.sqrt((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2+(pos1[2]-pos2[2])**2))

# Return [x,y,z] vector from Ra, Dec
# Ra, Dec in radians
# Output is [x,y,z] numpy array
# Hyomin Jeong
# 20180117
def radec2xyz(Ra, Dec):
	x = np.sin(0.5*np.pi-Dec)*np.cos(1.*Ra)
	y = np.sin(0.5*np.pi-Dec)*np.sin(1.*Ra)
	z = np.cos(0.5*np.pi-Dec)
	return(np.array([x,y,z]))

# Return [Ra, Dec] in radian from [X,Y,Z] float vector
# X, Y, Z can be numpy array
# Hyomin Jeong
# 20180130
def xyz2radec(X, Y, Z):
	R = np.sqrt(X**2 + Y**2 + Z**2)
	theta = np.arccos(Z / R)
	phi = np.arctan2(X, Y)
	Ra = phi
	Dec = 0.5 * np.pi - theta
	return Ra, Dec

# Return solid angle of each Ra, Dec bin
# ra_rad, dec_rad, delta_rad in radian unit
# Hyomin Jeong
# 20180130
def dohm(ra_rad, dec_rad, delta_rad):
	theta_rad = 0.5 * np.pi - dec_rad
	phi_rad = ra_rad

	if theta_rad > np.pi:
		print "Strange Dec!", ra_rad, dec_rad, delta_rad

	result = delta_rad * \
			 (-np.cos(theta_rad+0.5*delta_rad) + np.cos(theta_rad-0.5*delta_rad))

	if result < 0:
		print "Strange Result!", ra_rad, dec_rad, delta_rad, result
	return result

# Return homogeneous Ra, Dec sets
# nOfEvts : integer
# Output is [[Ra sets], [Dec sets]] numpy array in radian unit
# Hyomin Jeong
# 20180131
def CRHradec(nOfRand, show=0):
	# Ra from 0 to 360
	ra_random = np.random.rand(nOfRand)*2*np.pi#*360.0*u.degree

	# Dec from -90 to 90
	#dec_random=(np.random.rand(nOfRand)*2 -1)*2*np.pi#*180.0-90.0)*u.degree
	theta_random=np.arccos(np.random.rand(nOfRand)*2-1)
	dec_random=(0.5 * np.pi - theta_random)# * 180 / np.pi *u.degree

	###########################################################
	# Show polar scatter plot of Ra, Dec = (0, 0) and (0, 90) #
	###########################################################
	if show == 1:
		fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, subplot_kw=dict(projection='polar'))

		# Vernal equinox plot
		theta = np.arctan2(np.sin(dec_random), np.cos(dec_random)*np.sin(ra_random))
		r = np.sqrt((np.sin(ra_random)**2 + (np.sin(dec_random)**2*(np.cos(ra_random)**2))))

		ax1.scatter(theta.flatten(), r.flatten(), s=0.5, alpha= 0.3, marker='.', color='black')
		ax1.set_title("Center: Vernal equinox")

		# Polar plot
		theta=ra_random
		r = np.cos(dec_random)

		ax2.scatter(theta.flatten(), r.flatten(), s=0.5, alpha= 0.3, marker='.', color='black')
		ax2.set_title("Center: North Pole")

		# Y axis plot
		theta = np.arctan2(np.sin(dec_random),- np.cos(dec_random)*np.sin(ra_random))
		r = np.sqrt((np.sin(ra_random)**2 + (np.sin(dec_random)**2*(np.cos(ra_random)**2))))

		ax3.scatter(theta.flatten(), r.flatten(), s=0.5, alpha= 0.3, marker='.', color='black')
		ax3.set_title("Center: Y axis")


		plt.show()

	return ra_random, dec_random

# Return position delta histogram from Ra, Dec sets
# Ra, Dec sets are array in radian unit
# Output is numpy histogram of position delta in radian unit
# Hyomin Jeong
# 20180131
def histPosDelta(Ras, Decs, nbin):
	print "Getting position deltas"
	print "Total number of events is", len(Ras)
	rad_deltas = []
	for n in range(len(Ras)):
		for i in range(n):
			pos1 = radec2xyz(float(Ras[i]),float(Decs[i]))
			pos2 = radec2xyz(float(Ras[n]),float(Decs[n]))
			#print pos1, pos2, i, n
			len_delta = posDelta(pos1, pos2)
			#print 2.*np.arcsin(0.5*len_delta)
			rad_deltas.append(2.*np.arcsin(0.5*len_delta))

	# Draw histogram of 
	hist_rad_deltas, binEdge_rad_deltas = np.histogram(rad_deltas, bins=nbin)
	return hist_rad_deltas, binEdge_rad_deltas