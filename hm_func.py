import numpy as np

#def dohm(theta_rad, phi_rad, delta_rad):
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

#print dohm(85.,5.)
#dohms = [dohm(5+10*n, 5) for n in range(9)]
#print sum(dohms) * 36
