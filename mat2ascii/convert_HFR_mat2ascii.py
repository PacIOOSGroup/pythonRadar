
import scipy.io
import math
import numpy as np
from numpy import array
from datetime import datetime
import sys


def main():
	f = scipy.io.loadmat('radar_in.mat')
	numVars = len(scipy.io.whosmat('radar_in.mat'))

	try:
		year,month,day,hours,minutes,seconds,place = tuple(sys.argv[1:])

	except ValueError:
		print("Not enough command line argumentsarguments \n"
		+ "Format is year, month, day, hour, minute, second, location")
		quit()
#	print(scipy.io.whosmat('radar_in.mat'))

	print("Reading Matlab binary file...")
	#get the values from the table
	slon=f['slon']
	slat=f['slat']
	X = f['X']
	Y = f['Y']
	VELO = f['Ur']
	theta = f['theta']
	r = f['r']
	EACC = f['Acc']
	MANU = f['Manufacturer']
	SITE = f['Site']
	RRKM=f['RangeResolutionKMeters']
	TCFM=f['TransmitCenterFreqMHz']
	DRHB=f['DopplerResolutionHzPerBin']
	CVLM=f['CurrentVelocityLimit']
	TSRH=f['TransmitSweepRateHz']
	TBKH=f['TransmitBandwidthKHz']
	SDCS=f['SpectraDopplerCells']

	im, jm = np.shape(X)
	#flatten the 2d arrays to 1d so they're easier to work with
	X = X.flatten(order='F')
	Y = Y.flatten(order='F')
	r = r.flatten(order='F')
	
	#initialize arrays so python doesn't complain
	BEAR = []
	RNGE = []
	XLON = []
	YLAT = []			
	XDST = []
	YDST = []
	VELU = []
	VELV = []
	VEL = []
	ACC = []
	HEAD = []
	VELO=VELO.flatten(order='F')
	EACC=EACC.flatten(order='F')
	#print(np.shape(X))
	iomax = 0
	ij = 0
	print("Performing calculations...")
	for j in range(0, jm):
		for i in range (0, im):
			xc1 = X[ij]
			xc2 = Y[ij]
			xc3 = theta[0,j]
			xc4 = r[i]
			xc5 = VELO[ij]
			xc6 = EACC[ij]	
			#don't want to print any NaNs
			if not (math.isnan(xc1) or math.isnan(xc2) or math.isnan(xc3)
			or math.isnan(xc4) or math.isnan(xc5) or math.isnan(xc6)):
				BEAR.insert(iomax, theta[0,j])
				RNGE.insert(iomax, r[i])
				XLON.insert(iomax, X[ij])
				YLAT.insert(iomax, Y[ij])

				if(r[i] >= 0):
					XDST.insert(iomax, r[i]*np.sin(theta[0,j] * np.pi / 180.0))
					YDST.insert(iomax, r[i]*np.cos(theta[0,j] * np.pi / 180.0))
					VELU.insert(iomax, VELO[ij] * np.sin(theta[0,j] * np.pi / 180) * 100)
					VELV.insert(iomax, VELO[ij] * np.cos(theta[0,j] * np.pi / 180) * 100)
					VEL.insert(iomax,-1.0 * VELO[ij] * 100)
					ACC.insert(iomax, (EACC[i]*100))
					rad_bear = theta[0,j] * np.pi / 180.0
					rad_lat = slat * np.pi / 180.0
					rad_lon = slon * np.pi / 180.0
					rang = 1000.0 * r[i]
					
					rad_head = direct(rad_lat, rad_lon, rad_bear, rang)
					HEAD.insert(iomax, rad_head * 180.0 / np.pi)
					iomax += 1	
			ij = ij+1

	f = open('RDL_' + place + '_' + year + '_' + month + '_' + day + '_' + hours + minutes + seconds + '.ruv', "wt")
	print("Writing to " + str(f.name))
	
	#print header data
	f.write("%CTF: 1.00\n")
	f.write("%FileType: LLUV rdls\n")
	f.write("%Manufacturer: " + str(MANU[0]) + "\n")
	f.write("%LLUVSpec: 1.00 2007 12 06\n")
	f.write("%Site: " + str(SITE[0]) + "\n")
	f.write("%TimeStamp: " + year +" "+ month + " " + day + " " + hours + " " + minutes + " " + seconds + "\n")
	#current datetime= datetime.now().strftime('%Y-%m-%d %H:%M:%S' + "\n"))
	f.write("%TimeZone: UTC + 0.00 0 \n")
	f.write("%Origin: " + str(slat[0]) + "," + str(slon[0]) + "\n")
	f.write("%GreatCircle: \"WGS84\" 6378137.000 298.257223562997\n")
	f.write("%GeodVersion: \"Vincenty (1979)\" 2.0 2002 10 01\n")
	f.write("%RangeResolutionKMeters: " + str(RRKM[0]) + "\n")
	f.write("%TransmitCenterFreqHMz: " + str(TCFM[0]) + "\n")
	f.write("%DopplerResolutionHzPerBin: " + str(DRHB[0]) + "\n")
	f.write("%CurrentVelocityLimit: " + str(CVLM[0]) + "\n")
	f.write("%MergedCount: 0001\n")
	f.write("%TransmitSweepRateHz: " + str(TSRH[0]) + "\n")
	f.write("%TransmitBandwidthKhz: " + str(TBKH[0]) + "\n")
	f.write("%SpectraDopplerCells: " + str(SDCS[0]) + "\n")
	f.write("%TableType: LLUV RDL1\n")
	f.write("%TableColumns: 11\n")
	f.write("%TableColumnTypes: LOND LATD VELU VELV EACC XDST YDST YRNGE BEAR VELO HEAD\n")
	f.write("%TableRows: " + str(len(XDST)) + "\n")
	f.write("%TableStart:\n")
	f.write("%% %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s" % ("Longitude","Latitude","U comp","V comp","Accuracy","X Distance",
		 "Y Distance","Range","Bearing","Velocity","Direction\n"))
	f.write("%% %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s" % ("(deg)", "(deg)", "(cm/s)", "(cm/s)", "(cm/s)", "(km)", "(km)", "(km)", "(deg NCW)", "(cm/s)",
			"(deg NCW)\n"))
	#print the actual values
	for i in range(len(XDST)):
		f.write('{0:15.7f} {1:15.7f} {2:15.3f} {3:15.3f} {4:15.3f} {5:15.4f} {6:15.4f} {7:15.3f} {8:15.0f} {9:15.3f} {10:15.4f}'.format(XLON[i], YLAT[i], VELU[i], VELV[i], ACC[i], XDST[i], YDST[i], RNGE[i], BEAR[i], VEL[i], HEAD[i]))
		f.write("\n")	
	f.write("%TableEnd:\n")
	f.write("%End:")
	f.close()
					
# fancy geomancy!
# solution of the geodetic direct problem after t.vincenty
# modified rainsford's method with helmert's elliptical terms
#  effective in any azimuth and at any distance short of antipodal
#  semiMajorAxis is the semi-major axis of the reference ellipsoid
#  flattening is the flattening of the reference ellipsoid
#  latitudes and longitudes in radians positive north and east
#  azimuths in radians clockwise from north
#  geodesic distance s assumed in units of semi-major axis a
# programmed for cdc-6600 by lcdr l.pfeifer ngs rockville md 20feb75
# modified for system 360 by john g gergen ngs rockville md 750608
# translated to python by owen g wilcox uh-soest manoa hi 180913
def direct(glat1, glon1, bear, rang):
	eps = 0.5-13
	rad = 180.0 / np.pi
	semiMajorAxis = 6378137.0
	flattening = 1.0 / 298.257223563

	r = 1.0 - flattening
	tu = r * np.sin(glat1) / np.cos(glat1)
	sf = np.sin(bear)
	cf = np.cos(bear)
	head = 0.0
	if(cf != 0.0):
		head = np.arctan2(tu, cf) * 2.0
	cu = 1.0 / np.sqrt(tu * tu + 1.0)
	su = tu * cu
	sa = cu * sf
	c2a = -sa * sa + 1.0
	x = np.sqrt((1.0/r/r - 1.0) * c2a + 1.0) + 1.0
	x = (x - 2.0) / x
	c = 1.0 - x
	c = (x*x / 4.0 + 1 ) / c
	d = (0.375 * x * x - 1.0) * x
	tu = rang / r / semiMajorAxis / c
	y = tu
	iflag = 1
	while True:
		sy = np.sin(y)
		cy = np.cos(y)
		cz = np.cos(head+y)
		e = cz * cz * 2.0 - 1.0
		c = y
		x = e * cy
		y = e + e - 1.0
		y = ( ( ( sy * sy * 4.0 - 3.0 ) * y * cz * d / 6.0 + x ) * d / 4.0 - cz ) * sy * d + tu
		iflag = iflag + 1
		if not (np.abs(y-c) > eps and iflag < 20):
			break
	head = cu * cy * cf - su * sy
	c = r * np.sqrt(sa*sa + head*head)
	d = su * cy * cf - su * sy
	glat2 = np.arctan2(sy * sf, c)
	c = cu * cy + cu * sy * cf
	x = np.arctan2(sy*sf, c)
	c = ( ( -3.0 * c2a + 4.0 ) * flattening + 4.0 ) * c2a * flattening / 16.0
	d = ( ( e * cy * c + cz ) * sy * c + y ) * sa
	glon2 = glon1 + x - (1.0 - c) * d * flattening
	head = np.arctan2(sa, head) + np.pi
	return head[0,0]

if __name__ == "__main__":
	main()