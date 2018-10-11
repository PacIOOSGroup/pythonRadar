from netCDF4 import Dataset
import scipy.io
import time
import datetime
import numpy as np
from numpy import array
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

	#get the values from the table
	slon=f['slon']
	slat=f['slat']
	X = f['X']
	Y = f['Y']
	VELO = f['Ur']
	theta = f['theta']
	theta = theta.astype(float) #if this is left as an unsigned 2byte integer, can't set missing to -999.0 later 
	r = f['r']
	EACC = f['Acc']
	POW=f['Pow']
	varian=f['Var']

	imax = 245
	jmax = 125
	xmiss = -999.0
	im, jm = np.shape(X)
	if(place == "kok"):
		title = "Results from Koko Head HFR, raw 15-minute data provided by SOEST RADLAB group, P. Flament, PI." \
		" Data are radial velocities\n x-origin = -157.70330\n y-origin =  21.26139"
	elif(place == "kak"):
		title = "Results from Kakaako HFR, raw 15-minute data provided by SOEST RADLAB group, P. Flament, PI." \
		" Data are radial velocities\n x-origin = -157.85830\n y-origin =  21.29170"
	elif(place == "kal"):
		title = "Results from Kalaeloa HFR, raw 15-minute data provided by SOEST RADLAB group, P. Flament, PI." \
		" Data are radial velocities\n x-origin = -157.08360\n y-origin =  21.29750"
	elif(place == "kna"):
		title = "Results from Kaena Pt HFR, raw 15-minute data provided by SOEST RADLAB group, P. Flament, PI." \
		" Data are radial velocities\n x-origin = -158.2897\n y-origin =  21.57500"
	elif(place == "kkh" ):
		title = "Results from Keaukaha HFR, raw 15-minute data provided by SOEST RADLAB group, P. Flament, PI." \
		" Data are radial velocities\n x-origin = -155.04830\n y-origin =  19.73200"
	elif(place == "ppk" ):
		title = "Results from Pepeekeo HFR, raw 15-minute data provided by SOEST RADLAB group, P. Flament, PI." \
		" Data are radial velocities\n x-origin = -155.08280\n y-origin =  19.84720"
	elif(place == "kap" ):
		title = "Results from Kapolei HFR, raw 15-minute data provided by SOEST RADLAB group, P. Flament, PI." \
		" Data are radial velocities\n x-origin = -158.1169\n y-origin =  21.3114"
	else:
		print("You've entered an invalid location code, or this program hasn't been updated to process every radar location.")
		sys.exit()

	dataset = Dataset('RDL_' + place + '_' + year + '_' + month + day + '_' + hours + minutes + '.nc' , "w", format="NETCDF4")
	print(dataset.file_format)
	
	#check for NaN values and change them to the NaN value that netcdf recognizes
	for i in range(im):
		for j in range(jm):
			if np.isnan(X[i][j]):
				X[i][j] = xmiss
			if np.isnan(Y[i][j]):
				Y[i][j] = xmiss	
			if np.isnan(VELO[i][j]):
				VELO[i][j] = xmiss
			if np.isnan(EACC[i][j]):
				EACC[i][j] = xmiss
			if np.isnan(POW[i][j]):
				POW[i][j] = xmiss
			if np.isnan(varian[i][j]):
				varian[i][j] = xmiss
	
	print("Defining Dimensions")

	dataset.createDimension("range", jmax)
	dataset.createDimension("bearing", imax)
	dataset.createDimension("z", 1)
	dataset.createDimension("time", None)	

	print("Defining Variables")
	#0d
	xorig = dataset.createVariable("x-origin", "f4")
	xorig.long_name="longitude of HFR at zero range"
	xorig.short_name="slon"
	xorig.units="degrees_east"

	yorig = dataset.createVariable("y-origin", "f4")
	yorig.long_name="latitude of HFR at zero range"
	yorig.short_name="slat"
	yorig.units="degrees_north"

	#1d
	rangevar = dataset.createVariable("range", "f4", ("range",),fill_value=xmiss)
	rangevar.long_name="range to cell"
	rangevar.units="km"

	bearingvar = dataset.createVariable("bearing", "f4", ("bearing",),fill_value=xmiss)
	bearingvar.long_name="bearing to cell clockwise from north"
	bearingvar.units="degrees"

	zvar = dataset.createVariable("z", "f4", ("z",))
	zvar.long_name="depth below mean sea level"
	zvar.standard_name="depth"
	zvar.axis = "z"
	zvar.units="meters"

	timevar = dataset.createVariable("time", "f4", ("time",))
	timevar.long_name = "time"
	timevar.standard_name="time"
	timevar.axis="t"
	timevar.units="minutes since 2008-01-01 00:00:00"
	timevar.calendar="gregorian"

	#2d
	latvar = dataset.createVariable("lat", "f4", ("range", "bearing",),fill_value=xmiss)
	latvar.long_name = "latitude"
	latvar.standard_name="latitude"
	latvar.axis="y"
	latvar.units="degrees_north"
	

	lonvar = dataset.createVariable("lon", "f4", ("range", "bearing",),fill_value=xmiss)
	lonvar.long_name="longitude"
	lonvar.standard_name="longitude"
	lonvar.axis="x"
	lonvar.units = "degrees_east"
	

	#4d
	accuvar = dataset.createVariable("accu", "f4", ("time", "z", "range", "bearing",),fill_value=xmiss)
	accuvar.long_name="accuracy of the Bragg estimate in cell"
	accuvar.units="m/s"

	powrvar = dataset.createVariable("powr", "f4", ("time", "z", "range", "bearing",),fill_value=xmiss)
	powrvar.long_name="relative power of cell"
	powrvar.units="dB"

	uradvar = dataset.createVariable("urad", "f4", ("time", "z", "range", "bearing",),fill_value=xmiss)
	uradvar.long_name="radial current in cell"
	uradvar.standard_name="radial_sea_water_velocity_away_from_instrument"
	uradvar.units="m/s"

	varivar = dataset.createVariable("vari", "f4", ("time", "z", "range", "bearing",),fill_value=xmiss)
	varivar.long_name="variance of grid interpolation"
	varivar.units="unknown"

	dataset.title = title
	shapeX, shapeY = np.shape(Y)

	# time since 2008 in minutes

	currentDate = datetime.datetime.now()
	firstOf2008 = datetime.datetime.strptime('2008-01-01 00:00:00', '%Y-%m-%d %H:%M:%S')
	#thanks, SO: https://stackoverflow.com/questions/2788871/date-difference-in-minutes-in-python
	currentDateSeconds = time.mktime(currentDate.timetuple())
	firstOf2008Seconds = time.mktime(firstOf2008.timetuple())

	minutesSince = (currentDateSeconds - firstOf2008Seconds) / 60

	theta = np.reshape(theta, 121)

	Y = np.pad(Y, [(0,125-shapeX),(0,245-shapeY)], mode='constant', constant_values=-999.0)
	X = np.pad(X, [(0,125-shapeX),(0,245-shapeY)], mode='constant', constant_values=-999.0)
	r = np.pad(r, [(0,125-101),(0,0)], mode='constant', constant_values=-999.0)	
	EACC = np.pad(EACC,[(0,125-shapeX),(0,245-shapeY)], mode='constant', constant_values=-999.0)
	VELO = np.pad(VELO, [(0,125-shapeX),(0,245-shapeY)], mode='constant', constant_values=-999.0)
	POW = np.pad(POW, [(0,125-shapeX),(0,245-shapeY)], mode='constant', constant_values=-999.0)
	varian = np.pad(varian, [(0,125-shapeX),(0,245-shapeY)], mode='constant', constant_values=-999.0)
	theta = np.pad(theta,[0,245-121], mode='constant', constant_values=-999.0)

	print("Writing data to output")

	latvar[:,:] = Y
	lonvar[:,:] = X
	rangevar[:] = r
	bearingvar[:] = theta
	timevar[:] = minutesSince
	zvar[:] = 0.0
	xorig[:] = slon
	yorig[:] = slat

	accuvar[:,0,0,0] = minutesSince
	accuvar[0,:,0,0] = 0.0
	accuvar[0,0,:,:] = EACC

	powrvar[:,0,0,0] = minutesSince
	powrvar[0,:,0,0] = 0.0
	powrvar[0,0,:,:] = POW

	uradvar[:,0,0,0] = minutesSince
	uradvar[0,:,0,0] = 0.0
	uradvar[0,0,:,:] = VELO

	varivar[:,0,0,0] = minutesSince
	varivar[0,:,0,0] = 0.0
	varivar[0,0,:,:] = varian

if __name__ == "__main__":
	main()