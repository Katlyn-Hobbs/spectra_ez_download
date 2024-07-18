################################
#                              #
# CORRECTING DRS VELOCITIES    #
# Adapted from ACC by AM       #
#                              #
#    v1 - 31 Aug 2018          #
# Use JPL files in km and s    #
# Use harpn_extract and mean   # 
# Use 2006 procedure for       #
#  calculating GMST            #
# RA interpolation cyclic      #
#    v2 - Oct 2018             #
# updated FWHM/C correction    #
# updated extinction           #
#  correction                  #
# Need Qualresults.csv file    #
# Daily mean data removed      #
#    v3 - Mar 2019             #
# updated extinction           #
#   correction                 #
# updated faulty timings       #
# Added the binned data        #
#    v4 - Sep 2020             #
# New FWHM correction          #
# factor to be applied to      #
# the new data                 #
#    v5 - June 2021            #
# Now reads output of          #
# GetExtinctionCoefficients.py #
# Applying own extinction      #
# coefficient removing the     #
# need for Qualresults file.   #
#                              #
################################

import numpy as np
import matplotlib.pyplot as plt
import os

############################
# Define tunable constants #
############################

download = 0 # Turn to 1 if JPL files need to be downloaded (using wget in terminal)
             # Note these are different than the ones ACC's Mathematica uses

madclip = 500000.0 # A MAD clip is performed after 'the SS is removed'.
              # Turn to a ridiculous high number (e.g. 5000000000000000000)
              # if you do not want MADclipping

protsun = 25.38 # Solar rotation period in days

fact = 0.0  # Tuning parameter for scaling the contribution of vsini to overall FWHM (between 1.0 and 1.2 usually)
            # Note this is DRS version dependent!!
            # Leave it to 0 to define the best one based on lowest MAD

vrotsun = 2000.0 # Rotational velocity of the Sun in m/s

uu = 0.6 # limb darkening linear coefficient

##################################################
# Define constants                               #
# Source for solar/Earth values: NASA fact sheet #
##################################################

AUinkm = 149597870.7  # AU in kilometers
dins = 24.0*60.0*60.0 # days in seconds
clight = 299792.458   # Speed of light in km/s
tngwlon = (17.0 + 53.0/60. + 20.6/3600.0)*np.pi/180.0 # Western longitude of TNG in radians
tnglat = (28.0 + 45.0/60. + 14.4/3600.0)*np.pi/180.0 # TNG latitude in radians
sunpolera = 286.13*np.pi/180.0 # Solar rotation pole RA, in radians
sunpoledec = 63.87*np.pi/180.0 # Solar rotation pole DEC, in radians
eps = 23.439291111*np.pi/180.0 # Ecliptic obliquity, in radians
ecc = 0.0167 # eccentricity of Earth orbit
yrins = 365.256*dins # Year (earth orbital period) in seconds
rsun = 696000.0 # solar radius in km
protsun = protsun*dins # Solar rotation period in s
ASU = 2.9722*10**(-6) # Snodgrass & Ulrich 1990 sidereal differential rotation law
BSU = -0.484*10**(-6) # Snodgrass & Ulrich 1990 sidereal differential rotation law
CSU = -0.361*10**(-6) # Snodgrass & Ulrich 1990 sidereal differential rotation law

#####################
# Load all H-N data #
#####################

datatable = "Sun_harpn_2-2-3_alldata_quality_extinct.csv"

data = np.loadtxt(datatable,delimiter=',',skiprows=1,dtype='str')

time,rv,err,fwhm,efwhm,contrast,econtrast,airmass,kext = np.loadtxt(datatable,delimiter=',',skiprows=1,unpack=True,  usecols=(1,2,3,4,5,6,7,13,22))

######################
# Option to download # 
# new JPL files      #
######################

if download:

	## Sun-tng

	os.system("echo \"https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND='10'&MAKE_EPHEM='YES'&CENTER='Z19'&TABLE_TYPE='VECTORS'&OUT_UNITS='KM-S'&START_TIME='2015-07-01'&STOP_TIME='2022-09-27'&STEP_SIZE='01%20d'&CSV_FORMAT='YES'\" > tmpi")
	os.system('wget -i tmpi -O sun-tng-kms.csv')
	os.system('rm tmpi')

	## SSB-Sun
	## Needs to be same start and stop time and stepsize as Sun-tng

	os.system("echo \"https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND='0'&MAKE_EPHEM='YES'&CENTER='@sun'&TABLE_TYPE='VECTORS'&OUT_UNITS='KM-S'&START_TIME='2015-07-01'&STOP_TIME='2022-09-27'&STEP_SIZE='01%20d'&CSV_FORMAT='YES'\" > tmpi")
	os.system('wget -i tmpi -O ssb-sun-kms.csv')
	os.system('rm tmpi')

	## Sun-tng-obs

	os.system("echo \"https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND='10'&MAKE_EPHEM='YES'&CENTER='Z19'&TABLE_TYPE='OBSERVER'&START_TIME='2015-07-01'&STOP_TIME='2022-09-27'&STEP_SIZE='01%20d'&CSV_FORMAT='YES'&QUANTITIES='1,13,17,30'&CAL_FORMAT='JD'&ANG_FORMAT='DEG'\" > tmpi")
	os.system('wget -i tmpi -O sun-tng-obs.csv')
	os.system('rm tmpi')

#########################
# Load JPL sun-tng data #
# and scale values      #
#########################

jdb, xsun, ysun, zsun, sundist = np.genfromtxt('sun-tng-kms.csv',skip_header = 50, skip_footer = 50, delimiter=',', unpack=True, usecols=(0,2,3,4,9))
xsun = -xsun/sundist
ysun = -ysun/sundist
zsun = -zsun/sundist
tdelay = sundist/clight/dins # light travel time delay Sun-TNG in days

#########################
# Load JPL ssb-sun data #
#########################

vxbary, vybary, vzbary = np.genfromtxt('ssb-sun-kms.csv',skip_header = 32, skip_footer = 50, delimiter=',', unpack=True, usecols=(5,6,7))

#################################
# Load JPL sun-tng-obs data     #
# and convert values to radians #
# and days                      #
#################################

jd, ra, dec, angdia, nppa, deltat = np.genfromtxt('sun-tng-obs.csv',skip_header = 59, skip_footer = 105, delimiter=',', unpack=True, usecols=(0,3,4,5,6,8))
ra *= np.pi/180.0 # Sun RA
dec *= np.pi/180.0 # Sun DEC
angdia *= np.pi/180.0/3600. # angular diameter
nppa *= np.pi/180.0 # North Pole position angle
deltat /= dins # Difference between TT and UT1

##############################
# Correct DRS velocities for #
# solar barycentric velocity #
# away from TNG in m/s       #
##############################

rvbary = (xsun*vxbary + ysun*vybary + zsun*vzbary) # solar barycentric velocity
rvbary_corr = np.interp(time,jdb,rvbary) # solar barycentric velocity at H-N times

rv_notss = rv - rvbary_corr # Heliocentric HN radial velocity

################################
# Remove outliers using median #
# and mad (option madclip).    #
# Add SS back in afterwards to #
# have outlier-free original   #
# data as well                 #
################################

mrvc = np.median(rv_notss)
sdrvc = np.std(rv_notss)
madrvc = np.median(np.abs(mrvc-rv_notss))
cond = ( np.abs(rv_notss - mrvc) < madclip * madrvc )

data = data[cond]
time,rv,rvbary_corr, rv_notss,err,fwhm,efwhm,contrast,econtrast,airmass,kext = time[cond],rv[cond],rvbary_corr[cond], rv_notss[cond],err[cond],fwhm[cond],efwhm[cond], contrast[cond],econtrast[cond],airmass[cond],kext[cond]

#############################
# Approximate correction to #
# FWHM for solar ecliptic   #
# obliquity and Earth orbit #
# eccentricity.             #
#############################

# Calculate Solar rotation pole in ecliptic coordinates

sinsunpolelat = np.sin(sunpoledec)*np.cos(eps) - np.cos(sunpoledec)*np.sin(eps)*np.sin(sunpolera)
cossunpolelat = np.sqrt(1.-sinsunpolelat*sinsunpolelat)

sinsunpolelong = (np.sin(sunpoledec)*np.sin(eps) + np.cos(sunpoledec)*np.cos(eps)*np.sin(sunpolera))/cossunpolelat
cossunpolelong = np.cos(sunpoledec)*np.cos(sunpolera)/cossunpolelat

# Transform to Cartesian ecliptic coordinates

xsunpole = cossunpolelat*cossunpolelong
ysunpole = cossunpolelat*sinsunpolelong
zsunpole = sinsunpolelat

# Get inclination angle between Sun centre and pole (via dot product)

cosi = xsun*xsunpole + ysun*ysunpole + zsun*zsunpole
sini = np.sqrt(1.-cosi*cosi)

# Calculate Earth angular velocity wearth = dtheta/dt
# Using Kepler second law: pi a b / P = dA / dt = 0.5 r^2 dtheta/dt

wearth = 2.*np.pi*AUinkm*AUinkm*np.sqrt(1.-ecc*ecc)/yrins/sundist/sundist

# interpolate sini and wearth to observation times

siniobs = np.interp(time,jdb,sini)
cosiobs = np.interp(time,jdb,cosi)
wearthobs = np.interp(time,jdb,wearth)

# correct solar angular velocity for earth angular velocity

wsun = 2.*np.pi/protsun
wsunobs = wsun - wearthobs

# Observed and sidereal vsini for the Sun in km/s

vsiniobs = rsun*wsunobs*siniobs
vsinisun = rsun*wsun

# Define the correction factor
if fact == 0.0:

	factors = np.arange(0.90,1.30,0.001)
	rmsfwhm = np.zeros(len(factors))

	def mad(x):
		med = np.median(x)
		mad = np.median(np.abs(x-med))
		return mad

	for ind,fact in enumerate(factors):
		fwhmcor = np.sqrt(fwhm*fwhm - fact*fact*vsiniobs*vsiniobs + fact*fact*vsinisun*vsinisun)
		rmsfwhm[ind] = mad(fwhmcor)

	p = np.polyfit(factors,rmsfwhm,2) # 2-d polynomial fit
	fact = -p[1]/2./p[0] # minimum of that fit

# Correct FWHM with observed vsini

fwhmcor = np.sqrt(fwhm*fwhm - fact*fact*vsiniobs*vsiniobs + fact*fact*vsinisun*vsinisun)

############################
# Correct contrast so that #
# EW is conserved          #
############################

contrastcor = contrast*fwhm/fwhmcor

###################################
# Perform differential extinction #
# and solar parallactic angle     #
# correction procedure            #
###################################

# Correct times for light travel time
# Not strictly correct as we are adding travel time to Sun centre rather
# than the barycentre, but discrepancy only up to 2 or 3 seconds.

timejd = time + np.interp(time, jdb, tdelay)

# Get ra, dec, angular diameter, north pole phase angle, and
# deltat at the right times

ra2 = np.interp(timejd, jd, np.unwrap(ra))%(2.*np.pi)
dec2 = np.interp(timejd, jd, dec)
angdia2 = np.interp(timejd, jd, angdia)
nppa2 = np.interp(timejd, jd, nppa)
deltat2 = np.interp(timejd, jd, deltat)

# GMST procedure adapted from ERFA/gmst06
# based on SOFA version "20160503_a"
# Capitaine, N., Wallace, P.T. & Chapront, J., 2005, A&A 432, 355

dayfrac = (timejd)%1
era00 = 2.*np.pi * (dayfrac + 0.7790572732640 + 0.00273781191135448 * (time - 2451545.0)) # Earth rotation angle - in radians

ttjd = timejd + deltat2
bigt = (ttjd - 2451545.0)/36525.0

gmstarcsec = 0.014506 + (4612.156534 + (1.3915817 + (-0.00000044 + (-0.000029956 + (-0.0000000368)* bigt) * bigt) * bigt) * bigt) * bigt

gmstrad = (gmstarcsec/3600.0/180.0*np.pi + era00)%(2.*np.pi)

'''
# Old version from 1982.
bigt = (timejd - 2451545.0) / 36525.0
gmstout = 24110.54841 + 8640184.812866 * bigt + 0.093104 * bigt * bigt - 0.0000062 * bigt * bigt * bigt
gmstout = gmstout / 86400.0
gmstout = gmstout%1
dayfrac = (time - 0.5)%1
gmst = gmstout + 1.0027379093 * dayfrac
gmst = gmst%1 * 24.0
gmstrad2 = gmst*np.pi/12.0
'''

# Get hour angle of Sun

lmst = gmstrad - tngwlon # Local mean sidereal time - in rad

ha = (lmst - ra2 + np.pi)%(2.*np.pi) - np.pi # hour angle - in radians

# Get airmass

cosz = np.sin(dec2)*np.sin(tnglat) + np.cos(dec2)*np.cos(tnglat)*np.cos(ha)
sinz = np.sqrt(1.-cosz*cosz)

air = 1./cosz
airyoung = air*(1.-0.0012*(air*air-1.)) # Correction for high airmass following Young 1994

# Get epsilon

zangle = np.arccos(cosz) # zenith angle

zp = zangle + angdia2/2. # zenith angle bottom
zm = zangle - angdia2/2. # zenith angle top

airp = 1./np.cos(zp) # airmass at bottom
airm = 1./np.cos(zm) # airmass at top

airyp = airp*(1.-0.0012*(airp*airp-1.)) # Young correction
airym = airm*(1.-0.0012*(airm*airm-1.)) # Young correction

dair = airyp - airym # difference in airmass between bottom and top (dair > 0)

dmag = kext * dair / 2. # diference in magnitude due to airmass and extinction (dmag > 0)

eps = dmag * np.log(10.0) / 2.5

# Get cospp and sinpp (position angle corrected for NPole position angle)

sinp = np.sin(ha)*np.cos(tnglat) / sinz
cosp = np.sqrt(1.-sinp*sinp)

sinpa = np.sin(nppa2)
cospa = np.cos(nppa2)

sinpp = sinp*cospa - cosp*sinpa # pp = p - pa, apply trig rules
cospp = cosp*cospa + sinp*sinpa


# get the velocity bias due to differential extinction and parallactic angle

##vbias = vrotsun * (7.*uu-15.0) * eps * sinpp / (20.0*(3.0-uu)) # Older simplified version

vbias = eps * sinpp * rsun * ( ((ASU-wearthobs)*(7*uu-15))/(20*(3-uu)) + (BSU*(19*uu-35+cosiobs*cosiobs*(3*uu-35)))/(280*(3-uu)) + (CSU*(187*uu-315-cosiobs**2.*(630-118*uu)-105*cosiobs**4.*(uu-1)))/(6720*(3-uu)) )

# Correct velocities

rv_final = rv_notss - vbias

# Save the data

with open("Sun_harpn_2-2-3_alldata_quality_extinct_corrections.csv", 'w') as f:

	f.write('root,bjd,rvraw,rvheli,rvfinal,erv,fwhm,fwhmcor,efwhm,contrast,contrastcor,econtrast,asym,easym,bisspan,ebisspan,exptime,airmass,sn10,sn50,sn60,expmetermean,expmetermin,expmetermax,expmeterstd,quality,k60,rvSS,rvext\n')

	for ind, value in enumerate(rv_final):

		f.write(data[ind,0]+','+data[ind,1]+','+data[ind,2]+','+str(rv_notss[ind])+','+str(value)+','+data[ind,3]+','+data[ind,4]+','+str(fwhmcor[ind])+','+data[ind,5]+','+data[ind,6]+','+str(contrastcor[ind])+','+data[ind,7]+','+data[ind,8]+','+data[ind,9]+','+data[ind,10]+','+data[ind,11]+','+data[ind,12]+','+data[ind,13]+','+data[ind,14]+','+data[ind,15]+','+data[ind,16]+','+data[ind,17]+','+data[ind,18]+','+data[ind,19]+','+data[ind,20]+','+data[ind,21]+','+data[ind,22]+','+str(rvbary_corr[ind])+','+str(vbias[ind])+'\n')

## Outliers removed

plt.subplot(5,1,1)

plt.plot(time,rv,'ko',ms=1)
plt.xlabel('Julian date (JDB - 2450000)')
plt.ylabel('RV (km/s)')

## Corrected for SS

plt.subplot(5,1,2)

plt.plot(time,rv_notss,'ko',ms=1)
plt.xlabel('Julian date (JDB - 2450000)')
plt.ylabel('RV (km/s)')

## Like plot 3, but phasebinned over a day
plt.subplot(5,1,3)
plt.scatter((time-0.5)%1, rv_notss, c=time, s=4, alpha=0.2)
plt.axhline(np.median(rv_notss))
p = np.polyfit((time-0.5)%1, rv_notss, 1)
xfit = np.arange(0.3,0.85,0.05)
yfit = p[0]*xfit + p[1]
plt.plot(xfit,yfit,'r-')
plt.xlim(0.3,0.8)
plt.xlabel('Fraction of day')
plt.ylabel('RV (km/s)')

## Velocity bias - phasebinned over a day
plt.subplot(5,1,4)
plt.scatter((time-0.5)%1, vbias + np.median(rv_notss), c=time, s=1)
plt.axhline(np.median(rv_notss))
p = np.polyfit((time-0.5)%1, vbias + np.median(rv_notss), 1)
xfit = np.arange(0.3,0.85,0.05)
yfit = p[0]*xfit + p[1]
plt.plot(xfit,yfit,'r-')
plt.xlim(0.3,0.8)
plt.xlabel('Fraction of day')
plt.ylabel('RV (km/s)')

## Final corrected velocities - phasebinned over a day
plt.subplot(5,1,5)
plt.scatter((time-0.5)%1, rv_final, c=time, s=4, alpha=0.2)
plt.axhline(np.median(rv_final))
p = np.polyfit((time-0.5)%1, rv_final, 1)
xfit = np.arange(0.3,0.85,0.05)
yfit = p[0]*xfit + p[1]
plt.plot(xfit,yfit,'r-')
plt.xlim(0.3,0.8)
plt.xlabel('Fraction of day')
plt.ylabel('RV (km/s)')


plt.show()
plt.close()




