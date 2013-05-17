################################ LICENSE ##################################
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.                                                     #
#                                                                          #
# Redistribution and use in source and binary forms, with or without       #
# modification, are permitted provided that the following conditions       #
# are met:                                                                 #
#                                                                          #
#     * Redistributions of source code must retain the above copyright     #
#       notice, this list of conditions and the following disclaimer.      #
#     * Redistributions in binary form must reproduce the above copyright  #
#       notice, this list of conditions and the following disclaimer       #
#       in the documentation and/or other materials provided with the      #
#       distribution.                                                      #
#     * Neither the name of the South African Astronomical Observatory     #
#       (SAAO) nor the names of its contributors may be used to endorse    #
#       or promote products derived from this software without specific    #
#       prior written permission.                                          #
#                                                                          #
# THIS SOFTWARE IS PROVIDED BY THE SAAO ''AS IS'' AND ANY EXPRESS OR       #
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED           #
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE   #
# DISCLAIMED. IN NO EVENT SHALL THE SAAO BE LIABLE FOR ANY                 #
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL       #
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS  #
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)    #
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      #
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN #
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE          #
# POSSIBILITY OF SUCH DAMAGE.                                              #
############################################################################

import time, datetime
from pyraf import iraf
from math import *
import numpy as np
from string import zfill
import calendar

def datetime2julian(date):
    """Converts a standard library `datetime' object to julian date float.
    The maximum precision is equal to that of the `datetime' object.
    """

    # January 1, 2000 at midday corresponds to JD = 2451545.0
    reference=datetime.datetime(year=2000,month=1,day=1,hour=12,minute=0,second=0,microsecond=0)

    temp=date-reference

    return 2451545+temp.days+(temp.seconds+temp.microseconds*1.e-6)/(24*3600)

# Convert UTC to Julian Date

def juliandate(year,month,day,hour,min,sec):

    tmp = float(day) + float(hour) / 24. + float(min) / 1440. + sec / 86400.
    if (month < 3):
        year -= 1
        month += 12
    a = year / 100
    b = 2 - a + a / 4
    c = int(365.25 * float(year))
    d = int(30.6001 * float(month + 1))
    julian = float(b + c + d) + tmp + 1720994.5 - 2453371.5 # midnight 2005/01/01

    return julian

# -----------------------------------------------
# Convert UTC to seconds

def utc2sec(hour,min,sec):

    time = float(hour) * 3600 + float(min) * 60 + sec

    return time

# -----------------------------------------------
# Convert TIME-OBS to hours

def time_obs2hr(time_obs):

    try:
        hour = float(time_obs[0:2]) + float(time_obs[3:5]) / 60 + float(time_obs[6:12]) / 3600
    except ValueError:
        hour = 99.999

    return hour

# -----------------------------------------------
# Convert from sexigesimal to decimal formats

def sex2dec(sex_value):
    value_arr=sex_value.split(':')
    sign=1
    if float(value_arr[0])<0: sign=-1
    dec_value=abs(float(value_arr[0]))+float(value_arr[1])/60.00+float(value_arr[2])/3600.0
    return sign*dec_value


# -----------------------------------------------
# Convert from decimal to sexigesimal formats

def dec2sex(dec_value):
    """Calculate the sexigesimal string from a decimal time

       return  string
    """

    sign=''
    if dec_value< 0: sign='-'
    dec_value=abs(dec_value)
    h=int(dec_value)
    m=int((dec_value-h)*60.0)
    s=(dec_value-h-m/60.0)*3600.0
    h='%02i' % h
    m='%02i' % m
    s='%06.3f' % s
    sex_value='%s:%s:%s' % (h, m, s)
    return sex_value



# -----------------------------------------------
# Convert TIME-OBS to YYYYMMDD

def date_obs2yyyymmdd(date_obs):
    monthName = {'01':'Jan', '02':'Feb', '03':'Mar', '04':'Apr', '05':'May', '06':'Jun',
              '07':'Jul','08':'Aug','09':'Sep','10':'Oct','11':'Nov','12':'Dec'}

    try:
        yyyymmdd = date_obs[0:4] + date_obs[5:7] + date_obs[8:10]
        year = date_obs[0:4]
        month=monthName[date_obs[5:7]]
        day = date_obs[8:10]
    except:
        year = '????'
        month = '???'
        day = '??'
        yyyymmdd='????????'

    date = day + ' ' + month + ' ' + year

    return yyyymmdd, date

# -----------------------------------------------
# Convert SALT date-time format to the DATETIME for mysql

def datatimeobs2DateTime(date):
    monthName = {'Jan':'01', 'Feb':'02', 'Mar':'03', 'Apr':'04', 'May':'05', 'Jun':'06',
              'Jul':'07','Aug':'08','Sep':'09','Oct':'10','Nov':'11','Dec':'12'}
    date=date.split()
    try:
        day=date[2]
        if int(day) < 10: day='0'+day
        mon=monthName[date[1]]
        year=date[4]
        time=date[3]
        newdate=year+'-'+mon+'-'+day+' '+time
    except:
        newdate=''

    return newdate

def currentobsdate():
   """Return the current obsdate"""
   yesterday = time.localtime(time.time())
   year = str(yesterday[0])
   month = str(yesterday[1])
   if yesterday[1] < 10:
       month = '0' + month
   day = str(yesterday[2])
   if yesterday[2] < 10:
       day = '0' + day
   return year + month + day


def breakdate(date):
    """For a date in the YYYYMMDD format, return day, month, year as integers
    """
    day=int(date[6:8])
    month=int(date[4:6])
    year=int(date[0:4])
    return day, month, year


def numberofdays(date):
    """Return the number of days for that month"""
    day, month, year=breakdate(str(date))
    return  np.array(calendar.monthcalendar(year,month)).max()
 

def getnextdate(date):
   """Give a date in YYYYMMDD format, it will supply the next YYYYMMDD date"""
   day,month,year=breakdate(str(date))
   tdate = datetime.datetime(year, month, day)
   tdate=tdate+datetime.timedelta(1)
   year=zfill(tdate.year, 4)
   month=zfill(tdate.month, 2)
   day=zfill(tdate.day, 2)
   return  year+month+day




def numLeapSeconds(year):
    # First need to calculate number of leap seconds. This is needed for JD and MJD.
    # The table for leap seconds is located at ftp : // maia.usno.navy.mil/ser7/tai - utc.dat
    # ***must be updated when new leap seconds are added!
    # Here is the portion we use (so no dates < 1999 are valid!):
    #1999 JAN 1 = JD 2451179.5 TAI - UTC = 32.0 S + (MJD - 41317.) X 0.0 S
    #2006 JAN 1 = JD 2453736.5 TAI - UTC = 33.0 S + (MJD - 41317.) X 0.0 S
    #2009 JAN 1 = JD 2454832.5 TAI - UTC = 34.0 S + (MJD - 41317.) X 0.0 S
    if 1999<=year<2006:
        leapSec=32.0
    elif 2006<=year<2009:
        leapSec=33.0
    elif 2009<=year:
        leapSec=34.0
    elif year<1999:
        print "Invalid time.  Please enter a date post-1999 for accurate leap year calculations."
    
    return leapSec
# -----------------------------------------------
# Return number of leap seconds. Valid only for dates post-1999. **MUST BE UPDATED WHEN NEW SECONDS ARE ADDED.

def convertUTtoJDUTC(date_obs):
    #This function converts UTC to Julian Date. This is the number of days since 01 January 4713 BC.
    #The formula used here is to convert Gregorian Calendar to JD from http://scienceworld.wolfram.com/astronomy/JulianDate.html.
    # The time fraction is coordinated universal time (UTC) timescale, which is slower than Terrestrial & Atomic times
    # (TT and TAI).
    
    y=float(date_obs[0:4])
    m=float(date_obs[5:7])
    d=float(date_obs[8:10])
    t=float(date_obs[10:12])+float(date_obs[13:15])/60+float(date_obs[16:])/3600
    
    jDUTC= 367*y-floor(7*(y + floor((m + 9)/12))/4) - floor(3*(floor((y + (m + 9)/7)/100)+1)/4) + floor(275*m/9)+d+1721028.5+t/24
    
    return jDUTC
    
# -----------------------------------------------
# Convert SALT date-time format (DATE-OBS+TIME-OBS) to JD (UTC)

def convertUTtoJD(date_obs):
    #This function converts UTC to Juilan Date and time fraction in Terrestrial Time, JD(TT).
    # TT is recommeded by IAU XXIII resolution B.1.
    jdUTC=convertUTtoJDUTC(date_obs)
    jdTT=jdUTC + (numLeapSeconds(float(date_obs[0:4]))+32.184)/(24*60*60) # The 32.184 sec is to convert from TAI to TT.
    
    return jdTT
# -----------------------------------------------
# Convert SALT date-time format (DATE-OBS+TIME-OBS) to JD (TT)

def convertUTtoMJD(date_obs):
    # This function converts UTC into Modified Julian Date and time fraction Terrestrial time.
    # This corresponds to days since midnight on 17 November 1858, or 2400000.5 JD.
    # We use TT as the timescale to be consistent with the IAU recommendation.
    jdUTC=convertUTtoJDUTC(date_obs)
    mjdUTC=jdUTC-2400000.5
    mjdTT=mjdUTC + (numLeapSeconds(float(date_obs[0:4]))+32.184)/(24*60*60) # The 32.184 sec is to convert from TAI to TT.
    
    return mjdTT
# -----------------------------------------------
# Convert SALT date-time format (DATE-OBS+TIME-OBS) to MJD (TT)

def hmsToRadians(raString):
    return 15*2*pi/360*(float(raString[0:2])+float(raString[3:5])/60+float(raString[6:])/3600)    
# -----------------------------------------------
# Convert SALT Right Ascension header format to radians

def dmsToRadians(decString):
    fractionalHour=(float(decString[1:3])+float(decString[4:6])/60+float(decString[7:])/3600)
    sign=decString[0]
    if sign=="-":
        outsign=-1
    elif sign=='+':
        outsign=1

    return outsign*2*pi/360*fractionalHour
# -----------------------------------------------
# Convert SALT Declination header format to radians

def ephInterp(ephfile,inputTime_JDTT):
    # This function reads in an ephemeris, interpolates for the input time(s),
    # and returns the the X,Y, and Z positional values. The times in the eph.
    # JD (TT) so the input time must match that format.
     
    # Currently, this function is only used for convertUTtoBJD. So the format
    # of the epemerides are specific to that function, and we return only
    # the required positional information. This function can be
    # modified in the future if other ephemeris types or outputs are needed.
    
    # To generate new eph. files in the proper format, send the following text
    # (edited to your email address and desired dates) to
    # horizons@ssd.jpl.nasa.gov with the word "job" (sans quotes) in the
    # subject line.  The object is altered by the value of
    # COMMAND: for example, '10' is the Sun and '399' is the Earth.
    # The ephemeris is assumed to be DE405. Note that the ephemerides
    # should NOT be light-time corrected.

    #!$$SOF (ssd) JPL/Horizons Execution Control VARLIST
    # !Oct 30,2002
    # !ftp://ssd.jpl.nasa.gov/pub/ssd/horizons_batch_example.brief
    # !!+++++++++++++++++++++++++++++++++++++++++++++++++++
    # !NOTE:First line in this file must start!$$SOF
    # !Last line in this file must start!$$EOF
    # !Assigned values should be in quotes
    # !+++++++++++++++++++++++++++++++++++++++++++++++++++
    # EMAIL_ADDR='amanda@saao.ac.za'
    # COMMAND='10'
    # OBJ_DATA='YES'
    # MAKE_EPHEM='YES'
    #   TABLE_TYPE='VECTORS'
    #     CENTER='500@0'
    #     REF_PLANE='FRAME'
    #     START_TIME='1999-Jan-01 00:00'
    #     STOP_TIME='2012-Jan-01 00:00'
    #     STEP_SIZE='1 day'
    #     REF_SYSTEM='J2000'
    #     OUT_UNITS='KM-S'
    #     VECT_TABLE='3'
    #     VECT_CORR='NONE'
    #     TIME_ZONE='+00:00'
    #     TIME_DIGITS='FRACSEC'
    #     RANGE_UNITS='AU'
    #     CSV_FORMAT='YES'
    #    VEC_LABELS='NO'
    #     R_T_S_ONLY='NO'
    #!$$EOF~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #The reply to this message will contain header and footer information as well
    #as the ephemerides.  *DELETE* the header and footer information (retaining
    #only the lines between $$SOE and $$EOE).  Then save the files in UNIX text
    #format as Object.eph in XXdirectory.
    
    ephData=np.loadtxt(ephfile,dtype=({'names':['JDTT','GregorianTime','X',\
                 'Y','Z','VX','VY','VZ','LightTime','Range','RangeRate'],\
                 'formats':[np.float,'S100',np.float,np.float,np.float,np.float,\
                 np.float,np.float, np.float,np.float,np.float]}),delimiter=',')
     
    (Xpos,Ypos,Zpos)=\
                        (np.interp(inputTime_JDTT,ephData['JDTT'],ephData['X']),\
                        np.interp(inputTime_JDTT,ephData['JDTT'],ephData['Y']),\
                        np.interp(inputTime_JDTT,ephData['JDTT'],ephData['Z']))
                
    return Xpos,Ypos,Zpos

# -----------------------------------------------
# Read in JPL ephemeris file and return interpolated X-Y-Z position at inputTime

def convertUTtoHJD(date_obs,ra,dec):
    # The Heliocentric Julian Date (HJD) is the Julian Date adjusted to the center of the Sun.
    # It depends on the JD of the observation (from which we get the position of Sun, in terms
    # of the solar longitude and the obliquity of the ecliptic), and the object RA and DEC.
    # The intent of HJD is to make a first-order accounting of the paralactic time shift
    # between the positions of the Earth and the Sun. 
    
    # The obliquity of the ecliptic calculation here is from the 2009 Astronomical Almanac, pg. C5,
    #"Low precision formulas for the Sun". These equations give the apparent coordinates of the Sun to
    # a precision of 0.01 degrees and the equation of time to 0.01 min between 1950 and 2050.
    
    jdt=convertUTtoJD(date_obs) # Use TT, so the final HJD is referenced to TT.
    obliqEclip=23.439-0.0000004*(jdt-2451545.0)
     
    # The longitude of the sun and formula for the distance/time the light must travel between the
    # Earth and Sun are from _Observational Astronomy_ by D.S Birney, Cambridge Univ. Press, 1991,
    # pages 248-251. Note that there is a typo in the final equation on p251.
    # This is supposedly accurate to 0.008s, the correction for the Earth's orbital
    # eccentricity; however, the obilquity of the ecliptic does have some error as well.
    # Tests with obliquity of +/- 0.01 deg show a change in resulting HJD of 0.01s.
    # Comparison with an online calculator shows ~7 sec. difference.

    #Convert ra and dec into radians:
    alpha= hmsToRadians(ra)
    delta= dmsToRadians(dec)
    
    T=(jdt-2451545.0)/36525
    L=280.460+36000.772*T
    M=357.528+35999.050*T
    solarLong=L+(1.915-0.0048*T)*sin(M) + 0.020*sin(2*M)
    
    degree=2*pi/360
    
    cosU=sin(delta)*sin(solarLong*degree)*sin(obliqEclip*degree)+cos(delta)*cos(alpha)*cos(solarLong*degree)+cos(delta)*sin(alpha)*sin(solarLong*degree)*cos(obliqEclip*degree)
    deltaT=8.3168775*cosU 
    hjdTT=jdt-deltaT/(60*24) #value is in minutes, so convert to days for subtraction from JD
    
    return hjdTT
# -----------------------------------------------
# Convert SALT date-time format (DATE-OBS+TIME-OBS) and object coordinates to Heliocentric JD (TT)

def convertUTtoBJD(date_obs,ra_string,dec_string):
    #    Barycentric JD is relative to the dynamical center-of-mass ("barycenter")
    # of the Solar System. There is not an easy formula to convert between this
    # time system and others. Here we follow the basic formula from C. Markwardt's
    # description at http://lheawww.gsfc.nasa.gov/users/craigm/bary/.
    # Unlike this page, we use the numpy interpolation function to interpolate the
    # ephemerides.
    #
    #     While we do perform 2 relativistic corrections (order of 10^-5/10^-6 sec),
    # the accuracy is not expect to be extremely high -- tests indicate
    # the accuracy is on the order of 10^(-5) seconds.
    #
    #       IMPORTANT NOTES:
    #       (1) The location is assumed to be the center of the Earth: no observatory
    #     information is included.
    #       (2) There is NO dispersion time correction.
    
    #      The total timing correction from our source is as follows:
    # tb = tobs+(clock)-(dispersion)+(geometric)+("Einstein")-("Shapiro"),
    # where  tobs - is the observed time of arrival;
    # tb - is barycentric arrival time;
    # (clock) - are the corrections which convert the local clock time to geocentric TT or TDT;
    # (dispersion) - are dispersion corrections of the form D/fb^2 where fb is the barycentric frequency.  D=0 for X-rays;
    # (geometric) - geometric time-delay in the solar system, aka "Romer Delay";
    # ("Einstein") - Einstein corrections - relativistic corrections of clock time to SSB (here just TDB-TT);
    # ("Shapiro") - Shapiro time delay due to photon bending in the potential of the solar system.
    #
    #     The "clock" correction is already done by using JD(TT) throughout the analysis,
    # we neglect the "dispersion" correction, and the remaining factors are calculated below.
    
    # Get input into useable formats, including date/time in JD (TT).
    y=float(date_obs[0:4])
    m=float(date_obs[5:7])
    d=float(date_obs[8:10])
    t=float(date_obs[10:12])+float(date_obs[13:15])/60+float(date_obs[16:])/3600
    
    ra=hmsToRadians(ra_string)
    dec= dmsToRadians(dec_string)
    
    inputJD=convertUTtoJD(date_obs)
    
    # Read Sun and Earth ephemerides and interpolate for the input time.
    sunFile=iraf.osfn("pysalt$data/ephem/sun.eph")
    earthFile=iraf.osfn("pysalt$data/ephem/earth.eph")
    (earthX,earthY,earthZ)=ephInterp(sunFile,inputJD) 
    (sunX,sunY,sunZ)=ephInterp(earthFile,inputJD) 
    
    # The geometric correction uses the components of the object coordinate
    # unit vector, and corrects for distance in terms of speed of light.

    xObj=cos(dec)*cos(ra)
    yObj=cos(dec)*sin(ra)
    zObj=sin(dec)

    speedOlight=2.99792458E8 # from IAU best estimates:http://maia.usno.navy.mil/NSFA/CBE.html
    geometricCorr=(earthX*xObj+earthY*yObj+earthZ*zObj)/(speedOlight/1000)
    
    # The Einstein correction is determined by the function to transform between
    # Terrestrial Time and Barycentric Dynamical Time (TDB) from the
    # 2009 Astronomical Almanac pg. B7.  TT and TDB differ due to
    # variations in gravitational potential around Earth's orbit. They differ
    # by no more than 2 millisec. The following equation is good to +/-30 microsec
    # between 1950 and 2050.
    # An additional Einstein correction could be made for the observer's location
    # on the surface of the Earth; however, this effect is even smaller than the
    # Shapiro correcation (msec level) so we neglect it here.
    
    g=(357.53+0.98560028*(inputJD-2451545.0))*(2*pi)/360
    LminusLj=(246.11+0.90251792*(inputJD-2451545.0))*(2*pi)/360
    einsteincorr=tbdminusTT=0.001657*sin(g)+0.000022*sin(LminusLj)

    # The Shapiro correction: between 20 microsec and 0.1 msec in effect.
    sunDist=sqrt(pow((sunX-earthX),2)+pow((sunY-earthY),2)+pow((sunZ-earthZ),2))
    cosTheta=((earthX-sunX)*xObj+(earthY-sunY)*yObj+(earthZ-sunZ)*zObj)/sunDist
    
    gravConst=6.67428E-11 # from IAU best estimates:http://maia.usno.navy.mil/NSFA/CBE.html
    sunMass=1.98892E30   # from IAU best estimates:http://maia.usno.navy.mil/NSFA/CBE.html
    massFactor=sunMass*gravConst/pow(speedOlight,3)
    shapirocorr=2*massFactor*log(1-cosTheta)
    
    # Sum all corrections and return.
    bjdMinusUTC=geometricCorr+einsteincorr+shapirocorr
    
    return inputJD+bjdMinusUTC /86400.
    
# -----------------------------------------------
# Convert SALT date-time format (DATE-OBS+TIME-OBS) to Barycentric JD (TT)
