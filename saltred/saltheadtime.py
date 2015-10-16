#!/usr/bin/env python
#
# Author                 Version      Date
# -----------------------------------------------
# A. Gulbis (SAAO)        0.4          11 May 2011
# A. Gulbis (SAAO)        0.3          09 May 2011
# Amanda Gulbis (SAAO)    0.2          24 May 2010 
# Amanda Gulbis (SAAO)    0.1          17 May 2010              
#
# SALTTIMES is a tool to convert header times (UTC)
# into other date/time systems.
# After reading in the arguments, the algorithm is the
# following:
# 1.  Open the fits file(s)
# 2.  Read header time and target positions 
# 3.  Convert time into selected system
# 4.  Display header time and converted time
# 5.  Write new time to header [optional]
#
# All time conversion equations have been entered into the
# library file of salttime.py.
#
# Updates:
# ----------------------------------------------------------
# 0.4:  Revised to update all fits extension headers for slotmode data.
# 0.3:  Updated header keywords UTC-OBS-> TIME-OBS and OBJEPOCH-> EPOCH.
# 0.2: Added BJD conversion. Renamed saltheadtime from salttime, since the latter is a library.
#

# Ensure Python 2.5 compatibility
from __future__ import with_statement

#Python modules
import sys
import os

# Pyraf modules
from astropy.io import fits
from pyraf import iraf 
from pyraf.iraf import pysalt


# Project modules
import saltsafeio
from salttime import convertUTtoJDUTC, convertUTtoJD,convertUTtoMJD,convertUTtoHJD,convertUTtoBJD
from saltsafelog import logging
from salterror import SaltError, SaltIOError

# start main code
#-------------------------------------------------------------------------------
def saltheadtime(images,timetype,writetoheader,clobber,logfile,verbose,debug):
    # Start log.
    with logging(logfile,debug) as log:

        # is the input file specified?
        saltsafeio.filedefined('Input',images)

        # if the input file is a list, does it exist?
        if images[0] == '@':
            saltsafeio.listexists('Input',images)

        # parse list of input files
        infiles=saltsafeio.listparse('Raw image',images,'','','')

        # check input files exist
        saltsafeio.filesexist(infiles,'','r')
        
        # Loop over the input files.
        for file in infiles:
            print '---------------------------------------------------------\n\
            '
            print 'Reading file ', file
            openfile= saltsafeio.openfits(file)
            
            #Get information from header, display it, and close file.
            headers=openfile[0].header
            dateTimeString=[]
            ra=headers['RA']
            dec=headers['DEC']
            objEpoch=headers['EPOCH']
            obsEquinox=headers['EQUINOX'] 
                  
            if headers['DETMODE']=='SLOT':
               numExtensions=headers['NEXTEND']
               print 'Detector mode is slot. Number of extensions in this file is '+str(numExtensions)+'.'
               for num in range(numExtensions+1):
                    extheader=openfile[num].header
                    dateTimeString.append(extheader['DATE-OBS']+extheader['TIME-OBS'])       
            else:
                dateTimeString.append(headers['DATE-OBS']+headers['TIME-OBS'])                
            
            openfile.close()

            #Convert times. All equations are in library, salttime.py.
            newTime=[]
            if timetype=="JD(UTC)":
                for num in range(len(dateTimeString)):
                    newTime.append(convertUTtoJDUTC(dateTimeString[num]))
                    keyword='JDUT-OBS'
                    comment='Julian Date ref. to UTC.'
                    print 'The date and time from header extension '+ str(num) +' are:', dateTimeString[num], "UTC"
                    print 'Converted time ('+ timetype +') ='+str(newTime[num])+'\n\
                    '
            elif timetype=="JD(TT)":
                for num in range(len(dateTimeString)):
                    newTime.append(convertUTtoJD(dateTimeString[num]))
                    keyword='JD-OBS'
                    comment='Julian Date ref. to TT (Terrestrial timescale)'
                    print 'The date and time from header extension '+ str(num) +' are:', dateTimeString[num], "UTC"
                    print 'Converted time ('+ timetype +') ='+str(newTime[num])+'\n\
                    '
            elif timetype=="MJD(TT)":
                for num in range(len(dateTimeString)):
                    newTime.append(convertUTtoMJD(dateTimeString[num]))
                    keyword='MJD-OBS'
                    comment='Modified JD ref. to TT (Terrestrial timescale)'
                    print 'The date and time from header extension '+ str(num) +' are:', dateTimeString[num], "UTC"
                    print 'Converted time ('+ timetype +') ='+str(newTime[num])+'\n\
                    '
            elif timetype=="HJD(TT)":
                for num in range(len(dateTimeString)):
                    newTime.append(convertUTtoHJD(dateTimeString[num],ra,dec))
                    keyword='HJD-OBS'
                    comment='Heliocentric JD ref. to TT (Terrestrial timescale)'
                    print 'The date and time from header extension '+ str(num) +' are:', dateTimeString[num], "UTC"
                    print 'Converted time ('+ timetype +') ='+str(newTime[num])+'\n\
                    '
            elif timetype=="BJD(TDB)":
                for num in range(len(dateTimeString)):
                    newTime.append(convertUTtoBJD(dateTimeString[num],ra,dec))
                    keyword='BJD-OBS'
                    comment='Barycentric JD ref. to TDB (Barycentric Dynamical Time)'
                    print 'The date and time from header extension '+ str(num) +' are:', dateTimeString[num], "UTC"
                    print 'Converted time ('+ timetype +') ='+str(newTime[num])+'\n\
                    '
           
            
            #Open file to write new time to header, if desired.
            if writetoheader:
                try:
                    filetowrite=fits.open(file,mode='update')
                except Exception, e:
                    print e
                
                for num in range(len(dateTimeString)):
                    hdr=filetowrite[num].header
                    if hdr.has_key(keyword)==1:
                        if clobber:
                            hdr.update(keyword,newTime[num],comment)
                            print 'Keyword ' +keyword +' has been updated in header extension '+str(num)+' to the following value: '+ str(newTime[num])
                        else:
                            print 'Keyword '+keyword+' has not been updated.'
                    else:   
                        hdr.update(keyword,newTime[num],comment,after='TIME-OBS')
                        print 'Keyword ' + keyword +' has been added to header extension '+str(num)+' as the following value: '+ str(newTime[num])
                
                    filetowrite.flush()
                filetowrite.close()

# -----------------------------------------------------------
# end main code 
if not iraf.deftask('saltheadtime'):
   parfile = iraf.osfn("saltred$saltheadtime.par") 
   t =iraf.IrafTaskFactory(taskname="saltheadtime",value=parfile,function=saltheadtime,pkgname='saltred')
