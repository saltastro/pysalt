################################# LICENSE ##################################
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

"""SLOTREADTIME is a tool to fix the UTC time.  The time in the binary files is 
the readtime and not the UTC time.  The actually time of the exposure is 
7xEXPTIME less than the time of the read out as the image is shifted down 7 
places on the CCD before being read out.

This program will read in a file, and for each extension correct the UTC-TIME 
keyword by the time it takes to shift the exposure so the UTC-TIME corresponds 
to the start of the exposure and adds a keyword with READTIME. 

It will not run on files with the READTIME header keyword already in place.

"""

# Ensure python 2.5 compatibility
from __future__ import with_statement

import os
import sys
import time, datetime

from pyraf import iraf
from iraf import pysalt

import saltsafekey as saltkey
import saltsafeio as saltio
import salttime
import slottool
from saltsafelog import logging, history
from salterror import SaltError, SaltIOError


debug=True

def slotreadtimefix(images,outimages, outpref, 
                    clobber=False, logfile='slot.log',verbose=True):
    
    with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       #verify that the input and output lists are the same length
       saltio.comparelists(infiles,outfiles,'Input','output')

       for img, oimg in zip(infiles, outfiles):
           #check to see if the out image already exists
           if not clobber and os.path.isfile(oimg):
              raise SaltIOError('%s alraedy exists' % oimg)

           #open the file
           struct=saltio.openfits(img)

           #log the message
           log.message('Updateing times in %s' % img, with_header=False, with_stdout=verbose)

           #now for each science frame, corrent the readtime
           #Assumes SALT format and that the first extension 
           #is empty
           for i in range(1,len(struct)):
              try:
                  struct[i]=readtimefix(struct[i])
              except SaltIOError,e :
                  raise SaltError('%s %s' % (img,e))

           #Add history keywords
           # housekeeping keywords
           fname, hist=history(level=1, wrap=False, exclude=['images', 'outimages', 'outpref'])
           saltkey.housekeeping(struct[0],"SLOTREAD", 'READTIME added', hist)


           #write to the output
           saltio.writefits(struct, oimg, clobber)

       return

def readtimefix(hdu, dsteps=7, transtime=4e-3):
    """Update the hdu with the correct time for when the exposure started 
       and add the READTIME keyword

       dsteps--the number of readouts to correct for
       transtime--the transfer time between each frame
    """

    #check for if the data has already been processed
    if saltkey.found('READTIME', hdu):
        raise SaltIOError(' has already been processed')

    #determine the UTC time 
    utctime=saltkey.get('UTC-OBS', hdu)
    timeobs=saltkey.get('TIME-OBS', hdu)
    dateobs=saltkey.get('DATE-OBS', hdu)
    exptime=float(saltkey.get('EXPTIME', hdu))

    #add the readtime header
    saltkey.new("READTIME",utctime,'Time of the readout of the frame', hdu)

    #correct the utctime--first switch to datetime to properly handle
    #dates around changes in hours
    y,m,d=dateobs.split('-')
    H,M,S=utctime.split(':')
    s=int(float(S))
    ms=int(1e6*(float(S)-s))
    newtime=datetime.datetime(int(y),int(m),int(d),int(H),int(M),s,ms)

    #correct the datetime
    dtime=dsteps*(exptime+transtime)
    s=int(dtime)
    ms=int(1e6*(dtime-s))
    newtime=newtime-datetime.timedelta(0, s, ms)

    #update the headkeywords
    hdu.header["UTC-OBS"]=str(newtime.time())
    saltkey.put("UTC-OBS", str(newtime.time()), hdu)
    saltkey.put("TIME-OBS", str(newtime.time()), hdu)
    saltkey.put("DATE-OBS", str(newtime.date()), hdu)

    return hdu

# -----------------------------------------------------------
# main code

parfile = iraf.osfn("slottools$slotreadtimefix.par")
t = iraf.IrafTaskFactory(taskname="slotreadtimefix",value=parfile,function=slotreadtimefix,pkgname='slottools')
