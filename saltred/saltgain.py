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

#!/usr/bin/env python
"""
SALTGAIN corrects the multiplies images by a constant factor appropriate for gain
correction. Each CCD in SALT's SALTICAM and RSS instruments has two
readout nodes, SALTICAM has a two-CCD mosaic and RSS has a three-CCD
mosaic. Each amplifier has a specific gain factor which varies slowly
over time but which is constant across the amplifiers. Gain values
depend on the readout speed and gain setting of the CCD. For all
possible permutations, gains are stored in an ascii table which is
updated periodically. Saltgain extracts gains from the ascii table and
applies them to raw data.

New gain and readout noise values will be written to the header
keywords of each HDU. Keyword writing can also occur without peforming
the gain correction itself. If the gain correction is performed, a
keyword, GAINMULT is added to the image extension with the value
1.0. If the gain correction is not performed, GAINMULT wil contain the
gain factor recorded in the ascii table. The purpose of the GAINMULT
keyword is to report what multiplicative factor is required to gain
correct an image.

Based on data in image keywords, e.g. gain setting, readout speed and
amplifier number, saltgain will extract the correct gain and readout
noise values from the ascii table and update keywords and optionally
perform the gain correction.  

Optionally, saltgain will extract the values from the data headers
and perform the gain correction using those values


Author                 Version      Date
-----------------------------------------------
Martin Still (SAAO)    1.0          31 Aug 2006
S. M. Crawford (SAAO)  2.0          20 May 2011

Updates
-----------------------------------------------
20110520 --the new error handling
         --able to read in values from the headers
         --able to write out to an outimages or outpref
         --ability to do non-linear gain correction


TODO
-----------------------------------------------
-Keywords and formats for non-linear gain correction
need to be finalized
-ability for config file to handle non-linear gains
"""

from __future__ import with_statement

from pyraf import iraf

import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging, history

from salterror import SaltError    

debug=True




# -----------------------------------------------------------
# core routine

def saltgain(images,outimages, outpref, gaindb=None,usedb=False, mult=True,
             clobber=True, logfile='salt.log',verbose=True):

   #start logging
   with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       #verify that the input and output lists are the same length
       saltio.comparelists(infiles,outfiles,'Input','output')

       # read in the database file if usedb is true
       if usedb:
           gaindb = gaindb.strip()
           dblist= saltio.readgaindb(gaindb)
       else:
           dblist=[]


       for img, oimg in zip(infiles, outfiles):
           #open the fits file
           struct=saltio.openfits(img)

           # identify instrument
           instrume,keyprep,keygain,keybias,keyxtalk,keyslot = saltkey.instrumid(struct)

           # has file been prepared already?
           if saltkey.found(keygain, struct[0]):
               message='SALTGAIN: %s has already been gain-corrected' % img
               raise SaltError(message)


           # gain correct the data
           struct = gain(struct,mult=mult, usedb=usedb, dblist=dblist, log=log, verbose=verbose)

           # housekeeping keywords
           fname, hist=history(level=1, wrap=False, exclude=['images', 'outimages', 'outpref'])
           saltkey.housekeeping(struct[0],keygain, 'Images have been gain corrected', hist)

           # write FITS file
           saltio.writefits(struct,oimg, clobber=clobber)
           saltio.closefits(struct)

# -----------------------------------------------------------
# gain correct data

def gain(struct,mult=True,usedb=False, dblist=None, ampccd=2, log=None, verbose=True):
   """gain processes a image hduList and gain corrects each amplifier.  It can
      either use gain settings in the header or those supplied in a config file
      which would be suppleid in the dblist (see helpfile for structure of the
      config file).  If variance frames exist, it will update those for changes 
      in the header value as well.   

      In the end, it will update the gain with a value of one signfing the data
      has been transformed into e- from ADU

      The program will look for the non-linear gain settings which are given by:
      e  = GAIN*(1 + GAIN1*E-6*ADU)*ADU

      mult--if true, multiple the gains
      usedb--use the values in the dblist, if false use the header values
      dblist--values for the gain and readnoise from the 
      ampccd--number of amplifiers per ccd

      dblist should have the following lists:  speed, rate, gain, noise, bias, amp
   """
   #get the infile name
   infile=saltkey.getimagename(struct[0])

   #how many science extensions
   nsciext = saltkey.get('NSCIEXT',struct[0])
   #how many data extensions
   nextend = saltkey.get('NSCIEXT',struct[0])

   # how many amplifiers?
   amplifiers = ampccd*saltkey.get('NCCDS',struct[0])

   #read the gain and rospeed for the image
   gainset = saltkey.get('GAINSET',struct[0])
   rospeed = saltkey.get('ROSPEED',struct[0])

   #loop through each amplifier and gain correct it
   if log:
        message = '%28s %6s %5s %3s %5s %5s' \
            % ('HDU','GAIN','SPEED','AMP','GAIN','NOISE')
        log.message('\n      ---------------------------------------------------', \
                    with_header=False, with_stdout=verbose)
        log.message(message, with_header=False, with_stdout=verbose)
        log.message('      ---------------------------------------------------', \
                    with_header=False, with_stdout=verbose)
   for i in range(nsciext):
       hdu = i + 1
       amp = i%amplifiers+1
       #get the gain and rdnoise values for the array 
       if usedb:
           gain, rdnoise=get_values(dblist, gainset, rospeed, amp)
           gain1=0
       else:
           gain = saltkey.get('GAIN',struct[hdu])
           rdnoise = saltkey.get('RDNOISE',struct[hdu])
           try:
               gain1=saltkey.get('GAIN1',struct[hdu])
           except:
               gain1=0
 

       if mult:
           #correct the gain
           gainmult=1
           try:
               data=struct[hdu].data
               struct[hdu].data=gain*data+gain1*data**2
           except Exception, e:
               msg='Cannot gain correct %s[%i] because %s' % (infile, hdu, e)
               raise SaltError(msg)
               
           #correct the variance frame
           if saltkey.found('VAREXT', struct[hdu]):
               vhdu=saltkey.get('VAREXT', struct[hdu])
               try:
                   vdata=struct[vhdu].data
                   struct[vhdu].data=vdata*gain*(1+2*gain1*1e-6*data)
               except Exception, e:
                    msg='Cannot update the variance frame in %s[%i] because %s' % (infile, vhdu, e)
                    raise SaltError(msg)
       else:
           gainmult=gain
 
       #update the headers
       if usedb:
           saltkey.put('GAIN',gain,struct[hdu])
           saltkey.put('RDNOISE',rdnoise,struct[hdu])

       #add a keyword indicating what action was taken
       saltkey.new('GAINMULT',gainmult,'Gain multiplication', struct[hdu])

       #if logging is true, then print out the following information
       if log:
           message = '%25s[%1d] %6s %5s %2s %6.2f %5.2f' \
             % (infile,hdu,gainset,rospeed,amp, gain, rdnoise)
           log.message(message, with_header=False, with_stdout=verbose)

   #just to make it look pretty 
   if log:
       log.message('',  with_header=False, with_stdout=verbose)
          
   return struct

def get_values(dblist, gainset, rospeed, amp):
   """Get values for gain and rdnoise from the dblist.  The input of the dblist should be:
      dblist should have the following lists:  speed, rate, gain, noise, bias, amp
   """
   #default values in case none are found
   gain=1
   rdnoise=5
   found=False
   dbspeed, dbrate, dbgain, dbnoise, dbbias, dbamp=dblist

   #checkf or the gain
   if gainset not in dbrate:
       raise SaltError('GAINSET=%s does not match any vaue in gaindb file' % gainset)
   if rospeed not in dbspeed:
       raise SaltError('ROSPEED=%s does not match any vaue in gaindb file' % rospeed)
     
   #loop through the values and find the corresponding gain/rdnoise
   for i in range(len(dbspeed)):
       if dbspeed[i].strip().upper()==rospeed.strip().upper():
           if dbrate[i].strip().upper()==gainset.strip().upper():
               if int(dbamp[i])==int(amp):
                  gain=float(dbgain[i])
                  rdnoise=float(dbnoise[i])
                  found=True
                  return gain, rdnoise

   if not found:
       msg='Could not find a corresponding setting in the gaindb file for gainset=%s, rospeed=%s, amp=%i' \
            % (gainset, rospeed, amp)
       raise SaltError(msg)
                  
   return rdnoise, gain

# -----------------------------------------------------------
# main code
if not iraf.deftask('saltgain'):
   parfile = iraf.osfn("saltred$saltgain.par")
   t = iraf.IrafTaskFactory(taskname="saltgain",value=parfile,function=saltgain, pkgname='saltred')
