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

# Author                 Version      Date
# -----------------------------------------------
# Martin Still (SAAO)    1.0          25 Jul 2007
# S M Crawford (SAA0)    1.1          16 Jul 2009
# S M Crawford (SAA0)    1.2          19 Apr 2011

# saltobslog reads the critical header keywords of SALT FITS data and
# collates them into a FITS table. This can be used for the basis
# of an observation log or as a meta-table for pipeline processing.
#
# Updates:
# 16 Jul 2009   Handled a bug on how if data is missing in the headers
# 5 Apr 2011    Updated to handle new header key words and changes in existing 
#               keywords. We did not make it backwards compatible
#               Changes made:
#               TELTEMP->TELTEM
#               PAYLTEMP -> PAYLTEM
#               Removed DETSIZE
#               Added EPOCH
#               Removed UTC-OBS--now using TIME-OBS
#               FOCUS->CAMFOCUS
#               INTFR->TELFOCUS
# 19 Apr 2011   Updated to handle the new error handling and so it returns the 
#               contents without updating everything in the process
#               Converted it to using dictionary instead of a bunch of lists

from __future__ import with_statement

from pyraf import iraf
import os, glob, time 
from astropy.io import fits

import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging


from salterror import SaltError, SaltIOError


# -----------------------------------------------------------
# core routine

headerList=['FILENAME', 'PROPID', 'PROPOSER', 'OBJECT', 'RA', 'DEC', 'OBJEPOCH', 'EPOCH', 'EQUINOX', 'DATE-OBS', 'UTC-OBS', 'TIME-OBS', 'EXPTIME', 'OBSMODE', 'DETMODE', 'CCDTYPE', 'DETSIZE', 'NCCDS', 'CCDSUM', 'GAINSET', 'ROSPEED', 'INSTRUME', 'FILTER', 'CAMFOCUS', 'TELHA', 'TELRA', 'TELDEC', 'TELPA', 'TELAZ', 'TELALT', 'TRKX', 'TRKY', 'TRKZ', 'TRKPHI', 'TRKTHETA', 'TRKRHO', 'TELFOCUS', 'COLPHI', 'COLTHETA', 'TELTEM', 'PAYLTEM', 'CCDTEM', 'DEWTEM', 'AMPTEM', 'CENTEM', 'DETSWV', 'BLOCKID', 'BVISITID']
formatList=['32A', '50A', '20A', '100A', '12A', '12A', 'E', 'E', 'I', '10A', '12A', 
            '12A', 'D', '20A', '20A', '8A', '23A', 'I', '5A', '6A', '4A', '8A', '8A',
            'J', '11A', '11A', '12A', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E',
            'E', 'E', 'E', 'E', 'E', 'E', 'D', 'E', 'E', 'D', '16A', 'A', 'A' ]
rssheaderList=['DEWPRE', 'POSANG', 'LAMPID', 'CALFILT', 'CALND', 'TELRHO', 'PELLICLE', 'INSTPORT', 'CF-STATE', 'SM-STATE', 'SM-STA', 'SM-STEPS', 'SM-VOLTS', 'SM-STA-S', 'SM-STA-V', 'MASKID', 'MASKTYP', 'WP-STATE', 'HWP-CMD', 'HW-STEPS', 'HWP-STA', 'QWP-CMD', 'QW-STEPS', 'QWP-STA', 'QWP-ANG', 'HWP-ANG', 'SH-STATE', 'FO-STATE', 'FO-POS', 'FO-VOLTS', 'FO-POS-S', 'FO-POS-V', 'GR-STATE', 'GR-STA', 'GR-ANGLE', 'GM-STEPS', 'GM-VOLTS', 'GR-STA-S', 'GR-STA-V', 'GR-STEPS', 'GRATING', 'GRTILT', 'BS-STATE', 'FI-STATE', 'FI-STA', 'FM-STEPS', 'FM-VOLTS', 'FM-STA-S', 'FM-STA-V', 'AR-STATE', 'AR-STA', 'CAMANG', 'AR-STA-S', 'AR-ANGLE', 'COLTEMP', 'CAMTEMP', 'PROC', 'PCS-VER', 'WPPATERN']
rssformatList=['D', 'E', '8A', '8A', 'E', 'E', '8A', '8A', '20A', '20A', '8A', 'I', 'E', 'E', 'E', '16A', '16A', '20A', '16A', 'I', 'E', '16A', 'I', 'E', 'E', 'E', '20A', '20A', 'E', 'E', 'E', 'E', '20A', '10A', 'E', 'I', 'E', 'E', 'E', 'I', '8A', 'E', '24A', '20A', '7A', 'I', 'E', 'E', 'E', '24A', '16A', 'E', 'E', 'E', 'E', 'E', '20A', '4A', '20A']
scamheaderList=['FILPOS']
scamformatList=['I']

debug=True

def saltobslog(images,outfile,clobber=False,logfile='salt.log',verbose=True):
  """Create the observation log from the input files"""

  #start the logging
  with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)


       #create the header dictionary
       headerDict=obslog(infiles, log)

       #clobber the output if it exists
       if (clobber and os.path.isfile(outfile)):
           saltio.delete(outfile)

       #create the fits file
       struct=createobslogfits(headerDict)

       # close table file
       saltio.writefits(struct, outfile)

       #indicate the log was created
       log.message('\nSALTLOG -- created observation log ' + outfile)

def createobslogfits(headerDict):
   """Create the fits table for the observation log"""
   # define generic columns of output table
   col=[]
   for k, f in zip(headerList, formatList):
       print k,f, headerDict[k]
       col.append(fits.Column(name=k, format=f, array=headerDict[k]))
   for k, f in zip(scamheaderList, scamformatList):
       print k,f, headerDict[k]
       col.append(fits.Column(name=k, format=f, array=headerDict[k]))
   for k, f in zip(rssheaderList, rssformatList):
       print k,f, headerDict[k]
       col.append(fits.Column(name=k, format=f, array=headerDict[k]))

   # construct FITS table from columns
   table = fits.ColDefs(col)

   # write FITS table to output file
   struct = fits.BinTableHDU.from_columns(table)

   # name the table extension
   struct.header['EXTNAME'] = 'OBSLOG'
   struct.header['SAL_TLM'] = time.asctime(time.localtime())
   #saltkey.new('EXTNAME','OBSLOG','extension name', struct)
   #saltkey.put('SAL-TLM',time.asctime(time.localtime()), struct)

   # housekeeping keywords

   return struct

   


# -----------------------------------------------------------
# read keyword and append to list

def obslog(infiles, log=None):
   """For a set of input files, create a dictionary contain all the header 
      information from the files.  Will print things to a saltlog if log is
      not None
    
      returns Dictionary
   """
   #create the header dictionary
   headerDict={}
   for k in headerList: headerDict[k]=[]
   for k in scamheaderList: headerDict[k]=[]
   for k in rssheaderList: headerDict[k]=[]

   # interate over and open image files
   infiles.sort()
   for infile in infiles:

       #open the file
       struct = saltio.openfits(infile)

       # instrument
       scam = False
       rss = False
       instrume = saltkey.get('INSTRUME', struct[0])
       if (instrume=='RSS'): rss = True
       if (instrume=='SALTICAM'): scam=True

       #add in the image name
       headerDict['FILENAME'].append(os.path.basename(infile))

       # ingest primary keywords from files in the image list
       for k,f in zip(headerList[1:], formatList[1:]):
           default=finddefault(f)
           headerDict[k].append(getkey(struct[0], k, default=default, log=log, warn=True))

       # ingest scam specific primary keywords from files in the image list
       for k,f in zip(scamheaderList[1:], scamformatList[1:]):
           default=finddefault(f)
           headerDict[k].append(getkey(struct[0], k, default=default, log=log, warn=scam))

       # ingest rss specific primary keywords from files in the image list
       for k,f in zip(rssheaderList[1:], rssformatList[1:]):
           default=finddefault(f)
           headerDict[k].append(getkey(struct[0], k, default=default, log=log, warn=rss))

       # close image files
       saltio.closefits(struct)

       if log: log.message('SALTOBSLOG -- read %s' % infile, with_header=False)

   return headerDict

def finddefault(f):
   """return the default value given a format"""
   if f.count('A'): 
       default="UNKNOWN"
   elif f.count('I'):
       default=-999
   else: 
       default=-999.99
   return default

def getkey(struct,keyword,default,warn=True, log=None):
   """Return the keyword value.  Throw a warning if it doesn't work """

   try:
        value = saltkey.get(keyword, struct)
        if isinstance(default, str):  value=value.strip()
   except SaltIOError:
        value = default
        infile=struct._file.name
        message = 'WARNING: cannot find keyword %s in %s' %(keyword, infile)
        if warn and log: log.message(message, with_header=False)
   if (str(value).strip() == ''): value = default
   if (type(value) != type(default)):
        infile=struct._file.name
        message='WARNING: Type mismatch for %s for  %s in %s[0]' % (str(value), keyword, infile)
        message += '/n '+str(type(value)) + ' '+str(type(default))
        if warn and log: log.message(message, with_header=False)
        value=default

   return value


# -----------------------------------------------------------
# main code
if not iraf.deftask('saltobslog'):
   parfile = iraf.osfn("saltred$saltobslog.par")
   t = iraf.IrafTaskFactory(taskname="saltobslog",value=parfile,function=saltobslog, pkgname='saltred')
