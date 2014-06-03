################################# LICENSE ##################################
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.                                                     #
#                                                                          #
############################################################################

#!/usr/bin/env python

"""
SALT2IRAF--converts SALT MEF files to single extension fits file that 
are ready to work with basic IRAF tasks

 Author                 Version      Date
 -----------------------------------------------
 SM Crawford (SAAO)    1.0          26 Apr 2011

"""

# Ensure Python 2.5 compatibility
from __future__ import with_statement

from astropy.io import fits
from pyraf import iraf
from pyraf.iraf import pysalt

# Salt imports
import saltsafeio as saltio
from saltsafelog import logging

from salterror import SaltError

debug=True

# -----------------------------------------------------------
# core routine

def salt2iraf(images,outimages,outpref,ext=1, clobber=True,logfile='salt.log',
              verbose=True):

   with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       #verify that the input and output lists are the same length
       saltio.comparelists(infiles,outfiles,'Input','output')

       for img, oimg in zip(infiles, outfiles):
           message="SALT2IRAF--Converting %s to %s" % (img, oimg)
           log.message(message, with_header=False)
           try:
               convertsalt(img, oimg, ext=ext, clobber=clobber)
           except SaltError, e:
               log.message('%s' % e)

def convertsalt(img, oimg, ext=1, clobber=True):
   """Convert a SALT MEF file into a single extension fits file"""

   #open the image
   hdu=fits.open(img)

   #if len one, copy and return
   if len(hdu)==1:  
      hdu.writeto(oimg)
      return

   #create the new output image
   odu = fits.PrimaryHDU(data=hdu[ext].data, header=hdu[0].header)

   #combine the headers from the primary and first exention
   for c in hdu[ext].header.cards: 
       odu.header.set(c.keyword, c.value, c.comment)

   #write the data out
   saltio.writefits(odu, oimg, clobber=clobber)

# -----------------------------------------------------------
# main code

if not iraf.deftask('salt2iraf'):
   parfile = iraf.osfn("saltred$salt2iraf.par")
   t = iraf.IrafTaskFactory(taskname="salt2iraf",value=parfile,function=salt2iraf, pkgname='saltred')
