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
SALT2IRAF--converts SALT MEF files to single extension fits file that 
are ready to work with basic IRAF tasks

 Author                 Version      Date
 -----------------------------------------------
 SM Crawford (SAAO)    1.0          26 Apr 2011

"""

# Ensure Python 2.5 compatibility
from __future__ import with_statement

import pyfits
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
   hdu=saltio.openfits(img)

   #if len one, copy and return
   if len(hdu)==1:  
      hdu.writeto(oimg)
      return

   #create the new output image
   odu = pyfits.PrimaryHDU(data=hdu[ext].data, header=hdu[0].header)

   #combine the headers from the primary and first exention
   for c in hdu[ext].header.ascardlist(): odu.header.update(c.key, c.value, c.comment)

   #write the data out
   saltio.writefits(odu, oimg, clobber=clobber)

# -----------------------------------------------------------
# main code

if not iraf.deftask('salt2iraf'):
   parfile = iraf.osfn("saltred$salt2iraf.par")
   t = iraf.IrafTaskFactory(taskname="salt2iraf",value=parfile,function=salt2iraf, pkgname='saltred')
