#!/usr/bin/env python
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
"""
SPECSKY subtract the sky from a 2-D spectral image.

Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)    1.0       15 Nov 2010

TODO
----

LIMITATIONS
-----------

"""
# Ensure python 2.5 compatibility
from __future__ import with_statement

import os, sys
import time
import numpy as np

from pyraf import iraf 
import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging
import spectools as st
from spectools import SALTSpecError

from PySpectrograph.Spectra import apext

import pylab as pl

debug=True



# -----------------------------------------------------------
# core routine

def specsky(images,outimages,outpref, method='normal', section=None, 
            clobber=True,logfile='salt.log',verbose=True):

   with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       if method not in ['normal']: 
           msg='%s mode is not supported yet' % method
           raise SALTSpecError(msg)

       if section is None:
           section=saltio.getSection(section)
           msg='This mode is not supported yet'
           raise SALTSpecError(msg)
       else:
           section=saltio.getSection(section)

       # Identify the lines in each file
       for img, ofile in zip(infiles, outfiles):
          log.message('Subtracting sky spectrum in image %s into %s' % (img, ofile))
          #open the images
          hdu=saltio.openfits(img)

          #sky subtract the array
          hdu=skysubtract(hdu, method=method, section=section) 

          #write out the image
          if clobber and os.path.isfile(ofile): saltio.delete(ofile)
          hdu.writeto(ofile)



def skysubtract(hdu, method='normal', section=[]):
   """For a given image, extract a measurement of the sky from
      the image and then subtract that measurement from the 
      overall image
      and write the spectra to the output file

   """

   for i in range(len(hdu)):
       if hdu[i].name=='SCI':
           #set up the data, variance, and bad pixel frames
           #first step is to find the region to extract
           data_arr=hdu[i].data

           if saltkey.found('VAREXT', hdu[i]):
               var_ext=saltkey.get('VAREXT', hdu[i])
               var_arr=hdu[var_ext].data
           else:
               var_arr=None
               var_ext=None

           try:
               bpm_ext=saltkey.get('BPMEXT', hdu[i])
               bpm_arr=hdu[hdu[i].header['BPMEXT']].data
           except:
               bpm_ext=None
               bpm_arr=None
 
           #creat the xarr for the purposes of extracting the data
           xarr=np.arange(len(data_arr[0]))

           #TODO:  The variance array does not fully work at the moment
           var_arr=None

           #create the variance 
           ap=apext.apext(xarr, data_arr, ivar=var_arr)
  
           #extract the regions for the sky
           y1,y2=section
           nrows=abs(y2-y1)
           ap.flatten(y1,y2)
           ap.ldata=ap.ldata/nrows


           #subtract the data
           hdu[i].data=data_arr-ap.ldata

           #correct the variance frame
           if var_ext:
               hdu[var_ext].data=hdu[var_ext].data+ap.lvar/nrows

   return hdu

# main code 

parfile = iraf.osfn("saltspec$specsky.par") 
t = iraf.IrafTaskFactory(taskname="specsky",value=parfile,function=specsky, pkgname='saltspec')
