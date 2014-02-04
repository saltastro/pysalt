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

# Author                 Version      Date         Comment
# -----------------------------------------------------------------------
# S M Crawford (SAAO)    0.3          11 Oct 2013  

# hrsstack converts the MEF HRS files and stacks them together 
# into a single image 


import os, glob
import pyfits
import numpy as np
from pyraf import iraf
from pyraf.iraf import pysalt

import saltsafeio as saltio
import saltsafekey as saltkey
from saltsafelog import logging, history

from salterror import SaltError

debug=True

# -----------------------------------------------------------
# core routine

def hrsstack(images, outimages, outpref, clobber=True, logfile='salt.log',verbose=True):
    """Convert MEF HRS data into a single image.  If variance frames and BPMs, then 
       convert them to the same format as well.   Returns an MEF image but that is
       combined into a single frame

    """
 
    with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       #verify that the input and output lists are the same length
       saltio.comparelists(infiles,outfiles,'Input','output')

       #open the file and write out as a fits files in the output directory
       for img,oimg in zip(infiles, outfiles):
          hdu=saltio.openfits(img)
          hdu=stack(hdu)
          log.message('Stacking HRS %s to %s' % (img, oimg), with_header=False)
          saltio.writefits(hdu, oimg, clobber=clobber)
    
 
    return 

def stack(hdu):
    """Convert MEF HRS data into a single extension.   
    """

    #set the detector
    detname=saltkey.get('DETNAM', hdu[0])
    print detname
    if detname=='08443-03-01' or detname=='HRDET':
       detector='hrdet'
    elif detname=='04434-23-02' or detname=='HBDET':
       detector='hbdet'
    else:
       raise SaltError('%s is not an HRS detector' % detnam)
    print detector

    #get key parameters
    nccd=saltkey.get('CCDAMPS', hdu[0])

    #determine the shape of the CCD
    if detector=='hbdet':
       toty=hdu[1].shape[0]
       totx=hdu[1].shape[1]+hdu[2].shape[1]
    elif detector=='hrdet':
       toty=hdu[1].shape[0]+hdu[3].shape[0]
       totx=hdu[1].shape[1]+hdu[2].shape[1]
    data=np.zeros((toty,totx),np.float32)
    #embed the ccds
    for j in range(nccd):
        i=j+1
        y2,x2=hdu[i].data.shape
        if detector=='hbdet':
           if i==1:
              y1=0
              x1=0
              y2=y2
              x2=x2
           else:
              y1=0
              x1=x2
              x2=x1+x2
        elif detector=='hrdet':
           y1=0
           if i%2==1: x1=0
           if i%2==0: 
              x1=x2
              x2=x1+x2
           if i>2:
              y1=y2
              y2=y1+y2 
              
              
        data[y1:y2,x1:x2]=hdu[i].data
    
    ihdu = pyfits.ImageHDU(data)
    nhdu = pyfits.HDUList([hdu[0], ihdu])

    return nhdu


if not iraf.deftask('hrsstack'):
    parfile = iraf.osfn("salthrs$hrsstack.par")
    t = iraf.IrafTaskFactory(taskname="hrsstack",value=parfile,function=hrsstack, pkgname='salthrs')

