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

# salthrsprepare converts the raw hrs  files to fits and places HRS data
# into the standard SALT fits format and ready to run through 
# the rest of SALT pipeline tools


import os, glob
import pyfits
from pyraf import iraf
from pyraf.iraf import pysalt

import saltsafeio as saltio
import saltsafekey as saltkey
from saltsafelog import logging, history

from salterror import SaltError

debug=True

# -----------------------------------------------------------
# core routine

def hrsprepare(images, outimages, outpref, clobber=True, logfile='salt.log',verbose=True):
    """Convert .fit files to .fits files and place HRS data into 
       standard SALT FITS format

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
          hdu=prepare(hdu)
          log.message('Preparing HRS %s to %s' % (img, oimg), with_header=False)
          saltio.writefits(hdu, oimg, clobber=clobber)
    
 
    return 

def prepare(hdu):
    """Prepare HRS data to be similar to other SALT file formats
       This includes splitting each amplifier into multi-extension
       formats 
    """

    #set the detector
    detname=saltkey.get('DETNAM', hdu[0])

    if detname=='08443-03-01':
       detector='hrdet'
    elif detname=='04434-23-02':
       detector='hbdet'
    else:
       raise SaltError('%s is not an HRS detector' % detnam)

    #get key parameters
    nccd=saltkey.get('CCDAMPS', hdu[0])
    xbin, ybin=saltkey.ccdbin(hdu[0])
    gain=saltkey.get('GAIN', hdu[0])
    gain=gain.split()
    rospeed=saltkey.get('ROSPEED', hdu[0])
    ccdshape=hdu[0].data.shape

    #create a multiexention fits file
    phdu=pyfits.PrimaryHDU()
    phdu.header=hdu[0].header
 
    nhdu = pyfits.HDUList(phdu)

    #these keywords need to be added
    nhdu[0].header.set('DETMODE', value='Normal', comment='Detector mode')
    nhdu[0].header.set('NCCDS', value=1, comment='Detector mode')
    nhdu[0].header.set('GAINSET', value='SLOW', comment='Detector mode')

    for i in range(nccd):
        j=i+1
        x1,x2,y1,y2=definesection(j, detector, ccdshape)
        nhdu.append(pyfits.ImageHDU(hdu[0].data[y1:y2,x1:x2]))
        #keywords that need to be added include
        #gain, rdnoise, xtalk, saturate, 

        nhdu[j].header.set('GAIN', value=float(gain[i]), comment='Nominal CCD gain (e/ADU)')
        nhdu[j].header.set('RDNOISE', value=0, comment='Nominal readout noise in e')
        nhdu[j].header.set('SATURATE', value=1, comment='Pixel saturation level in ADU')
        nhdu[j].header.set('XTALK',  value=0.0,  comment='Cross talk coefficient')

        #DETSIZE, BIASSEC, DATASEC, 
        #AMPSEC, CCDSEC, DETSEC
        ampsec='[%i:%i,%i:%i]' % (x1+1,x2,y1+1,y2)
        ampshape=nhdu[j].data.shape
        nhdu[j].header.set('DETSIZE',  value=getdetsize(ccdshape, xbin, ybin),  comment='Detector size')
        nhdu[j].header.set('BIASSEC',  value=getbiassec(j, ampshape, xbin, ybin),  comment='Bias section')
        nhdu[j].header.set('DATASEC',  value=getdatasec(j, ampshape, xbin, ybin),  comment='Data section')
        nhdu[j].header.set('AMPSEC',  value=ampsec,  comment='Amplifier section')
        nhdu[j].header.set('CCDSEC',  value=getdetsize(ccdshape, xbin, ybin),  comment='CCD section')
        nhdu[j].header.set('DETSEC',  value=getdetsize(ccdshape, xbin, ybin),  comment='Detector section')

        #BSCALE, BZERO
        #nhdu[j].header.set('BSCALE', value=1.0, comment='Val=BSCALE*pix+BZERO')
        #nhdu[j].header.set('BZERO',  value=32768.,  comment='Val=BSCALE*pix+BZERO')
     
        #ATM1_1, ATM2_2, 
        #not implimented

    return nhdu

def getdatasec(namp, ampshape, xbin, ybin):
    """Set the detector size basec on the binning and shape of the ccd
 
       ccdshape: tuple
          Shape of the detector

       xbin: int
          x binning     

       ybin: int
          y binning     

    """
    #bias section on left hand side of chip
    if namp%2==1:
       x1=int(50/xbin)+1
       x2=ampshape[1]
       y1=1
       y2=ampshape[0]
    #bias section on left hand side of chip
    if namp%2==0:
       x1=1
       x2=ampshape[1]-50
       y1=1
       y2=ampshape[0]
    return '[%i:%i,%i:%i]' % (x1,x2,y1,y2)


def getbiassec(namp, ampshape, xbin, ybin):
    """Set the detector size basec on the binning and shape of the ccd
 
       ccdshape: tuple
          Shape of the detector

       xbin: int
          x binning     

       ybin: int
          y binning     

    """
    #bias section on left hand side of chip
    if namp%2==1:
       x1=int(10/xbin)
       x2=int(48/xbin)
       y1=1
       y2=ampshape[0]
    #bias section on left hand side of chip
    if namp%2==0:
       x1=ampshape[1]-50
       x2=x1+40
       y1=1
       y2=ampshape[0]
    return '[%i:%i,%i:%i]' % (x1,x2,y1,y2)



def getdetsize(ccdshape, xbin, ybin):
    """Set the detector size basec on the binning and shape of the ccd
 
       ccdshape: tuple
          Shape of the detector

       xbin: int
          x binning	

       ybin: int
          y binning	

    """
    x1=1
    y1=1
    x2=xbin*ccdshape[1]
    y2=xbin*ccdshape[0]
    return '[%i:%i,%i:%i]' % (x1,x2,y1,y2)

def definesection(namp, detector, ccdshape):
    """Define the data section for each amplifer

       namp: int
          Extension number of interest

       detector: str
          Name of the detector

       ccdshape: tuple
          Shape of the detector
       
    """
    if detector=='hrdet':
       #divide the image into four sections
       if namp==1:
          x1=0
          y1=0
       elif namp==2:
          x1=int(0.5*ccdshape[1])
          y1=0
       elif namp==3:
          x1=0
          y1=int(0.5*ccdshape[0])
       elif namp==4:
          x1=int(0.5*ccdshape[1])
          y1=int(0.5*ccdshape[0])
       x2 = x1 + int(0.5*ccdshape[1])
       y2 = y1 + int(0.5*ccdshape[0])
    if detector=='hbdet':
       #divide the image into two section
       if namp==1:
          x1=0
          y1=0
       elif namp==2:
          x1=int(0.5*ccdshape[1])
          y1=0
       x2 = x1 + int(0.5*ccdshape[1])
       y2 = y1 + ccdshape[0]
    return x1,x2,y1,y2

if not iraf.deftask('hrsprepare'):
    parfile = iraf.osfn("salthrs$hrsprepare.par")
    t = iraf.IrafTaskFactory(taskname="hrsprepare",value=parfile,function=hrsprepare, pkgname='salthrs')

