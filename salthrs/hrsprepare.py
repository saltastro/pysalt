################################# LICENSE ##################################
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.   See LICENSE file for more details                 #
#                                                                          #
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

    if detname=='08443-03-01' or detname=='HRDET':
       detector='hrdet'
    elif detname=='04434-23-02' or detname=='HBDET':
       detector='hbdet'
    else:
       raise SaltError('%s is not an HRS detector' % detnam)

    #get key parameters
    try:
        nccd=saltkey.get('CCDAMPS', hdu[0])
    except:
        nccd=saltkey.get('NAMPS', hdu[0])
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
    saltkey.new('DETMODE', value='Normal', comment='Detector mode', hdu=nhdu[0])
    saltkey.new('NCCDS', value=1, comment='Detector mode', hdu=nhdu[0])
    saltkey.new('GAINSET', value='SLOW', comment='Detector mode', hdu=nhdu[0])
    value='OBJECT'
    if hdu[0].header['OBJECT']=='Bias': value='ZERO'
    if hdu[0].header['OBJECT']=='Arc': value='ARC'
    if hdu[0].header['OBJECT']=='Flat': value='FLAT'
    saltkey.new('CCDTYPE', value=value, comment='CCD Type', hdu=nhdu[0])
 
    if nccd==1:
        nhdu.append(pyfits.ImageHDU(hdu[0].data))
        j=1
        saltkey.new('GAIN', value=float(gain[0]), comment='Nominal CCD gain (e/ADU)', hdu=nhdu[j])
        saltkey.new('RDNOISE', value=0, comment='Nominal readout noise in e', hdu=nhdu[j])
        saltkey.new('SATURATE', value=1, comment='Pixel saturation level in ADU', hdu=nhdu[j])
        saltkey.new('XTALK',  value=0.0,  comment='Cross talk coefficient', hdu=nhdu[j])

        detsize=getdetsize(ccdshape, xbin, ybin) 
        ampshape=nhdu[j].data.shape
        ampsec=getdetsize(ccdshape, 1, 1)
        saltkey.new('DETSIZE',  value=detsize,  comment='Detector size', hdu=nhdu[j])
        saltkey.new('BIASSEC',  value=getbiassec(j, ampshape, xbin, ybin),  comment='Bias section', hdu=nhdu[j])
        saltkey.new('DATASEC',  value=getdatasec(j, ampshape, xbin, ybin),  comment='Data section', hdu=nhdu[j])
        saltkey.new('AMPSEC',  value=ampsec,  comment='Amplifier section', hdu=nhdu[j])
        saltkey.new('CCDSEC',  value=detsize,  comment='CCD section', hdu=nhdu[j])
        saltkey.new('DETSEC',  value=detsize,  comment='Detector section', hdu=nhdu[j])
      
        return nhdu


    for i in range(nccd):
        j=i+1
        x1,x2,y1,y2=definesection(j, detector, ccdshape)
        nhdu.append(pyfits.ImageHDU(hdu[0].data[y1:y2,x1:x2]))
        #keywords that need to be added include
        #gain, rdnoise, xtalk, saturate, 

        saltkey.new('GAIN', value=float(gain[i]), comment='Nominal CCD gain (e/ADU)', hdu=nhdu[j])
        saltkey.new('RDNOISE', value=0, comment='Nominal readout noise in e', hdu=nhdu[j])
        saltkey.new('SATURATE', value=1, comment='Pixel saturation level in ADU', hdu=nhdu[j])
        saltkey.new('XTALK',  value=0.0,  comment='Cross talk coefficient', hdu=nhdu[j])

        #DETSIZE, BIASSEC, DATASEC, 
        #AMPSEC, CCDSEC, DETSEC
        ampsec='[%i:%i,%i:%i]' % (x1+1,x2,y1+1,y2)
        ampshape=nhdu[j].data.shape
        saltkey.new('DETSIZE',  value=getdetsize(ccdshape, xbin, ybin),  comment='Detector size', hdu=nhdu[j])
        saltkey.new('BIASSEC',  value=getbiassec(j, ampshape, xbin, ybin),  comment='Bias section', hdu=nhdu[j])
        saltkey.new('DATASEC',  value=getdatasec(j, ampshape, xbin, ybin),  comment='Data section', hdu=nhdu[j])
        saltkey.new('AMPSEC',  value=ampsec,  comment='Amplifier section', hdu=nhdu[j])
        saltkey.new('CCDSEC',  value=getdetsize(ccdshape, xbin, ybin),  comment='CCD section', hdu=nhdu[j])
        saltkey.new('DETSEC',  value=getdetsize(ccdshape, xbin, ybin),  comment='Detector section', hdu=nhdu[j])

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

