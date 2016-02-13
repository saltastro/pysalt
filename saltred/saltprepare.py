################################# LICENSE ##################################
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.                                                     #
#                                                                          #
############################################################################

#!/usr/bin/env python
"""
SALTPREPARE--Prepares the SALT data for pipeline processing.  The main
purpose of this task is to verify the correct header key words and the
data file format.  

Saltprepare is currently a relatively trivial task that adds several
missing but required keywords to the keyword headers of FITS data
files. These keywords are not written to the file at the telescope but
are required in order to improve efficiency in the data reduction
pipeline and retain consistency with the IRAF-based tools.

The user may also select to create variance and bad pixel frames for 
the data by setting the createvar varialble.  If selected, variance
frames will be created based on the pixel variance and the readnoise.
The user can also supply a bad pixel frame, but if one is not supplied,
one can be created from the data itself assuming all pixels
are good.


Author                 Version      Date
-----------------------------------------------
Martin Still (SAAO)    1.0          21 Jul 2006
S. M. Crawford (SAAO)  2.0          18 May 2011

TODO
-----------------------------------------------
# The number of amplifies is hardwired at 2 into the code

UPDATES
-----------------------------------------------
18 May 2011  -Updated to use new error handling
             -updated to create variance and bad pixel frames
"""

from __future__ import with_statement

import numpy as np 

from pyraf import iraf
from pyraf.iraf import pysalt
import os, string, sys, glob, time
from astropy.io import fits

import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging, history

from salterror import SaltError    

debug=True


# -----------------------------------------------------------
# core routine

def saltprepare(images,outimages,outpref,createvar=False, badpixelimage=None, clobber=True,logfile='salt.log',verbose=True):

   with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       #verify that the input and output lists are the same length
       saltio.comparelists(infiles,outfiles,'Input','output')

       # open the badpixel image
       if saltio.checkfornone(badpixelimage) is None:
          badpixelstruct=None
       else:
           try:
               badpixelstruct = saltio.openfits(badpixelimage)
           except saltio.SaltIOError,e:
               msg='badpixel image must be specificied\n %s' % e
               raise SaltError(msg)

       # open each raw image file
       for img, oimg, in zip(infiles, outfiles):

           #open the fits file
           struct=saltio.openfits(img)

           #if createvar, throw a warning if it is using slotmode
           if saltkey.fastmode(saltkey.get('DETMODE', struct[0])) and createvar:
               msg='Creating variance frames for slotmode data in %s' % img
               log.warning(msg)
  
           # identify instrument
           instrume,keyprep,keygain,keybias,keyxtalk,keyslot = saltkey.instrumid(struct)

           # has file been prepared already?
           try:
               key = struct[0].header[keyprep]
               message = 'ERROR -- SALTPREPARE: File ' + infile
               message += ' has already been prepared'
               raise SaltError(message)
           except:
               pass


           # prepare file
           struct=prepare(struct,createvar=createvar, badpixelstruct=badpixelstruct)

           # housekeeping keywords
           fname, hist=history(level=1, wrap=False, exclude=['images', 'outimages', 'outpref'])
           saltkey.housekeeping(struct[0],keyprep, 'File prepared for IRAF', hist)

           # write FITS file
           saltio.writefits(struct,oimg, clobber=clobber)
           saltio.closefits(struct)

           message = 'SALTPREPARE -- %s => %s' % (img, oimg)
           log.message(message, with_header=False)



# -----------------------------------------------------------
# prepare FITS file for processing

def prepare(struct,createvar=False, badpixelstruct=None, namps=2):

   #set up the file name
   infile=saltkey.getimagename(struct[0])

   #include information about which headers are science extensoisn
   #and check to make sure each array is a data array
   nextend=len(struct)-1
   nsciext=nextend
   for i in range(1,nextend+1):
       try:
           struct[i].header['XTENSION'] == 'IMAGE'
       except:
           msg='Extension %i in %s is not an image data'
           raise  SaltError(msg)
       saltkey.new('EXTNAME','SCI','Extension name',struct[i])
       saltkey.new('EXTVER',i,'Extension number',struct[i])

   #check the number of amplifiers
   #TODO:  This is current hardwired at a value of 2
   nccds = saltkey.get('NCCDS',struct[0])
   if (nextend%(nccds*namps) != 0):
        message = 'ERROR -- SALTPREPARE: Number of file extensions and '
        message += 'number of amplifiers are not consistent in ' + infile
        raise SaltError(message)

   #Add the inv. variance frame if requested--assumes no variance frame
   #and the science frames are in the first n frames
   if createvar:
       #create the inv. variance frames
       for i in range(1, nsciext+1):
           try:
               hdu=CreateVariance(struct[i], i, nextend+i)
           except Exception, e: 
               msg='Cannot create variance frame in extension %i of  %s because %s' % (nextend+i, infile, e)
               raise SaltError(msg) 
           struct[i].header.update('VAREXT',nextend+i, comment='Extension for Variance Frame')
           struct.append(hdu)
       nextend+=nsciext

       #create the badpixelframes
       for i in range(1, nsciext+1):
           try:
               hdu=createbadpixel(struct, badpixelstruct, i, nextend+i)
           except Exception, e: 
               msg='Could not create bad pixel extension in ext %i of %s because %s' % (nextend+i, infile, e)
               raise SaltError(msg) 
           struct[i].header.update('BPMEXT',nextend+i, comment='Extension for Bad Pixel Mask')
           struct.append(hdu)
       nextend+=nsciext


   #update the number of extensions
   saltkey.new('NSCIEXT',nsciext,'Number of science extensions', struct[0])
   saltkey.new('NEXTEND',nextend,'Number of data extensions', struct[0])


   return struct

def CreateVariance(inhdu, sci_ext, var_ext):
   """Create a variance hdu from an input hdu"""
   rdnoise=inhdu.header['RDNOISE']
   gain=inhdu.header['GAIN']
   #create the variance array
   data=inhdu.data.copy()
   if (data <=0).any():
       j=np.where(data>0)
       min_pos=data[j].min()
       j=np.where(data<=0)
       data[j]=min_pos
   data=(data+(rdnoise/gain)**2)

   header=inhdu.header.copy()
   header['EXTVER'] = var_ext
   header['SCIEXT'] = (sci_ext,'Extension of science frame')
   return fits.ImageHDU(data=data, header=header, name='VAR')

def createbadpixel(inhdu, bphdu, sci_ext, bp_ext):
   """Create the bad pixel hdu from a bad pixel hdu"""
   if bphdu is None:
       data=inhdu[sci_ext].data*0.0
   else:
       try:
           infile=inhdu._HDUList__file.name
           bpfile=bphdu._HDUList__file.name
       except:
           infile=inhdu.filename
           bpfile=bphdu.filename
    
       
       if not saltkey.compare(inhdu[0], bphdu[0], 'INSTRUME', infile, bpfile):
           message = '%s and %s are not the same %s' % (infile,bpfile, 'INSTRUME')
           raise SaltError(message)
       for k in ['CCDSUM', 'NAXIS1', 'NAXIS2']:
           if not saltkey.compare(inhdu[sci_ext], bphdu[sci_ext], k, infile, bpfile):
                  message = '%s and %s are not the same %s' % (infile,bpfile, k)
                  raise SaltError(message)
       data=bphdu[sci_ext].data.copy()

   header=inhdu[sci_ext].header.copy()
   header['EXTVER']=bp_ext
   header['SCIEXT']=(sci_ext,'Extension of science frame')
   return fits.ImageHDU(data=data, header=header, name='BPM')



# -----------------------------------------------------------
# main code
if not iraf.deftask('saltprepare'):
   parfile = iraf.osfn("saltred$saltprepare.par")
   t = iraf.IrafTaskFactory(taskname="saltprepare",value=parfile,function=saltprepare, pkgname='saltred')
