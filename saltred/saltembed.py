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
SALTEMBED--For data that has been taken in windowed mode, re-combine
the data into the number of amplifiers that it should be.  This 
tasks will just change windowed data which will be on N*Namps extensions
into data that only has Namps extensions where N is the number of
windows in the data.

Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)  1.0          25 Jul 2011

TODO
-----------------------------------------------
-The number of amplifiers is hardwired

UPDATES
-----------------------------------------------


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

def saltembed(images,outimages,outpref, clobber=True,logfile='salt.log',verbose=True): 

   with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       #verify that the input and output lists are the same length
       saltio.comparelists(infiles,outfiles,'Input','output')

       # open each raw image file
       for img, oimg, in zip(infiles, outfiles):

           #open the fits file
           struct=saltio.openfits(img)

           #get the number of CCDs
           nccd=saltkey.get('NCCDS', struct[0])

           #set the number of amps per CCD--*TODO* Read from header
           namps=2

           #get the number of windows--try to read from header
           try:
              nwindows=saltkey.get('NWINDOWS', struct[0])
           except:
              nwindows=len(struct)/(nccd*namps)
 
           outstruct=embedimage(struct, nccd, namps, nwindows)

           saltio.writefits(outstruct, oimg, clobber=clobber)

           message = 'SALTEMBED -- %s => %s' % (img, oimg)
           log.message(message, with_header=False)

def embedimage(struct, nccd=2, namps=2, nwindows=1):
    """For given number of windows in struct, create an embeded image.   

    struct: an image structure

    nccd: number of ccds
  
    namps: number of amplifiers per ccd

    nwindows:  number of windows

    """

    #create the output structure
    outstruct=fits.HDUList(struct[0])

    #get the total size of the section
    decsect=saltio.getSection(saltkey.get('DETSIZE', struct[0]), iraf_format=True)
    xbin,ybin=saltkey.get('CCDSUM', struct[0]).strip().split()
    xbin=int(xbin)
    ybin=int(ybin)

    #determine the data size--this assumes that the CCDs are in the x direction and 
    #all of the CCDs are the same size
    xsize=decsect[3]/xbin/(nccd*namps)
    ysize=decsect[1]/ybin
    edata=np.zeros((ysize,xsize))
  
    #loop through the number of extensions
    for i in range(1,nccd*namps+1):

        #create the total size for the extension
        extdata=edata.copy()


        #loop throug each window and add it to the extension
        xmax=0
        for j in range(nwindows):
           winsect=saltio.getSection(saltkey.get('AMPSEC', struct[i+6*j]), iraf_format=True)
           y1=(winsect[0]-1)/ybin
           #y2=(winsect[1])/ybin+1
           y2=y1+len(struct[i+6*j].data)

           x1=(winsect[2]-1)/xbin
           x2=len(struct[i+6*j].data[0])
           #x2=(winsect[3])/xbin
         
           extdata[y1:y2,x1:x2]=struct[i+6*j].data
        
           xmax=max(xmax, x2)
           
 
        #add the correct header information *TODO*   
 
        #append to the hdu list
        exthdu=yfits.ImageHDU(extdata[:,0:xmax])
        outstruct.append(exthdu)


    #return the new extension
    return outstruct

# -----------------------------------------------------------
# main code

if not iraf.deftask('saltembed'):
    parfile = iraf.osfn("saltred$saltembed.par")
    t =  iraf.IrafTaskFactory(taskname="saltembed",value=parfile,function=saltembed, pkgname='saltred')
