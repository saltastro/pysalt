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

"""SALTFPSKYRING is a tool to remove skylines from images. 

Updates:

20130115
    * Updated the code to work in python only

TODO:
    * Account for big hulking galaxies in the middle of the chip
    * Account of not subtracting off where there is no data

"""

# Ensure python 2.5 compatibility
from __future__ import with_statement

import os
import sys

import numpy as np
import numpy.ma as ma
import pyfits

from pyraf import iraf
from pyraf.iraf import pysalt
import saltsafekey
import saltsafeio as saltio
import fpsafeio
from saltsafelog import logging

debug=True

def saltfpskyring(images,outimages,outpref, axc, ayc, arad, rxc, ryc, pmin, pmax, swindow=5,
                  clobber=False, logfile='salt.log', verbose=True):
    """Sky subtracts Fabry-Perot images"""


# start log now that all parameter are set up          

    with logging(logfile, debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       #verify that the input and output lists are the same length
       saltio.comparelists(infiles,outfiles,'Input','output')

       for img, oimg in zip(infiles, outfiles):

           #open up the file
           hdu=saltio.openfits(img)

           #determine the azimuthally averaged profile
           rpro, radial_profile=median_radial_profile(hdu[0].data, xc=axc, yc=ayc, rmax=arad, 
                      nbins=100, pmin=pmin, pmax=pmax)

           if swindow>1:
               radial_profile=np.convolve(radial_profile, np.ones(swindow), mode='same')/swindow

           # calculate the indices from the image
           y, x = np.indices(hdu[0].data.shape)
           radius = np.hypot(x - axc, y - ayc)
           
           #subtract off the sky data
           mdata = hdu[0].data-np.interp(radius, rpro, radial_profile)
           hdu[0].data=mdata

           # write FITS file
           saltio.writefits(hdu,oimg, clobber=clobber)
           saltio.closefits(hdu)

           message = 'SALTFPSKYRING -- Subtracted sky from  %s' % (img )
           log.message(message, with_header=False)



           

def median_radial_profile(data, xc=0, yc=0, rmax=500, nbins=50, pmin=0, pmax=100):
    """
    calculate the azimuthally median radial profile.

    data - 2D image
    xc - x center
    yc - y center
    rmax - maximum radius for the distribution
    nbins - number of bins to use
    maskval - threshold value for including data in the profile
    """

    # calculate the indices from the image
    y, x = np.indices(data.shape)
    radius = np.hypot(x - xc, y - yc)

    #create the radii to sample out to
    hbin=0.5*float(rmax)/nbins
    rpro=np.arange(0,rmax,float(rmax)/nbins)+hbin
    arr=np.zeros(len(rpro))

    #step through each radius and determine the median value for the bin
    for i,r in enumerate(rpro):
        mask=(radius>r-hbin)*(radius<r+hbin)*(data>pmin)*(data<pmax)
        arr[i]=np.median(data[mask])
    return rpro, arr 


    
# -----------------------------------------------------------
# main code

parfile = iraf.osfn("saltfp$saltfpskyring.par")
t = iraf.IrafTaskFactory(taskname="saltfpskyring",value=parfile,function=saltfpskyring,pkgname='saltfp')
