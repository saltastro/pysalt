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

"""SALTFPMASK is a tool to aperture mask Fabry-Perot images. Using T. Williams code.  The code assumes all the files are in the same directory. Also assumes that if there is a config file, it is also in the same directory as the data. Note that this config file is in the original FORTRAN code format so that the user does not have to write another file.

Updates:

20100615
    * First wrote the code

20110428
    * Removed the fortran and wrote in pure python

TODO:
* Replace the mean with the biweight estimator
"""

# Ensure python 2.5 compatibility
from __future__ import with_statement

import os
import sys
import numpy as np
#import pyfits

from pyraf import iraf
from pyraf.iraf import pysalt
import saltsafekey as saltkey
import saltsafeio as saltio
import fpsafeio
from saltsafelog import logging

# This reads the FORTRAN config file if it exists

from fortranfp import masking_wrapper
from fortranfp.masking_wrapper import getpfp

debug=True

def saltfpmask(images, outimages, outpref, axc,ayc, arad, maskmethod='c', maskvalue=0, \
               radi=None,rado=None,clobber=False, logfile='saltfp.log', verbose=True):  
    """Aperture masks Fabry-Perot images"""

    with logging(logfile, debug) as log:
       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       #verify that the input and output lists are the same length
       saltio.comparelists(infiles,outfiles,'Input','output')

       for img, oimg in zip(infiles, outfiles):
           #open the image
           hdu=saltio.openfits(img)

           #test that the data are in the first 
           try:
               data=hdu[0].data
           except Exception, e:
               message='SALTFPMASK--ERROR:  Could not access data in Primary exention of %s because %s' % (img, e)
               log.message(message)
          

           #determine the mask value if it is a region
           if maskmethod=='region': 
               maskvalue=calc_maskvalue(data, axc, ayc, radi, rado)

           #mask the image
           message="SALTFPMASK--Masking image %s with a value of %s" % (img, maskvalue)
           log.message(message, with_header=False)
           hdu[0].data= maskimage(data, axc, ayc, arad, maskvalue)
           try:
               pass
           except Exception, e:
               message='SALTFPMASK--ERROR:  Could not create mask for %s because %s' % (img, e)
               log.message(message)

           #write out the image
           saltio.writefits(hdu, oimg, clobber=clobber)

def calc_maskvalue(data, axc, ayc, radi, rado):
    y,x=np.indices((len(data), len(data[0])))
    r = ((x-axc)**2+(y-ayc)**2)**0.5
    rmask=(r>radi)*(r<rado)
    return data[rmask].mean()

def maskimage(data,  axc, ayc, arad, maskvalue):
    """provide a constant mask beyond a region of the image
    """
    y,x=np.indices((len(data), len(data[0])))
    r = ((x-axc)**2+(y-ayc)**2)**0.5
    rmask=(r>arad)
    data[rmask]=maskvalue 
    return data
    
    

# -----------------------------------------------------------
# main code

parfile = iraf.osfn("saltfp$saltfpmask.par")
t = iraf.IrafTaskFactory(taskname="saltfpmask",value=parfile,function=saltfpmask,pkgname='saltfp')
