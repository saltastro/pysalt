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

"""SALTFPCALRING is a tool to analyse one or more calibration ring images 
   and determine the radius and centre coordinates for each ring. The code 
   will then create an output file that can be used to calculate the 
   coefficients used for converting pixel position into wavelength for 
   SALT Fabry Perot data.

Updates:

20100621
    * First wrote the code
20130115
    * Updated the code to only used python.
    * Change the code so that the ring was found and the center was determined
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
import salttime
from saltsafelog import logging

from saltfpringfind import findrings

debug=True

def saltfpcalring(images,outfile, waves, thresh=5, niter=3, minsize=10, 
                  clobber=False,logfile='salt.log',verbose=True):  
   """Fits rings in Fabry-Perot ring images"""
   with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # if the outfile exists and not clobber, fail
       saltio.overwrite(outfile,clobber) 

       #open the output file
       fout = saltio.openascii(outfile, 'w')

       #make sure that the list of waves is convertable to a numpy array
       #convert to floats in case the wave is given as a string
       if isinstance(waves, str):  
          waves=[float(x) for x in waves.split(',')]

       try: 
           waves=np.array(waves)
       except:
           raise SaltError('%s is not convertable to a numpy array' % waves)

       #setup the output file
       fout.write('#Comment\n')
       fout.write('# radius   err     xc      yc      z      ut     wave   dn  file\n')

       # open each image and detect the ring
       for img,w  in zip(infiles, waves):

          #measure the ring in each file
          xc, yc, radius, err, z, ut=make_calring(img, thresh=5, niter=3, minsize=10)

          #output the results
          outstr=' %7.2f %6.2f %7.2f %7.2f %6.2f %7.4f %8.3f  0 %s\n' % (radius, err, xc, yc, z, ut, w, img)
          fout.write(outstr)
          log.message(outstr.strip(), with_stdout=verbose, with_header=False)

       fout.close()

def make_calring(image, thresh=5, niter=3, minsize=10):
   """Open each image and measure the position of the ring including its center and radius
     
      Return the information about the calibration ring
   """   

   #open the image
   hdu=saltio.openfits(image)
    
   #extract the time and convert to decimal hours
   utctime=saltkey.get( 'UTC-OBS', hdu[0])
   utctime=salttime.time_obs2hr((utctime.split()[-1]))
  
   #determine the correct etalon and information to extract
   etstate=saltkey.get( 'ET-STATE', hdu[0])

   if etstate.count('S2'):
        etz=saltkey.get('ET1Z', hdu[0])
   elif etstate.count('S3'):
        etz=saltkey.get('ET2Z', hdu[0])
   else:
       msg='This etalon state is not currently supported'
       raise SaltError(msg)

   #extract the ring
   ring_list=findrings(hdu[0].data, thresh=thresh, niter=niter, minsize=minsize)
  
   #TODO: Fit the data for a better deterimination of the central ring parameters
   #iterate to improve the result
 
   #assumes only one ring in the data set
   ring=ring_list[0]

   return ring.xc, ring.yc, ring.prad, 1.0, etz, utctime
  

# -----------------------------------------------------------
# main code

parfile = iraf.osfn("saltfp$saltfpcalring.par")
t = iraf.IrafTaskFactory(taskname="saltfpcalring",value=parfile,function=saltfpcalring,pkgname='saltfp')
