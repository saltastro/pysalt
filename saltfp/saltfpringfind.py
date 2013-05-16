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
SALTFPRINGFIND is a program to find rings from a Fabry Perot image.
The first step is to identify the approximate position, center,
radius, and FWHM of a ring.  

# Author                 Version      Date
# -----------------------------------------------
# S. M. Crawford (SAAO)    0.1       21 Jul 2011

"""

from __future__ import with_statement

import os
import pyfits
import math
import numpy as np

from pyraf import iraf 
from pyraf.iraf import pysalt
import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging
from saltstat import iterstat


from display import display
from PySpectrograph.Spectra import findobj

from salterror import SaltError

from fptools import findrings, findcenter

import pylab as pl

debug=True

# -----------------------------------------------------------
# core routine


def saltfpringfind(images, method=None, section=None, thresh=5, minszie=10, niter=5, conv=0.05,
                displayimage=True, clobber=False, logfile='salt.log',verbose=True):

   with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       #check the method
       method=saltio.checkfornone(method)

       # read in the section
       section=saltio.checkfornone(section)
       if section is None: 
           pass
       else:
           section=saltio.getSection(section)
 
       # open each raw image file
       for img in infiles:

          #open the fits file
	  struct=saltio.openfits(img)
          data=struct[0].data

          #determine the background value for the image
          if section is None:
             #if section is none, just use all pixels greater than zero
             bdata=data[data>0]
          else:
             y1,y2,x1,x2=section
             bdata=data[y1:y2,x1:x2]
          bmean, bmedian, bstd=iterstat(bdata, sig=thresh, niter=niter, verbose=False)
          message="Image Background Statistics\n%30s %6s %8s %8s\n%30s %5.4f %5.4f %5.4f\n" %  \
                ('Image', 'Mean', 'Median', 'Std',img, bmean, bmedian, bstd)
          log.message(message, with_stdout=verbose)
    
          mdata=data*(data-bmean>thresh*bstd)

          #prepare the first guess for the image
          ring_list=findrings(data, thresh=5, niter=5, minsize=10)

          #if specified, find the center of the ring
          if method is not None:
             for i in range(len(ring_list)):
                 ring_list[i]=findcenter(data, ring_list[i], method)
            

          #if one peak: no rings.  If two peaks: one ring, if two peaks: four rings
          if len(ring_list)==1:
             msg="One ring dected in image"
          else:
             msg="%i rings found in image" % len(ring_list)
          log.message(message, with_stdout=verbose)

          if displayimage:
             regfile=img.replace('.fits', '.reg')
             if clobber and os.path.isfile(regfile): fout=saltio.delete(regfile)
             fout=open(regfile, 'w')
             fout.write("""# Region file format: DS9 version 4.1
# Filename:  %s
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
physical
""" % img)
             for ring in ring_list:
                 fout.write('circle(%f, %f, %f)\n' % (ring.xc,ring.yc,ring.prad))
                 fout.write('circle(%f, %f, %f)\n' % (ring.xc,ring.yc,ring.prad-3*ring.sigma))
                 fout.write('circle(%f, %f, %f)\n' % (ring.xc,ring.yc,ring.prad+3*ring.sigma))

             fout.close()
             display(img, catname=regfile, rformat='reg')

	  message = 'Ring Parameters\n%30s %6s %6s %6s\n' % ('Image', 'XC', 'YC', 'Radius')
          log.message(message, with_stdout=verbose)
          for ring in ring_list:
              msg='%30s %6.2f %6.2f %6.2f\n' % (img, ring.xc, ring.yc, ring.prad)
              log.message(msg, with_header=False, with_stdout=verbose)

parfile = iraf.osfn("saltfp$saltfpringfind.par")
t = iraf.IrafTaskFactory(taskname="saltfpringfind",value=parfile,function=saltfpringfind,pkgname='saltfp')

