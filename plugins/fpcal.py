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
# DAMAGES (INCte: 2007/05/26
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)    #
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      #
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN #
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE          #
# POSSIBILITY OF SUCH DAMAGE.                                              #
############################################################################


#!/usr/bin/env python

"""
FPCAL-- FPCAL further reduces and quickly calibrates FP data

Author                 Version      Date
-----------------------------------------------
S M Crawford (SAAO)    0.1          20 Jul 2011
ERC (SAAO)             0.2          20 Sep 2013

"""
import os 
import pyfits

#pysalt imports
from pyraf import iraf
from pyraf.iraf import pysalt

from saltflat import saltflat
from saltfpprep import saltfpprep
from saltfpringfit import saltfpringfit


import saltsafeio


def fpcal(profile, flatimage=None, minflat=15000, bthresh=5, niter=5, displayimage=False, clobber=False, logfile='saltclean.log', verbose=True):
   """Flatfield and determine the ring center for SALT FP calibration data"""

   #run saltflat
   if flatimage:
       saltflat(profile, outimages='', outpref='p', flatimage=flatimage,minflat=minflat, clobber=clobber, logfile=logfile, verbose=verbose)

   #run saltfpprep
   saltfpprep('p'+profile, outimages='', outpref='f', full_reduce=False, clobber=clobber, logfile=logfile, verbose=verbose)

   #run saltfpringfit
   section='[600:800,300:500]' #*TODO* Update so it automatically determines this
   outfile='fpout.txt'
   saltfpringfit('fp'+profile, outfile, section=section, bthresh=bthresh, niter=niter,
                displayimage=displayimage, clobber=clobber,logfile=logfile, verbose=verbose)

