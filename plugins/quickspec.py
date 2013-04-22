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
QUICKSPEC -- QUICKSPEC  provides a plugin for saltfirst that provides quick 
reductions for spectroscopic data

Author                 Version      Date
-----------------------------------------------
S M Crawford (SAAO)    0.1           5 Jun 2011

"""
import os 
import pyfits

#pysalt imports
from pyraf import iraf
from pyraf.iraf import pysalt
from pyraf.iraf import saltspec

from specreduce import specreduce

import saltsafeio


def quickspec(profile, lampid=None, findobj=False, objsection=None, skysection=None, clobber=False, logfile='saltclean.log', verbose=True):
   """From mosaicked data, produce wavelength calibrated files"""
   struct=pyfits.open(profile)
   xlen=len(struct[1].data[0])
   ylen=len(struct[1].data)
   print xlen,ylen
   badpixelimage=None
   caltype='rss'
   function='polynomial'
   order=3
   skysub=True
   if saltsafeio.checkfornone(lampid): skysub=False

   if objsection is None:
      y1=0.5*ylen-0.05*ylen
      y2=0.5*ylen+0.05*ylen
      objsection='[%i:%i]' % (y1,y2)
   else:
      y1=int(objsection[1:].split(':')[0])
      y2=int(objsection[:-1].split(':')[1])

   if skysection is None:
      skysection='%i:%i' % (y2+0.1*ylen,y2+0.2*ylen)
   thresh=3.0
   logfile='saltspec.log'
   specreduce(images=profile,badpixelimage=badpixelimage,caltype=caltype,
                function=function, order=order, skysub=skysub, skysection=skysection,
                findobj=findobj, objsection=objsection, thresh=thresh,
                clobber=clobber, logfile=logfile, verbose=verbose)

   return y1,y2

def quickap(profile, objsection=None, logfile='salt.log', verbose=True):
   "run a quick aperture extraction for the data"
   outfile=profile[:-4]+'txt'
   print outfile
   saltspec.specextract(profile, outfile, method='normal', section=objsection, thresh=3.0,    \
                minsize=3.0, outformat='ascii', convert=True, clobber=True,    \
                logfile=logfile, verbose=verbose)
    
