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
RSSINFO -- RSSINFO is a plugin for saltfirst that returns infomation
about the spectrograph including the central wavelength, grating information,
resolution, and wavelength coverage.

Author                 Version      Date
-----------------------------------------------
S M Crawford (SAAO)    0.1          19 Jun 2011

"""

from PySpectrograph.Models import RSSModel
from spectools import getslitsize


def rssinfo(grating, grang, arang, slitname, xbin, ybin):
   """RSSINFO--Given informtion about the set up of the spectrograph,
      return information about the spectrograph
    
      grating--name of the grating used
      gratang--angle of the grating in degrees
      arang--articulation ange in degrees
      slitname--name of the slit being used
      xbin--binning in x-direction
      ybin--binning in y-direction

      returns central wavelength, blue wavelength, red wavelenth, resolution
              resolution element
   """

   #get the slitsize
   slit=getslitsize(slitname)

   #set up the rss model
   rss=RSSModel.RSSModel(grating_name=grating, gratang=grang, \
                         camang=arang,slit=slit, xbin=xbin, ybin=ybin)

   #calcuation the resolution element
   res=1e7*rss.calc_resolelement(rss.alpha(), -rss.beta())

   #calculate the central wavelength
   wcen=1e7*rss.calc_centralwavelength()

   #calculate the wavelength edges
   w2=1e7*rss.calc_bluewavelength()
   w1=1e7*rss.calc_redwavelength()
   print w1,w2

   #calculate the 
   R=rss.calc_resolution(wcen/1e7, rss.alpha(), rss.beta())
 
   return wcen, w1, w2, res, R, slit
