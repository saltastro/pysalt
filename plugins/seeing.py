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
CALCSEEING--Determine the seeing gvine an array of magnitudes
and an array of FWHM

Author                 Version      Date
-----------------------------------------------
S M Crawford (SAAO)    0.1          19 Jun 2011

"""

import math
import numpy as np

from saltfit import *

def hist_stats(arr, bins=200, range=(0,50)):
    h=np.histogram(arr, bins=bins, range=range)
    i=h[0][:-1].argmax()
    #histogram returns the lower edge and we want the middle
    meandiff=np.diff(h[1]).mean()
    return i,h[1][:-1]+meandiff, h[0]


def mirror_hist(p_i, X, data):
   """Remove all points higher than the peak, and replace with a mirror image of the lower points"""

   ndata=data.copy()
   for i in range(len(data[:])):
       if X[i]>X[p_i]:
           j=p_i - (i-p_i)
           if j>=0:
               ndata[i]=data[j]
           else:
               ndata[i]=0
   return ndata


def seeing_stats(arr, bins=200, hrange=(0,50)):
   """Fit_stats finds the peak in a distribution, mirrors the 
   data around that peak rejecting all values above the peak, and
   then returns the best Gaussian fit to that peak

   arr--array of FWHM
   bins--the number of bins for the histogram
   hrange--the range for the histogram
  
   """
   #create the histogram
   p_i, X, data=hist_stats(arr)
   #remove all the seeing values 'above' the peak
   ndata=mirror_hist(p_i, X, data)
   #set up the parameters to fit
   m = sum(X*ndata)/sum(ndata)
   width = (abs(sum((X-m)**2*ndata)/sum(ndata)))**0.5
   dmax=data.max()
   return m, width, dmax, X[p_i]

   #this is unreliable
   mu = Parameter(m)
   sigma = Parameter(width)
   height = Parameter(dmax)
   #return   m, width, max, X[p_i]

   #fit the data
   def f(x): return height() * np.exp(-0.5*((x-mu())/sigma())**2)

   info=fit(f, [mu, sigma, height], ndata)
   return  mu(), sigma(), height(), X[p_i]


def airmass(Z, units='degrees'):
    """Calculate the airmass for a source according to Airmass=sec(Z) where
       Z is the zenith distnce
 
       Z--zenith distance
       units--units of Z
    """
    if units=='degrees': Z=math.radians(Z)
    return 1.0/math.cos(Z)
