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
AutoIDENTIFY  is a program to automatically identify spectral lines in 
an arc image.  

Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)    1.0       21 Aug 2010

TODO
----


LIMITATIONS
-----------

"""
# Ensure python 2.5 compatibility

import os, sys
import time
import numpy as np

from pyraf import iraf 
import saltprint, saltio, saltkey
import saltsafekey
import saltsafeio
from saltsafelog import logging
from salterror import SaltError, SaltIOError


from PySpectrograph.Spectra import apext, detectlines


import spectools as st
from spectools import SALTSpecError

debug=True

#autoidentify_options=['Zeropoint', 'Matchlines', 'MatchZero', 'FullXCor']
autoidentify_options=['Zeropoint', 'Matchlines', 'MatchZero']

def AutoIdentify(xarr, specarr, slines, sfluxes, ws, method='Zeropoint',   \
                     rstep=1, istart=None, nrows=1, res=2, dres=0.1,         \
                     sigma=5, niter=5, mdiff=20, dc=20, ndstep=20, farr=None, \
                     oneline=False, verbose=True):
   """Automatically find the wavlength solution for the entire image.  The following 
      methods are used:

      Zeropoint--Assume that the form for the initial guess of the wavelength solution
                 is correct

      Matchlines--For each line, use the initial guess to match the lines and then find
                  the best fit

      MatchZero--First calculate the zeropoint, then match the lines.  

      FullXCor--Provides a full cross correlation for all the coefficients



   """
   ImageSolution={}

   #run it if only the zeropoint needs to be calculated
   if method=='Zeropoint':
      func=st.findzeropoint
      ImageSolution=runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=False, oneline=oneline, \
                  rstep=rstep, istart=istart, nrows=nrows, res=res, dres=dres, farr=farr,      \
                  dsigma=sigma, dniter=niter, verbose=verbose, dc=dc, ndstep=ndstep)

   #use a line matching algorithm to match the lines
   #in the image with those in the line list
   if method=='Matchlines':
      func=st.findwavelengthsolution
      #set wdiff
      try:
           wdiff=mdiff*ws.coef[1]
      except:
           wdiff=mdiff
      ImageSolution=runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=True, oneline=oneline,\
                  rstep=rstep, istart=istart, nrows=nrows, res=res, dres=dres, farr=farr,      \
                  dsigma=sigma, dniter=niter, verbose=verbose, mdiff=mdiff, wdiff=wdiff, sigma=sigma, niter=niter)

   #first fit a zeropoint, then match the lines, and then 
   #find the rest of the points by using only the zeropoint
   if method=='MatchZero':
      func=st.findzeropoint
      ws=runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=False, oneline=True, \
                  rstep=rstep, istart=istart, nrows=nrows, res=res, dres=dres, farr=farr,      \
                  dsigma=sigma, dniter=niter, verbose=verbose, dc=10, ndstep=20)

      func=st.findwavelengthsolution
      ws=runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=True, oneline=True,\
                  rstep=rstep, istart=istart, nrows=nrows, res=res, dres=dres, farr=farr,      \
                  dsigma=sigma, dniter=niter, verbose=verbose, sigma=sigma, niter=niter)
      func=st.findzeropoint
      ImageSolution=runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=False, oneline=oneline,\
                  rstep=rstep, istart=istart, nrows=nrows, res=res, dres=dres, farr=farr,      \
                  dsigma=sigma, dniter=niter, verbose=verbose, dc=dc, ndstep=ndstep)


   if method=='FullXcor':
      func=st.findxcor 
      dcoef=ws.coef*0.1
      dcoef[0]=dc
      ws=runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=True, oneline=True,\
                  rstep=rstep, istart=istart, nrows=nrows, res=res, dres=dres, farr=farr,      \
                  dsigma=sigma, dniter=niter, verbose=verbose, dcoef=dcoef)
      ImageSolution=runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=False, oneline=oneline, \
                  rstep=rstep, istart=istart, nrows=nrows, res=res, dres=dres, farr=farr,     \
                  dsigma=sigma, dniter=niter, verbose=verbose, dc=dc, ndstep=ndstep)



   return ImageSolution
   
def runsolution(xarr, specarr, slines, sfluxes, ws,  func, ivar=None,          \
                fline=True, oneline=False, farr=None, rstep=20,\
                istart=None, nrows=1, dsigma=5, dniter=5, res=2.0, dres=0.1, verbose=True, **kwargs):
   """Starting in the middle of the image, it will determine the solution
      by working its way out to either edge and compiling all the results into 
      ImageSolution

      xarr--Full range in x of pixels to solve for
   
      specarr--Input 2D flux

      func--function to use for the solution
 
      fline--whether the spectral lines are in list or array format

      oneline--whether to measure one line or all lines


   """
   #set up the variables
   ImageSolution={}

   #Setup the central line if it isn't specified
   if istart==None: istart=int(0.5*len(specarr))

   #set up the flux from the central line (or the line specified by the user in istart)
   if farr is None:
       specext=apext.apext(xarr, specarr, ivar=ivar)
       farr=apext.makeflat(specarr, istart, istart+nrows)
       farr=st.flatspectrum(xarr, farr, mode='poly', order=2)

   #detect the lines 
   cxp=detectlines.detectlines(xarr, farr, dsigma, dniter)
   nlines=len(cxp)

       
   #first set up the artificial spectrum
   swarr, sfarr=st.makeartificial(slines, sfluxes, farr.max(), res, dres)
   
   #find the solution for the central wavelegnth
   k=istart
   min_lines=0.1*len(cxp)
   if oneline:
     if fline:
       mws=solution(xarr, farr, slines, sfluxes, ws, func, \
                    min_lines=min_lines, dsigma=dsigma, dniter=dniter, **kwargs)
     else:
       mws=solution(xarr, farr, swarr, sfarr, ws, func, \
                    min_lines=min_lines,dsigma=dsigma, dniter=dniter, **kwargs)
     return mws
   else:
     ImageSolution[k]=ws

   #now loop through each step, and calculate the wavelengths for the given 
   for i in range(rstep,int(0.5*len(specarr)), rstep):
       for k in [istart-i, istart+i]:
           lws=getwsfromIS(k, ImageSolution)
           #set up the flux from the set of lines
           farr=apext.makeflat(specarr, k, k+nrows)
           try:
              farr=st.flatspectrum(xarr, farr, mode='poly', order=2)
           except:
              continue
           if fline:
               fws=solution(xarr, farr, slines, sfluxes, lws, func, \
                            min_lines=min_lines, dsigma=dsigma, dniter=dniter, **kwargs)
           else:
               fws=solution(xarr, farr, swarr, sfarr, lws, func, \
                            min_lines=min_lines, dsigma=dsigma, dniter=dniter, **kwargs)
           if fws is not None:
               ImageSolution[k]=fws

               if verbose:
                   p_new=i*100.0/(0.5*len(specarr))
                   #ctext='Percentage Complete: %d %d %f\r' % (i,p_new, time.clock()) #p_new
                   #sys.stdout.write(ctext)#
                   #sys.stdout.flush()
	           print "%5i %3i %3.2f" % (k, fws.func.mask.sum(), fws.sigma(fws.func.x, fws.func.y) )

   return ImageSolution

def solution(xarr, farr, sl, sf, ws, func, min_lines=2, dsigma=5, dniter=3, pad=50, **kwargs):
   """Extract a single line and calculate the wavelneght solution"""


   #check to see if there are any points
   xp=detectlines.detectlines(xarr, farr, dsigma, dniter)

   #print len(xp)
   if len(xp) > min_lines and ws:
       #make the artificial list
       wmin=ws.value(xarr.min())
       wmax=ws.value(xarr.max())
       smask=(sl>wmin-pad)*(sl<wmax+pad)

       #fit the function
       try:
           fws=func(xarr, farr, sl[smask], sf[smask], ws, **kwargs)
       except SALTSpecError, e:
           return None
       except TypeError, e:
           return None
       except IndexError, e:
           return None
       except Exception, e:
           print e
           return None
           
       return fws

   return None

def getwsfromIS(k, ImageSolution):
    """From the imageSolution dictionary, find the ws which is nearest to the value k

    """
    ISkeys=np.array(ImageSolution.keys())
    ws=ImageSolution[ISkeys[abs(ISkeys-k).argmin()]]
    if ws==None:
       dist=abs(ISkeys[0]-k)
       ws=ImageSolution[ISkeys[0]]
       for i in ISkeys:
           if ImageSolution[i] and abs(i-k)<dist:
               dist=abs(i-k)
               ws=ImageSolution[i]
    return ws
    

    


