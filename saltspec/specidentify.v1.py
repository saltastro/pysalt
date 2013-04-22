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
SPECIDENTIFY  is a program to read in SALT RSS spectroscopic arc lamps and 
determine the wavelength solution for that data.  The input data should be
a SALT arc lamp and a line list or another arc lamp image with high quality
wavelength solution. The lamp list can be either wavelengths, wavelengths
and fluxes, or an arc image with a high quality solution.  The line lamp
can also be left unspecified and the user will manually enter the data.

From there, the user has several different possible choices.  They can provide
a first guess of the coefficients for the wavelength solution or transformation, 
indicate the model for the spectragraph for the first guess, or provide an 
image with a solution already as the first guess.  

Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)    1.0       10 Oct 2009

TODO
----


LIMITATIONS
-----------
1. Currently assumes that the linelist is of the form of an ascii file with 
   either appropriate information in either one or two columns

"""
# Ensure python 2.5 compatibility
from __future__ import with_statement

import os, string, sys, glob, pyfits, math, time
import numpy as np

from pyraf import iraf 
import saltprint, saltio, saltkey
import saltsafekey
import saltsafeio
from saltsafelog import logging
from salterror import SaltError, SaltIOError


from PySpectrograph import RSSModel
from PySpectrograph import apext
from PySpectrograph import WavelengthSolution
from PySpectrograph.detectlines import detectlines


import spectools as st
from spectools import SALTSpecError
from InterIdentify import InterIdentify

debug=True



# -----------------------------------------------------------
# core routine

def specidentify(images,linelist, outfile, guesstype, guessfile, function, \
                 order, rstep, interact, clobber,logfile,verbose,status):

   with logging(logfile,debug) as log:

       #set up the variables
       infiles = []
       outfiles = []
       status = 0

       # Check the input images 
       infiles = saltsafeio.argunpack ('Input',images)

       # create list of output files 
       outfiles = saltsafeio.argunpack ('Input',outfile)

       #if outfiles is a single image, turn it into a list 
       #of the same length as infiles
       if len(outfiles)!=len(infiles):
          if len(outfiles)==1:
               outfiles=outfiles*len(infiles)
          elif len(outfiles)==0:
               outfiles=[None]*len(infiles)
          else:
               msg='Please enter an appropriate number of outfiles'
               raise SALTSpecError(msg)

       # open the line lists
       slines, sfluxes = readlinelist(linelist)

       # Identify the lines in each file
       for img, oimg in zip(infiles, outfiles):
           identify(img, oimg, slines, sfluxes, guesstype, guessfile, function, 
                    order, rstep, interact, clobber, verbose)



#------------------------------------------------------------------
# Find the solution for lines in a file

def identify(img, oimg, slines, sfluxes, guesstype, guessfile, function, order,  
             rstep, interact, clobber, verbose):
   """For a given image, find the solution for each row in the file.  Use the appropriate first guess and 
      guess type along with the appropriate function and order for the fit.

      Write out the new image with the solution in the headers and/or as a table in the multi-extension
      fits file

       returns the status
   """
   status=0
   ImageSolution={}

   #Open up the image
   hdu=saltsafeio.openfits(img)   

   #Read in important keywords

   #determine the central row and read it in
   try:
       data=hdu[1].data
       midline=int(0.5*len(data))
       xarr=np.arange(len(data[midline]))
       specarr=st.flatspectrum(xarr, data[midline, :], mode='poly', order=2)
   except Exception, e:
       message = 'Unable to read in data array in %s because %s' % (img, e)
       raise SALTSpecError(message)


   #determine the type of first guess.  Assumes none 
   if guesstype=='user':
       pass
   elif guesstype=='rss':
       dateobs=saltsafekey.get('DATE-OBS', hdu[0], img)
       utctime=saltsafekey.get('UTC-OBS', hdu[0], img)
       instrume=saltsafekey.get('INSTRUME', hdu[0], img)
       grating=saltsafekey.get('GRATING', hdu[0], img)
       grang=saltsafekey.get('GR-ANGLE', hdu[0], img)
       arang=saltsafekey.get('AR-ANGLE', hdu[0], img)
       filter=saltsafekey.get('FILTER', hdu[0], img)
       xbin, ybin = saltsafekey.ccdbin( hdu[0], img)
       if not instrume in ['PFIS', 'RSS']:
           msg='%s is not a currently supported instrument' % instrume
           raise SALTSpecError(msg)
       ws=useRSSModel(xarr, grating, grang, arang, xbin, function=function, order=order)
   elif guesstype=='image':
      pass
   else:
      ws=None

   #if interactive is selected launch the tool in interactive mode
   #if not, produce a wavelength solution from the information already
   #provided.
   xdiff=20
   ws= findwavelengthsolution(xarr, specarr, slines, sfluxes, ws, xdiff, function, 
                              order, zeropoint=False, interact=interact, verbose=verbose)
   ws=findwavelengthsolution(xarr, specarr, slines, sfluxes, ws, xdiff, 
                             function, order, zeropoint=True, dc=3, nstep=100,
                             interact=interact, verbose=verbose)
   ImageSolution[midline]=ws

   #now using that information, find the fit to all other data
   #repeat in the same way if interactive is selected
   #start in the middle and then go either way
   uws=WavelengthSolution.WavelengthSolution(ws.x_arr, ws.w_arr, order=order, function=function)
   uws.fit()
   lws=WavelengthSolution.WavelengthSolution(ws.x_arr, ws.w_arr, order=order, function=function)
   lws.fit()
   zpfind=True 
   dcstep=3
   nstep=50 
   for i in range(rstep,int(0.5*len(data)), rstep):
   #for i in range(1,10, rstep):
           specarr=st.flatspectrum(xarr, data[midline-i, :], mode='poly', order=2)
           lws=findwavelengthsolution(xarr, specarr, slines, sfluxes, lws, xdiff, 
                                      function, order, zeropoint=zpfind, dc=dcstep, nstep=nstep,
                                      interact=interact, verbose=verbose)
           ImageSolution[midline-i]=lws
           specarr=st.flatspectrum(xarr, data[i+midline, :], mode='poly', order=2)
           uws=findwavelengthsolution(xarr, specarr, slines, sfluxes, uws, xdiff, 
                                      function, order, zeropoint=zpfind, dc=dcstep, nstep=nstep,
                                      interact=interact, verbose=verbose)
           ImageSolution[midline+i]=uws
           if verbose:
              p_new=i*100.0/(0.5*len(data))
              ctext='Percentage Complete: %d %d %f\r' % (i,p_new, time.clock()) #p_new
              sys.stdout.write(ctext)
              sys.stdout.flush()


   #set up the list of solutions to into an array
   key_arr=np.array(ImageSolution.keys())
   arg_arr=key_arr.argsort()
   ws_arr=np.zeros((len(arg_arr),len(ws.coef)+1), dtype=float)

   #write the solution to an array 
   for j,i in enumerate(arg_arr):
      if  isinstance(ImageSolution[key_arr[i]], WavelengthSolution.WavelengthSolution):
          ws_arr[j,0]=key_arr[i]
          ws_arr[j,1:]=ImageSolution[key_arr[i]].coef

   #write the solution as an file 
   if outfile:
       #write header to the file that should include the order and function
       if os.path.isfile(outfile) and not clobber:
          dout=open(outfile, 'a')
       else:
          dout=open(outfile, 'w')

       msg='#WS: Wavelength solution for image %s\n' % img
       msg+= '#The following parameters were used in determining the solution:\n'
       msg+= '#name=%s\n' % img 
       msg+= '#time-obs=%s %s\n' % (dateobs, utctime )
       msg+= '#instrument=%s\n' % instrume
       msg+= '#grating=%s\n' % grating.strip()
       msg+= '#graang=%s\n' % grang   
       msg+= '#arang=%s\n' % arang 
       msg+= '#filter=%s\n' % filter.strip() 
       msg+= '#Function=%s\n' % function
       msg+= '#Order=%s\n' % order
       msg+= '#Starting Data\n'
       dout.write(msg)

       for i in range(len(ws_arr)):
           if ws_arr[i,0]:
               msg = '%5.2f ' % ws_arr[i,0]
               msg +=' '.join(['%e' % k for k in ws_arr[i,1:]])
               dout.write(msg+'\n')
       dout.write('\n')
       dout.close()
          

   #write the solution as an extension
   hdu.close()

   return 

def findwavelengthsolution(xarr, specarr, slines, sfluxes, ws, xdiff, function, order, 
                           zeropoint=True, dc=10, nstep=20, interact=False, verbose=True):
   """Find the wavlength solution.  Either find it automatically or find it interactively depending on what is selected

      return ws
   """
   if interact:
       iplot=InterIdentify(xarr, specarr, slines, sfluxes, ws, xdiff=xdiff, function=function, \
                           order=order, verbose=True)
       iplot.runplotdata()
       ws=iplot.ws
   else:
       try:
           ws = findsolution(xarr, specarr, slines, sfluxes, ws,  function=function, order=order, dc=dc, nstep=nstep, zeropoint=zeropoint, verbose=verbose)
       except Exception, e:
           print e
           return None

   return ws


def findzeropoint(xarr, specarr, slines, sfluxes, ws, swarr=None, sfarr=None, function='poly', order=3, dc=10, nstep=20, res=2.0, dres=0.1):
   """Uses cross-correlation to find the best fitting zeropoint""" 

   #if an initial solution, then cut the template lines to just be the length of the spectrum
   if ws==None: return ws


   lmax=specarr.max()
   wmin=ws.value(xarr.min())
   wmax=ws.value(xarr.max())
   mask=(slines>wmin)*(slines<wmax)
   sl=slines[mask]
   sf=sfluxes[mask]
   #if no lines after matching, then just exit
   if not sl.any(): return None

   if swarr==None and sfarr==None:
       swarr, sfarr=st.makeartificial(sl, sf, lmax, res, dres)

   #cross-correlate the spectral lines and the observed fluxes in order to refine the solution
   nws=WavelengthSolution.WavelengthSolution(ws.x_arr, ws.w_arr, order=order, function=function)
   nws.setcoef(ws.coef)

   #create the range of coefficents
   dcoef=ws.coef*0.0
   #dcoef[-2]=0.001
   dcoef[-1]=dc
   dlist=st.mod_coef(ws.coef, dcoef, 0, nstep)

   #loop through them and deteremine the best cofficient
   cc_arr=np.zeros(len(dlist), dtype=float)
   zp_arr=np.zeros(len(dlist), dtype=float)
   dp_arr=np.zeros(len(dlist), dtype=float)
   for i in range(len(dlist)):
       #set the coeficient
       nws.setcoef(dlist[i])

       #set the wavelegnth coverage 
       warr=nws.value(xarr)

       #resample the artificial spectrum at the same wavelengths as the 
       asfarr=st.interpolate(warr, swarr, sfarr, left=0.0, right=0.0)

       #calculate the correlation value
       zp_arr[i]=dlist[i][-1]
       dp_arr[i]=dlist[i][-2]
       cc_arr[i]=st.ncor(specarr, asfarr)
   
   #print cc_arr,  zp_arr
   nws.setcoef(dlist[cc_arr.argmax()])
   coef=np.polyfit(zp_arr, cc_arr, 2)
   nws.coef[-1]=-0.5*coef[1]/coef[0]
   #coef=np.polyfit(dp_arr, cc_arr, 2)
   #nws.coef[-2]=-0.5*coef[1]/coef[0]
   
   nws.setcoef(nws.coef)
   warr=nws.value(xarr)
   ws=WavelengthSolution.WavelengthSolution(xarr, warr, order=order, function=function)
   ws.fit()
   return ws

def findsolution(xarr, specarr, slines, sfluxes, ws, function, order, res=2.0, dres=0.1, dc=10, nstep=20, zeropoint=True, verbose=True):
   """Automatically finds the wavelength solution.   """ 

   ws=findzeropoint(xarr, specarr, slines, sfluxes, ws,  function=function, order=order, dc=dc, nstep=nstep, res=res, dres=dres)
   #if poor value return None
   if ws==None: return ws

   #if only the zeropoint needed to be found, return
   if zeropoint:  return ws
 
   #find the best match
   ws = st.findwavelengthsolution(xarr, specarr, slines, sfluxes, ws, verbose, function=function, order=order)

   #return the solution
   return ws

def useRSSModel(xarr,grating, gratang, camang, binx, order=1, function='poly', slit=1.0):
   """Returns the wavelength solution using the RSS model for the spectrograph
     

   """

   #set up the rss model
   rssmodel=RSSModel.RSSModel(grating_name=grating.strip(), gratang=gratang, camang=camang,slit=slit, xbin=binx)
   rss=rssmodel.rss
   
   #now for each position on the detector, calculate the wavelength at that position

   d=rss.detector.xbin*rss.detector.pix_size*(xarr-0.5*len(xarr))
   alpha=rss.gratang
   beta=rss.gratang-rss.camang
   dbeta=-np.degrees(np.arctan(d/rss.camera.focallength))
   y=1e7*rss.calc_wavelength(alpha, beta+dbeta)

   #for these models, calculate the wavelength solution
   ws=WavelengthSolution.WavelengthSolution(xarr, y, order=order, function=function)
   ws.fit()

   return ws


#------------------------------------------------------------------
# Read in the line list file

def readlinelist(linelist):
   """Read in the line lists.  Determine what type of file it is.  The default is
       an ascii file with line and relative intensity.  The other types are just line, 
       or a wavelenght calibrated fits file 

      return lines, fluxes, and status
   """
   slines=[]
   sfluxes=[]
   status=0

   #Check to see if it is a fits file
   #if not, then read in the ascii file
   if linelist[-4:]=='fits':
       try:
           slines, sfluxes=readfitslinelist(linelist)
       except Exception, e:
           message='Unable to read in the line list %s because %s' % (linelist, e)
           raise SALTSpecError(message)
   else:
       try:
           slines, sfluxes=readasciilinelist(linelist)
       except Exception, e:
           message='Unable to read in the line list %s because %s' % (linelist, e)
           raise SALTSpecError(message)

   #conver to numpy arrays   
   try:
       slines=np.asarray(slines)
       sfluxes=np.asarray(sfluxes)
   except Exception, e:
       message='Unable to create numpy arrays because %s' % (e)
       raise SALTSpecError(logfile, message)
 
   return slines, sfluxes

#------------------------------------------------------------------
# Read in the line list file

def readfitslinelist(linelist):
   """Read in the line lists from an fits file.  If it is a 2-D array
      it will assume that it is an image and select the central wavlength

      return lines, fluxes, and status
   """
   slines=[]
   sfluxes=[]
   
   #open the image
   shdu=pyfits.open(linelist)
   nhdu=len(shdu)
   #determine if it is a one or two-d image
   #if ndhu=0 then assume that it is in the zeroth image
   #otherwise assume the data is in the first extension
   #assumes the x-axis is the wavelength axis
   if nhdu==1:
      ctype1=shdu[0].header['CTYPE1']
      crval1=shdu[0].header['CRVAL1']
      cdelt1=shdu[0].header['CDELT1']
      if shdu[0].data.ndim==1:
           data=shdu[0].data
           wave=crval1+cdelt1*np.arange(len(shdu[0].data))
    
   #detect lines in the input spectrum and identify the peaks and peak values
   slines, sfluxes=st.findpoints(wave, data, 3, 5)
   """
   figure(figsize=(8,8), dpi=72)
   axes([0.1, 0.1, 0.8, 0.8])
   plot(wave, data, ls='-')
   plot(slines, sfluxes, ls='', marker='o')
   xlim(4220,4900)
   show()
   """

   return slines, sfluxes


#------------------------------------------------------------------
# Read in the line list file

def readasciilinelist(linelist):
   """Read in the line lists from an ascii file.  It can either be a 
       file with one or two columns.  Only read in lines that are not 
       commented out.

      return lines, fluxes, and status
   """
   slines=[]
   sfluxes=[]

   #read in the file
   f=open(linelist)
   lines=f.readlines()
   f.close()

   #for each line, 
   for l in lines:
        l=l.strip()
        if not (l and l.startswith('#')):
           l=l.split()
           slines.append(float(l[0]))
           try:
              sfluxes.append(float(l[1]))
           except IndexError:
              sfluxes.append(-1)
   return slines, sfluxes

# main code 

parfile = iraf.osfn("saltspec$specidentify.par") 
t = iraf.IrafTaskFactory(taskname="specidentify",value=parfile,function=specidentify, pkgname='saltspec')
