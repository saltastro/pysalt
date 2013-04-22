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
SPECSKYNORMALIZE is a program to apply a correction to the data according
to the response of sky lines in the verticle direction.   The response
is either calculated via the sky lines in the image or they can also
be supplied via a file

Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)    1.0       12 Feb 2012

TODO
----


LIMITATIONS
-----------

"""
# Ensure python 2.5 compatibility
from __future__ import with_statement

import os, sys, math, time
import numpy as np

from pyraf import iraf 
import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging
from salterror import SaltError, SaltIOError
from scipy import ndimage as nd

from PySpectrograph.WavelengthSolution import WavelengthSolution
from PySpectrograph.Spectra import apext

import spectools as st
import mostools as mt

from spectools import SALTSpecError
from AutoIdentify import getwsfromIS



debug=True



# -----------------------------------------------------------
# core routine

def specskynormalize(images, outfile, outpref, response=None,
                  startext=0, clobber=False,logfile='salt.log',verbose=True):

   with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       #read in the response function
       if response is not None:
           norm_arr=readresponse(response)

       # Identify the lines in each file
       for img, ofile in zip(infiles, outfiles):

           #open the image
           hdu=saltio.openfits(img)

           for i in range(startext, len(hdu)):
               if hdu[i].name=='SCI':
                   log.message('Normalizing extension %i in  %s' % (i, img))
                   #things that will change for each slit
 
                   #set up the data for the source
                   try:
                       data=hdu[i].data
                   except Exception, e:
                       message = 'Unable to read in data array in %s because %s' % (img, e)
                       raise SALTSpecError(message)

                   #set up the center row
                   if rstart=='middlerow':
                       ystart=int(0.5*len(data))
                   else:
                       ystart=rstart

                   #set up the xarr array based on the image
                   xarr=np.arange(len(data[ystart]), dtype='int64')

                   #calculate the transformation
                   ImageSolution=arcstraight(data, xarr, ystart, ws=None, function=function, order=order,dcoef=dcoef,
                                             rstep=rstep, nrows=nrows, ndstep=ndstep, log=log, verbose=verbose)


                   if outfile and len(ImageSolution):
                       writeIS(ImageSolution, outfile, dateobs=dateobs, utctime=utctime, instrume=instrume, 
                               grating=grating, grang=grang, grasteps=grasteps, arsteps=arsteps, \
                               arang=arang, rfilter=rssfilter, slit=slit, xbin=xbin,      \
                               ybin=ybin, objid=objid, \
                               filename=img, log=log, verbose=verbose)


#------------------------------------------------------------------
# Find the solution for lines in a file

def arcstraight(data, xarr, istart, ws=None, function='poly', order=3,  
             rstep=1, nrows=1, dcoef=None, ndstep=50, log=None, verbose=True):
   """For a given image, assume that the line given by istart is the fiducial and then calculate
      the transformation between each line and that line in order to straighten the arc

      returns Wavlenght solution
   """
   ImageSolution={}

   #extract the central row
   oxarr=xarr.copy()
   ofarr=data[istart]
   print function, order
   ws=WavelengthSolution.WavelengthSolution(xarr, xarr, function, order)
   ws.fit()
   print ws.coef
   ImageSolution[istart]=ws
   if dcoef is None: 
       docef=ws.coef*0.0
       dcoef[0]=10.0
   else:
       dcoef=np.array(dcoef)
   print dcoef

   data=nd.gaussian_filter(data, 3)

   #now step around the central row
   for i in range(rstep,int(0.5*len(data)), rstep):
       for k in [istart-i, istart+i]:
           lws=getwsfromIS(k, ImageSolution)
           xarr=np.arange(len(data[k]))
           farr=apext.makeflat(data, k, k+nrows)
           nws=st.findxcor(xarr, farr, oxarr, ofarr, lws, dcoef=dcoef, ndstep=ndstep, best=True, inttype='interp', debug=False)
           ImageSolution[k]=nws
           print k, nws.coef

   return ImageSolution

def writeIS(ImageSolution, outfile, dateobs=None, utctime=None, instrume=None,       \
                       grating=None, grang=0.0, grasteps=None, objid=None,     \
                       arang=0.0, arsteps=None,rfilter=None, slit=None, xbin=2, ybin=2,\
                       filename=None, log=None, verbose=False):



   #set up the list of solutions to into an array
   key_arr=np.array(ImageSolution.keys())
   arg_arr=key_arr.argsort()

   #set up the wavelength solution
   ws=ImageSolution[key_arr[0]]
   ws_arr=np.zeros((len(arg_arr),len(ws.coef)+1), dtype=float)

   #write the solution to an array 
   for j,i in enumerate(arg_arr):
      if  isinstance(ImageSolution[key_arr[i]], WavelengthSolution.WavelengthSolution):
          function=ImageSolution[key_arr[i]].function
          order=ImageSolution[key_arr[i]].order
          ws_arr[j,0]=key_arr[i]
          ws_arr[j,1:]=ImageSolution[key_arr[i]].coef

   #write header to the file that should include the order and function
   if os.path.isfile(outfile):
      dout=open(outfile, 'a')
   else:
      dout=open(outfile, 'w')

   msg='#WS: Wavelength solution for image %s\n' % filename
   msg+= '#The following parameters were used in determining the solution:\n'
   msg+= '#name=%s\n' % filename
   msg+= '#time-obs=%s %s\n' % (dateobs, utctime )
   msg+= '#instrument=%s\n' % instrume
   msg+= '#grating=%s\n' % grating.strip()
   msg+= '#graang=%s\n' % grang   
   msg+= '#gratilt=%s\n' % grasteps
   msg+= '#arang=%s\n' % arang 
   msg+= '#artilt=%s\n' % arsteps
   msg+= '#filter=%s\n' % rfilter.strip() 
   if objid:  msg+= '#slitid=%s\n' % objid
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
          
   return 



# main code 

parfile = iraf.osfn("saltspec$specskynormalize.par") 
t = iraf.IrafTaskFactory(taskname="specskynormalize",value=parfile,function=specskynormalize, pkgname='saltspec')
