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

#!/usr/bin/env python

"""
SALTCRCLEAN provides cosmic ray cleaning for multi-dimensial fits files.  The
user can select three different methods for cleaning the data including fast, 
median, and edge.  The fast detection algorithm locates peaks in the 
distribution and removes pixels which are substantial above their neighbors.  
The median method is similar to iraf.imred.crutil.crmedian.   The edge
detection method invokes the algorithm outlined in van Dokkum (2001)


Author                 Version      Date
-----------------------------------------------
S M Crawford (SAAO)    0.3          10 Apr 2008

Updates
----------------------------------------------
9 Jun 11  SMC   -Updated to new error handling
                -Updated to be compatible with new libraries
                -Improved crmedian to run faster
                -re-factored so it now each task only returns the crarr
"""

from __future__ import with_statement


import os, time
import numpy as np

try:
   import multiprocessing as mp
except:
   mp=None

from pyraf import iraf

import saltstat, salttran
import saltsafekey as saltkey
import saltsafeio as saltio

from salterror import SaltError
from saltsafelog import logging, history


#from scipy.signal import fftconvolve as conv2d
#from scipy.signal import convolve2d as conv2d
from scipy.ndimage.filters import convolve as conv2d
from scipy.ndimage.filters import generic_filter

debug=True

# -----------------------------------------------------------
# core routine

def saltcrclean(images,outimages,outpref,crtype='fast',thresh=5,mbox=3,         \
                bthresh=3, flux_ratio=0.2, bbox=11, gain=1, rdnoise=5, fthresh=5,\
                bfactor=2, gbox=3, maxiter=5, multithread=False, update=True,
                clobber=True,  logfile='salt.log', verbose=True):


   with logging(logfile,debug) as log:


       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       #verify that the input and output lists are the same length
       saltio.comparelists(infiles,outfiles,'Input','output')

       #check to see if multithreading is available
       if mp:
          pass
       else:
          multithread=False
          log.warning('multiprocessing module is not available.  Setting multiththread=False')


       # Begin processes each file
       for infile,outfile in zip(infiles,outfiles):

           #open the infile
           struct=saltio.openfits(infile)

           #clean the cosmic rays
           if multithread and len(struct)>1:
               struct=multicrclean(struct, crtype, thresh, mbox, bbox, bthresh, flux_ratio, \
                          gain, rdnoise, bfactor, fthresh, gbox, maxiter, log, verbose=verbose)
           else:
               struct=crclean(struct, crtype, thresh, mbox, bbox, bthresh, flux_ratio, \
                          gain, rdnoise, bfactor, fthresh, gbox, maxiter, update, log, verbose=verbose)
          
           #log the call
           #log.message('Cleaned %i cosmic rays from %s using %s method' % (totcr, infile, crtype), with_header=False)
           log.message('', with_header=False, with_stdout=verbose)

           #add house keeping keywords
           saltkey.put('SAL-TLM',time.asctime(time.localtime()), struct[0])
   
           #add the history keyword
           fname, hist=history(level=1, wrap=False)
           saltkey.history(struct[0],hist)

           #write out the file
           saltio.writefits(struct, outfile, clobber=clobber)

           #close the image
           saltio.closefits(struct)

def multicrclean(struct, crtype='fast', thresh=5, mbox=5, bbox=11, bthresh=3, flux_ratio=0.2, \
              gain=1, rdnoise=5, bfactor=2, fthresh=5, gbox=0, maxiter=1, update=True, log=None, verbose=True):
   """MULTICRCLEAN cleans SALT-like data of cosmic rays.  The user has 
      three different choices for the type of cosmic ray cleaning being
      fast, median, and edge.   The process is set up to use multithreading for 
      quick processing of the data.  

           crytpe--type of cosmic ray cleaning. Either fast, median, or 
           thresh--threshold for detecting cosmic rays
           mbox--box for median cleaning
           bbox--background box for median measurement
           bthresh--threshold for iterating on background calculation
           flux_ratio--ratio of fluxes for 'fast' method
           gain--gain of images--set to None to read in from header
           rdnoise--read noise of images--set to None to read in from header
           bfactor--block replication factor for 'edge' method
           fthresh--threshold for excluding compact sources (edge only)
           gbox--Window size to grow sources.  gbox=0 for no growth of cosmic rays
           maxiter--maximum number of iterations

           return struct
   """
   #setup the image name
   infile=saltkey.getimagename(struct[0])

           
   #count the CR
   totcr=0

   #print out the header for the log
   if log:
       message = '%28s %11s' % ('HDU','COSMICRAYS')
       log.message('\n      ---------------------------------------------------', \
            with_header=False, with_stdout=verbose)
       log.message(message, with_header=False, with_stdout=verbose)
       log.message('      ---------------------------------------------------', \
            with_header=False, with_stdout=verbose)

   
   task_list=[]
   nproc=0
   for hdu in struct:
       if hdu.name=='SCI':
           task_list.append((hdu.data, crtype, thresh, mbox, bbox, bthresh, flux_ratio, gain, rdnoise, bfactor , fthresh, gbox, maxiter))
           nproc+=1

   #set up the multi-thread
   p=mp.Pool()


   results=[p.apply_async(cleancosmicrays, i) for i in task_list]
   p.close()


   for i, hdu in enumerate(struct):
       if hdu.name=='SCI':
           #set up the cosmic ray array
           crarr=results[i-1].get()

           #update the frame for the various values
           mask=(crarr>0)
           if update: struct[i].data[mask]=crarr[mask]

           #track the number of cosmic rays
           ncr=mask.sum()
           totcr += ncr
           
           #if verbose print out information
           if log:
               message='%25s[%1d]  %i' % (infile, i, ncr)
               log.message(message, with_header=False, with_stdout=verbose)

           #correct the BPM frame
           if saltkey.found('BPMEXT', struct[i]):
                b_i=saltkey.get('BPMEXT', struct[i])
                try:
                   struct[b_i].data[mask]=1
                except Exception, e:
                   msg='Cannot update the BPM frame in %s[%i] because %s' % (infile, b_i, e)
                   raise SaltError(msg)
   return struct

def crclean(struct, crtype='fast', thresh=5, mbox=5, bbox=11, bthresh=3, flux_ratio=0.2, \
              gain=1, rdnoise=5, bfactor=2, fthresh=5, gbox=0, maxiter=1, update=True, log=None, verbose=True):
   """CRCLEAN cleans SALT-like data of cosmic rays.  The user has 
           three different choices for the type of cosmic ray cleaning being
           fast, median, and edge.  

           crytpe--type of cosmic ray cleaning. Either fast, median, or 
           thresh--threshold for detecting cosmic rays
           mbox--box for median cleaning
           bbox--background box for median measurement
           bthresh--threshold for iterating on background calculation
           flux_ratio--ratio of fluxes for 'fast' method
           gain--gain of images--set to None to read in from header
           rdnoise--read noise of images--set to None to read in from header
           bfactor--block replication factor for 'edge' method
           fthresh--threshold for excluding compact sources (edge only)
           gbox--Window size to grow sources.  gbox=0 for no growth of cosmic rays
           maxiter--maximum number of iterations

           return struct
   """
   #setup the image name
   infile=saltkey.getimagename(struct[0])

           
   #count the CR
   totcr=0

   #print out the header for the log
   if log:
       message = '%28s %11s' % ('HDU','COSMICRAYS')
       log.message('\n      ---------------------------------------------------', \
            with_header=False, with_stdout=verbose)
       log.message(message, with_header=False, with_stdout=verbose)
       log.message('      ---------------------------------------------------', \
            with_header=False, with_stdout=verbose)

   #cosmic ray clean each extension
   for i in range(len(struct)):
          
       #only clean the cosmic rays if it is a SCI extension or a single extension
       if struct[i].name=='SCI' or len(struct)==1:           
           #for the edge method, get the gain and rdnoise from the fits header
           #if they are not set
           if crtype=='edge':
               if gain==None: gain=saltkey.get('GAIN',struct[i])
               if rdnoise==None: rdnoise=saltkey.get('RDNOISE',struct[i])

           #get all the cosmic rays from an array
  
           crarr=cleancosmicrays(struct[i].data, crtype, thresh, mbox, bbox, bthresh, flux_ratio, \
                                gain, rdnoise, bfactor, fthresh, gbox, maxiter)
           #update the frame for the various values
           mask=(crarr>0)
           if update: struct[i].data[mask]=crarr[mask]

           #track the number of cosmic rays
           ncr=mask.sum()
           totcr += ncr
           
           #if verbose print out information
           if log:
               message='%25s[%1d]  %i' % (infile, i, ncr)
               log.message(message, with_header=False, with_stdout=verbose)

           #correct the BPM frame
           if saltkey.found('BPMEXT', struct[i]):
                b_i=saltkey.get('BPMEXT', struct[i])
                try:
                   struct[b_i].data[mask]=1
                except Exception, e:
                   msg='Cannot update the BPM frame in %s[%i] because %s' % (infile, b_i, e)
                   raise SaltError(msg)

   return struct

#-----------------------------------------------------------
# Clean cosmic rays

def cleancosmicrays(arr, crtype='fast', thresh=5, mbox=3, bbox=5, bthresh=5, flux_ratio=0.2, gain=1, \
                    rdnoise=5, b_factor=2, fthresh=3, gbin=0, max_iter=1):
   """Clean the cosmic rays from an input array using three different types of methods
   """
   niter=0
   npix=1
   onpix=0

   #makesure the data array is in floating point
   try:
       arr=arr*1.0
   except:
       message='Array cannot be convert to the correct data type'
       raise SaltError(message)

   #measure the mean values for the image--leaving this outside 
   #because it doesn't need to be done on each iter loop
   if crtype=='fast':
       mean, midpt, bsigma=saltstat.iterstat(arr, bthresh, max_iter)

   #set up the array
   crarr=arr*0.0

   while niter < max_iter and npix > onpix:
       #apply the correct function to idenitfy cosmic rays
       if crtype=='median':
           crarr +=crmedian(arr, thresh, mbox, bbox)
       elif crtype=='fast':
           crarr +=crfast(arr, mean, bsigma, thresh, mbox, flux_ratio)
       elif crtype=='edge':
           crarr +=credge(arr, thresh, gain, rdnoise, b_factor, fthresh)
       else:
           message='CRTYPE %s is not a valid value' % crtype
           raise SaltError(message)

       #update the array
       mask=(crarr>0)
       arr[mask]=crarr[mask]

       #count the number of cosmic rays identified
       onpix=npix
       npix=crarr.sum()
       niter +=1
       
   #grow the result
   if gbin>0: crarr=crgrow(crarr, grad=gbin) 
        

   return crarr

#-----------------------------------------------------------
# set up a boxes dimensions

def setbox(x, y, mbox, xlen, ylen):
   """set up the edges of the array"""
   y1=max(0, y-mbox)
   y2=min(y+mbox+1, ylen-1)
   x1=max(0, x-mbox)
   x2=min(x+mbox+1, xlen-1)

   return x1, x2, y1, y2


#-----------------------------------------------------------
# Fast identification of a cosmic ray

def idcray(arr, flux_ratio, y, x, ylen, xlen, mbox):

    #set up the box dimensions
    x1, x2, y1, y2 = setbox(x, y, mbox, xlen, ylen)

    #is the source the brightest in the box
    if arr[y,x] < arr[y1:y2,x1:x2].max(): return y, x, None

    #If the pixel to median flux  is less than
    # the flux ratio, reject the objects as a comsic ray
    median_pix=saltstat.median2d(arr[y1:y2,x1:x2])
    if abs(median_pix)/arr[y,x] > flux_ratio: return y, x, None

    return y,x,median_pix

#-----------------------------------------------------------
# Fast clean of the cosmic rays

def crfast(arr, mean, bsigma, thresh, mbox, flux_ratio):
    """Provide fast cosmic ray rejection based on identifying isolated peaks 


       crarr--array with zero for non-cosmic rays and positive value 
              equal to the median of the surrounding pixels if a cosmic ray

       returns crarr
    """

    #set up the size of the box.  Minimimum size is a 3x3 box around
    # the central pixel
    mbox=max(3,mbox)
    mbox=int(mbox/2)

    ylen=len(arr)
    xlen=len(arr[0])

    #set up the cr identification array
    crarr=0.0*arr
 
    #Detect all possible objects on the image
    rarr=arr-mean
    darr= np.where(rarr>bsigma*thresh)
          
    #Out of those detections, eliminate those that aren't cosmic rays
    for i in range(len(darr[0])):
        y=darr[0][i]
        x=darr[1][i]
        y, x, mpix=idcray(rarr, flux_ratio, y, x, ylen, xlen, mbox)
        if mpix is not None: crarr[y,x]=mpix+mean

    return crarr

def sigma_func (arr):
    return arr.std()
    arr=arr.ravel()
    arr=np.sort(arr)
    larr=len(arr)
    x1=int(0.159*larr)
    x2=int(0.841*larr)
    std=(arr[x2]-arr[x1])/2.0
    return std

def rolling_window(a, window):
   from numpy.lib.stride_tricks import as_strided
   shape = (arr.shape[0] - window + 1, window, arr.shape[1])
   strides = (strides[0],) + strides
   return as_strided(arr, shape=shape, strides=strides)



#-----------------------------------------------------------
# Median clean the cosmic rays

def crmedian(arr, thresh, mbox, bbox):
    """Identify cosmic rays through median technique.  Similar implimentation to
       crmedian in iraf.imred.crutil.crmedian
    """
    #Check to make sure the background box is an appropriate size
    #If it is too small, then insufficient statistics are generated
    if bbox <  10: bbox=10
    cbox=int(0.5*bbox)

    #make the background image
    #barr=generic_filter(arr,sigma_func, size=(bbox,bbox))
    barr=arr*0.0+arr.std()
    xlen=len(arr[0])
    ylen=len(arr)
    for i in range(mbox,xlen,bbox):
     for j in range(mbox,ylen,bbox):
         x1,x2,y1,y2=setbox(i,j,cbox, xlen, ylen)
         barr[y1:y2,x1:x2]=sigma_func(arr[y1:y2,x1:x2])

    #Median smooth the image
    marr=saltstat.median_image(arr, mbox)

    #Find the residual image
    rarr=(arr-marr)/barr

    #identify all sources
    crarr=marr*(rarr > thresh)

    return crarr

#-----------------------------------------------------------
# Edge clean the cosmic rays

def credge(arr, thresh, gain, rdnoise, b_factor, fthresh):
    """Identify cosmic rays following the van Dokkem (2001) prescription
    """

    min_limit=0.01

    #set up the convolution kernel
    f_conv=np.array([[0,-1,0],[-1,4,-1],[0,-1,0]])

    #set up the array
    shape=arr.shape
    sdata=arr

    #rebin the data
    newshape=(b_factor*shape[0],b_factor*shape[1])
    ldata=salttran.rebin_factor(sdata,newshape)

    #convolve with f_conv
    #ldata=conv2d(ldata,f_conv,mode='same')
    ldata=conv2d(ldata,f_conv)
    mask = (ldata >= 0)
    ldata = ldata * mask

    #return to the original binning
    ldata = salttran.blockave(ldata,shape)

    #create noise model
    meddata =saltstat.median_image(sdata,5)
    noise = abs(meddata)
    noise = (gain*noise+rdnoise**2)**0.5/gain

    #create S/N image
    odata = ldata / noise / b_factor

    #remove extended objects
    odata = odata - saltstat.median_image(odata,5)

    #select objects
    masks = (odata>thresh)

    #remove compact bright sources
    mdata = saltstat.median_image(sdata,5)
    fdata = mdata - saltstat.median_image(mdata,7)
    fdata = fdata / noise

    # set a minimum value for all pixels so no divide by zero problems
    indf = np.where(fdata < min_limit)
    fdata[indf]=min_limit

    fdata = odata * masks/fdata
    maskf = (fdata > fthresh)

    #make the list of cosmic rays
    mask =  masks*maskf

    #identify the cosmic rays
    crarr=mdata*mask

    return crarr



def crgrow(crarr, grad=3):
   """If a cosmic ray has been identified in the source, delete the pixels around
      it according to gradius.  
   """

   darr=np.where(crarr>0)
   xlen=len(crarr[0])
   ylen=len(crarr)

   #set it up as a box
   if grad < 3: grad=3
   mbox=int(grad/2)

   for i in range(len(darr[0])):
        y=darr[0][i]
        x=darr[1][i]
        x1, x2, y1, y2 = setbox(x, y, mbox, xlen, ylen)
        crarr[y1:y2,x1:x2]=crarr[y,x]

   return crarr

#-----------------------------------------------------------
# Function to calculate the sigma at each pixel


# -----------------------------------------------------------
# main code
if not iraf.deftask('saltcrclean'):
   parfile = iraf.osfn("saltred$saltcrclean.par")
   t = iraf.IrafTaskFactory(taskname="saltcrclean",value=parfile,function=saltcrclean, pkgname='saltred')
