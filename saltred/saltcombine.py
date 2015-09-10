#!/usr/bin/env python

"""

 saltcombine is a tool for combining images into
 a master image.  The program has a number of options for how the
 images are combine including the combination method and any
 rejection algorithm.  

 The task assume that all input images are flats and these are
 the images that you want combined.  Options for combining the 
 sources are average or median.  Options for rejection include
 none, minmax, ccdclip, sigclip, and avsigclip.  The program
 should be able to handle image masks, zeropoints, weights, 
 and scales

 In order to process the rejection method or the median method
 it will require all of the data to be read in.  But a fast 
 non-memory intensive method to combine the data can be accomplished
 using only the average method  with only limited rejection.

Author                 Version      Date
-----------------------------------------------
Steve Crawford (SAAO)    1.0       21 Aug 2009

TODO
-----------------------------------------------------------------

Updates
-----------------------------------------------------------------
 
"""

from __future__ import with_statement

import os
import numpy as np
from astropy.io import fits
from pyraf import iraf 

import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging, history

from salterror import SaltError

import gc
debug=True

# -----------------------------------------------------------
# core routine

def saltcombine(images,outimage, method='average', reject=None, mask=True,     \
                weight=True, blank=0, scale=None, statsec=None, lthresh=3,    \
                hthresh=3, clobber=False, logfile='salt.log',verbose=True):

   with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       #set reject as None
       reject=saltio.checkfornone(reject)
       if reject is not None: reject=reject.lower()

       #set scale 
       scale=saltio.checkfornone(scale)
       if scale is not None: scale=scale.lower()

       #set statsec
       statsec=saltio.checkfornone(statsec)
       statsec=saltio.getSection(statsec, iraf_format=True)

       #Fast combine the images
       outstruct=imcombine(infiles, method=method, reject=reject, mask=mask, \
                           weight=weight, blank=blank, scale=scale,          \
                           statsec=statsec, lthresh=lthresh, hthresh=hthresh)

       # housekeeping keywords
       fname, hist=history(level=1, wrap=False)
       saltkey.housekeeping(outstruct[0],'SCOMBINE', 'File Combined by SALTCOMBINE', hist)

       # write FITS file               
       saltio.writefits(outstruct, outimage)

def imcombine(infiles, method='average', reject=None, mask=True,  weight=True, \
              blank=0, scale=None, statsec=None, lthresh=3,  hthresh=3):
   """Fast combine is for special circumstances where the files and the combination 
       can be down one by one and not all of the data need to be read in all at once.
       It computes the average of the files input in infile

       Input Variables:
       infiles:  List of files to combine

       method:  Combination method to use, either average or median

       reject:  Reject outliers using different methods
                -ccdclip--reject the pixels according to the CCD parameters
                -sigclip--reject the pixels according to a sigma clipping
                
       mask:  Whether to correct the data using the BPM frame

       weight:  If true, weight by the inverse of the variant frame

       blank:  Value if zero pixels are combined

       scale:  Method to scale the images by

       statsec:  Region for determining statistics if images are to be scaled
 
       lthresh: low rejection threshold
      
       hthresh: high rejection threshold

       Return Variables:
       
       hdustruct:  Output struct 
   """
   hdu_list=[]

   #read in all of the data
   for infile in infiles:
       #hdu_list.append(saltio.openfits(infile))
       hdu_list.append(fits.open(infile, memmap=True))

   #Copy the first image and create the HDU for the other images
   try:
       outhdu=hdu_list[0]
   except Exception, e:
       message='Cannot create output hduList because %s' % e
       raise SaltError(message)
   #determine the size of a single extension.  Assumes all arrays are the same size
   hdu=hdu_list[0]

   #For each science extension, construct the array of arrays to go into that
   #Three arrays need to be constructed including one for the data array,
   #one for the inv variance array, and one for the BPM


   #This assumes all the data are the same shape and size
   nimages=len(hdu_list)

   mem_lim=0.2
   if reject or weight or mask: 
      mem_lim=mem_lim/4.0
  
   for i in range(len(outhdu)):
       if outhdu[i].name=='SCI':
           nb=outhdu[i].data.nbytes
           outhdu[i].data *= 0.0
           #let's set the limit at 0.2 GB of memory
           if nb*nimages/1e9>mem_lim: 
              imod=int((nb*nimages/1e9)/mem_lim)+1
              y2,x2=outhdu[i].data.shape
              ystep=y2/imod
              xstep=x2/imod
              for j in range(imod):
                for k in range(imod):
                  datasec=[j*ystep, min(y2,(j+1)*ystep), k*xstep, min(x2, (k+1)*xstep)]
                  outhdu=hducombine(hdu_list, outhdu, ext=i, method=method, datasec=datasec,
                         reject=reject, mask=mask,weight=weight, scale=scale, statsec=statsec, blank=blank,
                         lthresh=lthresh, hthresh=hthresh) 
              #this step is needed to remove data arrays that are read into memory
              #--not sure what object they are linked to, but this works
              try:
                 for j in range(nimages): del hdu_list[j][i].data
              except:
                 pass
           else:  
              outhdu=hducombine(hdu_list, outhdu, ext=i, method=method, datasec=None, 
                         reject=reject, mask=mask,weight=weight, scale=scale, statsec=statsec, 
                         lthresh=lthresh, hthresh=hthresh)
           

   #Add any header frames
   saltkey.new('NCOMBINE',nimages,'Number of images combines', outhdu[0])

   return outhdu


def hducombine(hdu_list, outhdu, ext, method='average', datasec=None, reject=None, mask=False, weight=False, scale=None, statsec=None, blank=0, lthresh=3, hthresh=3):
   """Combine a set of images in imlist 

   """
   #set up i as the extionsion variable as  shortcut
   i=ext
   nimages=len(hdu_list)
   #set all the data arrays
   data_list=[]
   gain=np.zeros(nimages)
   rdnoise=np.zeros(nimages)
   varext=None
   bpmext=None


   #set the datasec in case it is none as the full image
   if datasec is None:
      sh=outhdu[ext].data.shape
      y1,x1=(0,0)
      y2,x2=outhdu[i].data.shape
   else:
      y1,y2,x1,x2=datasec
        
   dshape=outhdu[i].data[y1:y2,x1:x2].shape
   dtype=outhdu[i].data[y1:y2,x1:x2].dtype
   wei_arr=outhdu[i].data[y1:y2,x1:x2]*0.0

   #check for variance frame
   if saltkey.found('VAREXT', outhdu[i]) and weight:
       varext=saltkey.get('VAREXT', outhdu[i])
       var_list=[]
   #check for bad pixel mask
   if saltkey.found('BPMEXT', outhdu[i]) and mask:
       bpmext=saltkey.get('BPMEXT',outhdu[i])
       bpm_list=[]

   #create the lists of arrays and scale the arrays if requests
   mean_scale=0
   data_arr=np.zeros((nimages, dshape[0], dshape[1]), dtype=dtype)
   for j in range(nimages):
       data_arr[j,:,:]=hdu_list[j][i].data[y1:y2,x1:x2]
       #calculate the scale 
       if scale:
          scale_val=CalculateScale(data_arr, scale, statsec)
          mean_scale += scale_val
       else:
          scale_val=1
          mean_scale+=1

       if varext:
           var_list.append(hdu_list[j][varext].data[y1:y2,x1:x2]/scale_val)
       if bpmext:
           bpm_list.append(hdu_list[j][bpmext].data[y1:y2,x1:x2])

       #get the gain and rdnoise
       if reject=='ccdclip':
           if saltkey.found('GAINMULT', hdu_list[j][i]):
               gain[j]=1
           else:
               gain[j]=saltkey.get('GAIN', hdu_list[j][i])
           rdnoise[j]=saltkey.get('RDNOISE', hdu_list[j][i])

   #convert the lists to arrays
   if varext:
       var_arr=np.array(var_list)
       ivar_arr=1.0/var_arr
   else:
       var_arr=None
       ivar_arr=None

   if bpmext:
       bpm_arr=np.array(bpm_list)
   else:
       bpm_arr=None
   
   #reject outliers if set
   bpm_arr=RejectArray(data_arr, reject=reject, var=var_arr, bpm=bpm_arr, \
                       lthresh=lthresh, hthresh=hthresh, gain=gain, rdnoise=rdnoise)

   #calculate the average values
   outdata, outwei=CombineArray(data_arr, method=method, ivar=ivar_arr, bpm=bpm_arr)
   outhdu[i].data[y1:y2,x1:x2]=outdata
   if scale is not None:
       mean_scale = mean_scale/nimages
       outhdu[i].data[y1:y2,x1:x2] *= mean_scale

   #create the combine variance frame
   if varext:
       outhdu[varext].data[y1:y2,x1:x2], tmp_arr=CombineArray(var_arr, method=method, ivar=ivar_arr, bpm=bpm_arr)
       del tmp_arr
       if scale is not None:
           outhdu[varext].data *= mean_scale

   #check to see if any of the pixels have no values and replace with blank
   #if wei_arr.any()==0:
   #    wmask=(wei_arr==0)
   #    outhdu[i].data[wmask]=blank


   #create the combine BPM frames
   if bpmext:
       outhdu[bpmext].data=1.0*(wei_arr==0)

   return outhdu
  
def CalculateScale(arr, scale=None, statsec=None):
  """Caculate the scale according to different methods

  Input Variables:
  
  arr:  Input data array

  scale:  Method for calculating the scacle
          -average: calculate the average
          -median: calculate the median

  statsec:  Section over which to calculate the scale

  """
 
  if statsec is None:
     y1,y2,x1,x2=[0,len(arr),0,len(arr[0])]
     data=arr[y1:y2,x1:x2]
  else:
     data=arr

  if scale=='average':
       return np.mean(data)
  if scale=='median':
       return np.median(data)
  else:
       msg='%s is not a type of scale' % scale
       raise SaltError(msg)



def RejectArray(arr, reject=None, var=None, bpm=None, lthresh=3, hthresh=3, gain=1, rdnoise=5):
   """reject outliers in an array.  
      Input Variables:
         arr: Input data array

         reject:  Method for reject bad pixels
            ccdclip--use information about the gain and variance of the array to reject outliers
            sigclip--calculates the average and then rejects either high or low values
        
      Returns:
         bpm:  Data array with bad pixesl set to zero
   """
   #return bpm if reject is None
   if reject is None: return bpm

   #reject the pixels
   if reject=='ccdclip':
      #create mean array
      mean_arr=arr.mean(axis=0)
      #create variance array
      if var is None:
         var=arr.copy()
         for i in range(len(var)):
             var[i]=var[i]*gain[i]+rdnoise[i]**2
      mask=np.ma.mask_or((arr-mean_arr<-lthresh*var),(arr-mean_arr>hthresh*var))
      return 1*np.ma.mask_or((bpm==1),mask)
   elif reject=='sigclip':
      mean_arr=arr.mean(axis=0)
      std_arr=arr.std(axis=0)
      mask=np.ma.mask_or((arr-mean_arr<-lthresh*std_arr),(arr-mean_arr>hthresh*std_arr))
      return 1*np.ma.mask_or((bpm==1),mask)
   else:
       msg='%s is currently not a supported rejection method' % reject
       raise SaltError(msg)

def CombineArray(arr, method='average', ivar=None, bpm=None):
   """Combine--Combine an array of arrays.  The 

    method: type of combination--either average or median
   
    ivar: inverse variance array

    bpm:  bad pixel mask array

   """
   if bpm is None:
       bpm=arr*0.0+1.0
       wei=None
   else:
       # correct the weights for the bad pixel mask
       if ivar is None: ivar=arr*0.0+1.0
       wei=ivar*(1-bpm)
       #TODO: need to add a check to make sure that there are is a place
       #to make sure that one of the weights is at least zero
       check_wei=wei.sum(axis=0)
       wei[0][check_wei==0]=wei.min()
  
   if method=='average':
       c_arr, s_arr=np.average(arr, axis=0, weights=wei, returned=True)
       return c_arr, s_arr
   elif method=='median':
       return np.median(arr, axis=0), bpm.sum(axis=0)
   else: 
       msg='%s is not a method for combining arrays' % method
       raise SaltError(msg)

# -----------------------------------------------------------
# main code 
if not iraf.deftask('saltcombine'):
   parfile = iraf.osfn("saltred$saltcombine.par") 
   t = iraf.IrafTaskFactory(taskname="saltcombine",value=parfile,function=saltcombine, pkgname='saltred')
