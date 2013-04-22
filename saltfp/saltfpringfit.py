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
SALTFPRINGFIT is a program to fit rings from a Fabry Perot image.
The first step is to identify the approximate position, center,
radius, and FWHM of a ring.  Once the first approximate
for the ring is determined it will fit the ring using all
of the points available to determine the ring parameters.

# Author                 Version      Date
# -----------------------------------------------
# S. M. Crawford (SAAO)    0.1       21 Jul 2011

"""

from __future__ import with_statement

import os
import pyfits
import math
import numpy as np
import scipy.ndimage as nd

from pyraf import iraf 
from pyraf.iraf import pysalt
import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging
from saltstat import iterstat


from display import display
from PySpectrograph.Spectra import findobj

from salterror import SaltError

from FPRing import FPRing

import pylab as pl

debug=True

# -----------------------------------------------------------
# core routine


def saltfpringfit(images, outfile, section=None, bthresh=5, niter=5,
                displayimage=True, clobber=True,logfile='salt.log',verbose=True):

   with logging(logfile,debug) as log:


       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # read in the section
       if section is None:
           section=saltio.getSection(section)
           msg='This mode is not supported yet'
           raise SaltError(msg)
       else:
           section=saltio.getSection(section)
       print section

       # open each raw image file
       for img in infiles:

          #open the fits file
	  struct=saltio.openfits(img)
          data=struct[0].data

          #only keep the bright pixels
          y1,y2,x1,x2=section
          bmean, bmedian, bstd=iterstat(data[y1:y2,x1:x2], sig=bthresh, niter=niter, verbose=False)
          message="Image Background Statistics\n%30s %6s %8s %8s\n%30s %5.4f %5.4f %5.4f\n" %  \
                ('Image', 'Mean', 'Median', 'Std',img, bmean, bmedian, bstd)
          log.message(message, with_stdout=verbose)
    
          mdata=data*(data-bmean>bthresh*bstd)

          #prepare the first guess for the image
          ring_list=findrings(data, thresh=5, niter=5, minsize=10)

          if displayimage:
             regfile=img.replace('.fits', '.reg')
             print regfile
             if clobber and os.path.isfile(regfile): fout=saltio.delete(regfile)
             fout=open(regfile, 'w')
             fout.write("""# Region file format: DS9 version 4.1
# Filename:  %s
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
physical
""" % img)
             for ring in ring_list:
                 fout.write('circle(%f, %f, %f)\n' % (ring.xc(),ring.yc(),ring.prad()))
                 fout.write('circle(%f, %f, %f)\n' % (ring.xc(),ring.yc(),ring.prad()-5*ring.sigma()))
                 fout.write('circle(%f, %f, %f)\n' % (ring.xc(),ring.yc(),ring.prad()+5*ring.sigma()))

             fout.close()
             display(img, catname=regfile, rformat='reg')

          #write out the result for viewing
          struct[0].data=mdata 
          saltio.writefits(struct, 'out.fits', clobber=True)

	  message = 'Ring Parameters'
          log.message(message)

def findpeaks(data, fpeak=0.8,minsize=10):
    """Find peakes median filters an image and finds the peaks in the image
    """
    #median filter the image
    mdata=nd.filters.median_filter(data, size=minsize)

    #take the 80% points and find the peaks
    mask=(mdata>fpeak*mdata.max())

    #find all the objects
    obj_arr, obj_num=nd.label(mask)
    
    print obj_arr, obj_num
    #ypeaks
    peaks=[]
    for i in range(obj_num):
        pid=np.where(obj_arr==i+1)[0]
        peaks.append((pid.min(), pid.max()))

    return peaks
    

def findrings(data, thresh=5, niter=5, minsize=10):
    """findrings makes a rough calculation for the parameters of the rings
       based on single line cuts through the data.  It returns a list of rings
    """
    ring_list=[]

    #first guess the middle is in the middle of the data
    xc=int(0.5*len(data[0]))
    yc=int(0.5*len(data))
    print xc,yc
    #take a look at the y cut through the data
    xdata=data[yc,:]
    #take a look through the xdata.  check for the same thing and make sure they are consistent
    ydata=data[:,xc]
    #xdata=nd.filters.median_filter(xdata, size=minsize)
    #ydata=nd.filters.median_filter(ydata, size=minsize)
    #get rid of all the lower points
    #find the peaks in the data
    print xdata.mean(), xdata.max(), xdata.std()
    #ypeak_list=findobj.findLines(ydata, method='median', thresh=thresh, niter=niter, minsize=minsize)
    #xpeak_list=findobj.findLines(xdata, method='median', thresh=thresh, niter=niter, minsize=minsize)
    ypeak_list=findpeaks(ydata, 0.4, 10)
    xpeak_list=findpeaks(xdata, 0.4, 10)
    print xpeak_list 
    print ypeak_list
    pl.figure()
    pl.plot(ydata)
    pl.savefig('out.png')

    if abs(len(ypeak_list)-len(xpeak_list))>1: 
       msg="Non-symmetrically rings in the image"
       #raise SaltError(msg)

    nrings=max(len(ypeak_list)/2, len(xpeak_list)/2)

    #if one peak: no rings.  If two peaks: one ring, if four peaks: four rings
    if nrings<1:
       msg="No rings detected in image"
       raise SaltError(msg)
    elif nrings==1:
       msg="One ring dected in image"
       print msg
    else:
       msg="%i rings found in image" % nrings
       print msg


    #loop through the image and determine parameters of rings
    for i in range(0,nrings,2):
       #determine the y-center
       try:
           y1,y2=ypeak_list[i]
           yarr=np.arange(x1,x2)
           ypa=y1+ydata[y1:y2].argmax()
           ysiga=(abs(np.sum((yarr-ypa)**2*ydata[y1:y2])/ydata[y1:y2].sum()))**0.5
           y1,y2=ypeak_list[i+1]
           ypb=y1+ydata[y1:y2].argmax() 
           ysigb=(abs(np.sum((yarr-ypb)**2*ydata[y1:y2])/ydata[y1:y2].sum()))**0.5
           yc=0.5*(ypa+ypb)
           ymax=max(ydata[ypa], ydata[ypb])
           yrad=0.5*abs(ypb-ypa)
           ysig=0.5*(ysiga+ysigb)
       except:
           yc=yc
           yrad=0
           ysig=0
           ymax=ydata.max()

       #determine the x-center
       try:
           x1,x2=xpeak_list[i]
           xarr=np.arange(x1,x2)
           xpa=x1+xdata[x1:x2].argmax()
           xsiga=(abs(np.sum((xarr-xpa)**2*xdata[x1:x2])/xdata[x1:x2].sum()))**0.5
           x1,x2=xpeak_list[i+1]
           xpb=x1+xdata[x1:x2].argmax() 
           xarr=np.arange(x1,x2)
           xsigb=(abs(np.sum((xarr-xpb)**2*xdata[x1:x2])/xdata[x1:x2].sum()))**0.5
           xc=0.5*(xpa+xpb)
           xmax=max(xdata[xpa], xdata[xpb])
           xsig=0.5*(xsiga+xsigb)
           xrad=0.5*abs(xpa-xpb)
       except:
           xc=yc
           xrad=0
           xsig=0
           xmax=xdata.max()

       #print the results
       print xc,yc,max(yrad,xrad), max(xmax,ymax), max(xsig,ysig)
        
       ring_list.append(FPRing(xc, yc, max(yrad,xrad), max(xmax,ymax), max(xsig,ysig)))
 


    #if not take the lower of the two
 
    return ring_list

# -----------------------------------------------------------
# main code 

parfile = iraf.osfn("saltfp$saltfpringfit.par") 
t = iraf.IrafTaskFactory(taskname="saltfpringfit",value=parfile,function=saltfpringfit, pkgname='saltfp')
