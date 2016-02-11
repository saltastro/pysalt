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
 SALTBIAS corrects the bias in SALT CCD data

 Author                 Version      Date
 -----------------------------------------------
 Martin Still (SAAO)    0.1          05 Sep 2006
 S. M Crawford          0.2          05 April 2011



 Updates:
 ------------------------------------------------------------
 12 Dec 2010   Updated the fitting code to use most recent version of interfit
 5 Apr 2011    Updated to use saltsafekey and to use TIME-OBS instead of
               using UTC-OBS, added a catch to using interfit
 10 Sep 2011   Removed master bias from here and moved to saltclean
               Updated the error handling to current form
 

"""

from __future__ import with_statement

import os, string, sys, glob, time
from pyraf import iraf
from pyraf.iraf import pysalt
import numpy as np
import saltstat, saltfit

import saltsafeio as saltio
import saltsafekey as saltkey
from salterror import SaltError, SaltIOError
from saltsafelog import logging, history

debug=True

# Make sure the plotting functions work with an older version of matplotlib
try:
    #import matplotlib
    #matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    #import matplotlib
    #matplotlib.use('Agg')
    import matplotlib.pylab as plt

from matplotlib import font_manager

# -----------------------------------------------------------
# core routine

def saltbias(images,outimages,outpref,subover=True,trim=True,subbias=False,
             masterbias='bias.fits', median=False, function='polynomial', 
             order=3, rej_lo=3, rej_hi=3, niter=10, plotover=False, 
             turbo=False, clobber=False, logfile='salt.log', verbose=True):

   status = 0
   ifil = 0
   ii = 0
   mbiasdata = []
   bstruct = ''
   biasgn = ''
   biassp = ''
   biasbn = ''
   biasin = ''
   filetime = {}
   biastime = {}
   for i in range(1,7):
        filetime[i] = []
        biastime[i] = []

   with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       # are input and output lists the same length?
       saltio.comparelists(infiles,outfiles,'Input','output')

       # Does master bias frame exist?
       # gain, speed, binning and instrument of master bias frame
       if subbias:
           if os.path.isfile(masterbias):
               bstruct = saltio.openfits(masterbias)
           else:
               message = 'Master bias frame %s does not exist' % masterbias
               raise SaltError(message)
       else:
           bstruct=None

       # open each raw image file
       for img, oimg in zip(infiles, outfiles):

           #open the file
           struct = saltio.openfits(img)

           #check to see if it has already been bias subtracted
           instrume,keyprep,keygain,keybias,keyxtalk,keyslot = saltkey.instrumid(struct)

           # has file been biaseded already?
           try:
               key = struct[0].header[keybias]
               message = 'File %s has already been de-biased ' % infile
               raise SaltError(message)
           except:
               pass

           #compare with the master bias to make sure they are the same
           if subbias:
               pass

           #subtract the bias
           struct=bias(struct,subover=subover, trim=trim, subbias=subbias, 
                       bstruct=bstruct, median=median, function=function,
                       order=order, rej_lo=rej_lo, rej_hi=rej_hi, niter=niter,
                       plotover=plotover, log=log, verbose=verbose)

           #write the file out
           # housekeeping keywords
           fname, hist=history(level=1, wrap=False, exclude=['images', 'outimages', 'outpref'])
           saltkey.housekeeping(struct[0],keybias, 'Images have been de-biased', hist)

           # write FITS file
           saltio.writefits(struct,oimg, clobber=clobber)
           saltio.closefits(struct)


def bias(struct,subover=True,trim=True, subbias=False, bstruct=None, 
         median=False, function='polynomial',order=3,rej_lo=3,rej_hi=3,niter=10,
         plotover=False, log=None, verbose=True):
   """Bias subtracts the bias levels from a frame.  It will fit and subtract the overscan
      region, trim the images, and subtract a master bias if required.

      struct--image structure 
      subover--subtract the overscan region
      trim--trim the image
      subbias--subtract master bias
      bstruct--master bias image structure
      median--use the median instead of mean in image statistics
      function--form to fit to the overscan region
      order--order for the function
      rej_lo--sigma  of low points to reject in the fit
      rej_hi--sigma of high points to reject in the fit
      niter--number of iterations
      log--saltio log for recording information
      verbose--whether to print to stdout 
   """
   infile=saltkey.getimagename(struct[0])

   # how many extensions?
   nsciext = saltkey.get('NSCIEXT',struct[0])
   nextend = saltkey.get('NEXTEND',struct[0])
   nccd = saltkey.get('NCCDS',struct[0])

   # how many amplifiers?--this is hard wired
   amplifiers = 2 * nccd


   #log the process
   if subover and log:
        message = '%28s %7s %5s %4s %6s' % \
            ('HDU','Overscan','Order','RMS','Niter')
        log.message('\n     --------------------------------------------------', 
                   with_header=False, with_stdout=verbose)
        log.message(message, with_header=False, with_stdout=verbose)
        log.message('     --------------------------------------------------', 
                   with_header=False, with_stdout=verbose)

   if (plotover): 
       plt.figure(1)
       plt.axes([0.1,0.1,0.8,0.8])
       plt.xlabel('CCD Column')
       plt.ylabel('Pixel Counts (e-)')
       plt.ion()


   #loop through the extensions and subtract the bias
   for i in range(1,nsciext+1):
     if struct[i].name=='SCI':

       #get the bias section
       biassec = saltkey.get('BIASSEC',struct[i])
       y1,y2,x1,x2 = saltio.getSection(biassec, iraf_format=True)
       #get the data section
       datasec = saltkey.get('DATASEC',struct[i])
       dy1,dy2, dx1, dx2 = saltio.getSection(datasec, iraf_format=True)
  
       #setup the overscan region
       if subover:
           yarr=np.arange(y1,y2, dtype=float)
           data=struct[i].data
           odata=struct[i].data[y1:y2,x1:x2]
           if median:
              odata=np.median((struct[i].data[y1:y2,x1:x2]),axis=1)
              olevel=np.median((struct[i].data[y1:y2,x1:x2]))
              saltkey.new('OVERSCAN','%f' % (olevel),'Overscan median value', struct[i])
           else:
              odata=np.mean((struct[i].data[y1:y2,x1:x2]),axis=1)
              olevel=np.mean((struct[i].data[y1:y2,x1:x2]))
              saltkey.new('OVERSCAN','%f' % (olevel),'Overscan mean value', struct[i])

           #fit the overscan region
           ifit=saltfit.interfit(yarr, odata, function=function, \
                                 order=order, thresh=rej_hi, niter=niter)
           try:
               ifit.interfit()
               coeffs=ifit.coef
               ofit=ifit(yarr)
               omean, omed, osigma=saltstat.iterstat((odata-ofit), sig=3, niter=5)
           except ValueError:
               #catch the error if it is a zero array
               ofit=np.array(yarr)*0.0
               osigma=0.0
           except TypeError:
               #catch the error if it is a zero array
               ofit=np.array(yarr)*0.0
               osigma=0.0

           #if it hasn't been already, convert image to
           #double format
           struct[i].data = 1.0 * struct[i].data
           try:
               struct[i].header.remove('BZERO')
               struct[i].header.remove('BSCALE')
           except:
               pass


           #subtract the overscan region
           for j in range(len(struct[i].data[0])):
               struct[i].data[y1:y2,j] -= ofit

           #report the information 
           if log:
                message = '%25s[%1d] %8.2f %3d %7.2f %3d' % \
                    (infile, i, olevel, order, osigma, niter)
                log.message(message, with_stdout=verbose, with_header=False)

           #add the statistics to the image header
           saltkey.new('OVERRMS','%f' % (osigma),'Overscan RMS value', struct[i])

           #update the variance frame
           if saltkey.found('VAREXT', struct[i]):
               vhdu=saltkey.get('VAREXT', struct[i])
               try:
                   vdata=struct[vhdu].data
                   #The bias level should not be included in the noise from the signal
                   for j in range(len(struct[i].data[0])):
		       vdata[y1:y2,j] -= ofit
                   #add a bit to make sure that the minimum error is the rednoise
                   rdnoise= saltkey.get('RDNOISE',struct[i])
                   vdata[vdata<rdnoise**2]=rdnoise**2
                   struct[vhdu].data=vdata+osigma**2

               except Exception, e:
                    msg='Cannot update the variance frame in %s[%i] because %s' % (infile, vhdu, e)
                    raise SaltError(msg)


           #plot the overscan region
           if plotover:  
              plt.plot(yarr, odata)
              plt.plot(yarr, ofit)

       #trim the data and update the headers
       if trim:
           struct[i].data=struct[i].data[dy1:dy2,dx1:dx2]
           datasec = '[1:'+str(dx2-dx1)+',1:'+str(dy2-dy1)+']'
           saltkey.put('DATASEC',datasec,struct[i])

           #update the variance frame
           if saltkey.found('VAREXT', struct[i]):
               vhdu=saltkey.get('VAREXT', struct[i])
               struct[vhdu].data=struct[vhdu].data[dy1:dy2,dx1:dx2]
               datasec = '[1:'+str(dx2-dx1)+',1:'+str(dy2-dy1)+']'
               saltkey.put('DATASEC',datasec,struct[vhdu])
           #update the BPM frame
           if saltkey.found('BPMEXT', struct[i]):
               bhdu=saltkey.get('BPMEXT', struct[i])
               struct[bhdu].data=struct[bhdu].data[dy1:dy2,dx1:dx2]
               datasec = '[1:'+str(dx2-dx1)+',1:'+str(dy2-dy1)+']'
               saltkey.put('DATASEC',datasec,struct[bhdu])

       #subtract the master bias if necessary
       if subbias and bstruct:
           struct[i].data -= bstruct[i].data

           #update the variance frame
           if saltkey.found('VAREXT', struct[i]):
               vhdu=saltkey.get('VAREXT', struct[i])
               try:
                   vdata=struct[vhdu].data
                   struct[vhdu].data=vdata+bstruct[vhdu].data
               except Exception, e:
                    msg='Cannot update the variance frame in %s[%i] because %s' % (infile, vhdu, e)
                    raise SaltError(msg)
 
       

   if plotover: 
       plt.ioff()
       plt.show()

   return struct

# -----------------------------------------------------------
# main code

if not iraf.deftask('saltbias'):
  parfile = iraf.osfn("saltred$saltbias.par")
  t = iraf.IrafTaskFactory(taskname="saltbias",value=parfile,function=saltbias, pkgname='saltred')
