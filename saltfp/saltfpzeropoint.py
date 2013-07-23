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

"""SALTFPZEROPOINT is a tool to add Fabry Perot coefficients to the data. The
   task reads in a calibration file and calculates the zeropoint shift for 
   each image.  It then updates the headers in the image with the FP
   coefficients

Updates:

20130115
    * First wrote the code
"""

# Ensure python 2.5 compatibility
from __future__ import with_statement

import os
import sys
import numpy as np
import pyfits
from datetime import datetime

from pyraf import iraf
from pyraf.iraf import pysalt
import saltsafekey as saltkey
import saltsafeio as saltio
import salttime
from saltsafelog import logging
from salterror import SaltError

from fptools import fpfunc


debug=True

def saltfpzeropoint(images,outimages, outpref, calfile, fpa, fpb, fpc, fpd, fpe, fpf, clobber=False,logfile='salt.log',verbose=True):  
   """Adds zeropoint information to each image"""

   with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create the output files
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       #verify that the input and output lists are the same length
       saltio.comparelists(infiles,outfiles,'Input','output')


       #calculate the zeropoint coefficients
       fpcoef=np.array([fpa,fpb,fpc,fpd,fpe,fpf])
       zpcoef,tstart=calc_zpcoef(calfile, fpcoef)

       # open each image and detect the ring
       for img, oimg  in zip(infiles, outfiles):
      
          hdu=saltio.openfits(img)

          #get the image time
          t=(get_datetime(hdu)-tstart).seconds

          #add the header values to the image
          saltkey.new('FPA',fpa+zpcoef[1]+zpcoef[0]*t,'FPA Coef',hdu[0])
          saltkey.new('FPB',fpb,'FPB Coef',hdu[0])
          saltkey.new('FPC',fpc,'FPC Coef',hdu[0])
          saltkey.new('FPD',fpd,'FPD Coef',hdu[0])
          saltkey.new('FPE',  0,'FPE Coef',hdu[0])
          saltkey.new('FPF',fpf,'FPF Coef',hdu[0])

          #write the file out
          saltio.writefits(hdu, oimg, clobber)

          #log the action
          msg='Updating the FP information in %s' % (img)
          log.message(msg, with_stdout=verbose, with_header=False)


def calc_zpcoef(calfile, fpcoef):    
   """Given the values given in the calibraiton file, calculate the
      offset to the FPA value and FPE value that would give an 
      appropriate calibration values to the frame

      The task returns the zeropoint coefficient and the initial
      start time

   """
   #read in the values for the calfile--assumes the format for the file
   #r,r_err,z,t,w,img=np.loadtxt(calfile, usecols=(0,1,4,5,6,8), unpack=True, dtype={8,str})
   data=np.loadtxt(calfile, dtype={'names': ('r','r_err', 'x', 'y', 'z', 't', 'w', 'dn', 'image'),'formats':('f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4','S30')})
   #set the start time and convert the time to a date time unit
   time_list=[]
   time_start=None

   if data.size==0:
      raise SaltError('%s needs to have at least one entry' % calfile)
   elif data.size==1:
      img_list=[str(data['image'])]
   else:
      img_list=data['image']

   for img in img_list:
      if not os.path.isfile(img): 
         raise SaltError('%s needs to be available to open' % img)
      t=get_datetime(saltio.openfits(img))
      if time_start is None: time_start=t
      time_list.append((t-time_start).seconds)
   time_arr=np.array(time_list)

   #calculate the coefficients
   wf=fpfunc(data['z'], data['r'], time_arr, coef=fpcoef)
   if data.size==1:
      coef=np.array([0, data['w']-wf])
   else:
      coef=np.polyfit(time_arr, data['w']-wf, 1)
   
   return coef, time_start

def get_datetime(hdu):
    """Determine the datetime of an observation"""
    d=hdu[0].header['DATE-OBS']
    t=hdu[0].header['TIME-OBS']
    return datetime.strptime('%s %s' % (d,t), '%Y-%m-%d %H:%M:%S.%f')

# -----------------------------------------------------------
# main code

parfile = iraf.osfn("saltfp$saltfpzeropoint.par")
t = iraf.IrafTaskFactory(taskname="saltfpzeropoint",value=parfile,function=saltfpzeropoint,pkgname='saltfp')
