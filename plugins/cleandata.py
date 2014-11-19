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
QUICKCLEAN -- QUICKCLEAN is a plugin for saltfirst that provides quick 
reductions for normal imaging (ie not slotmode).   The tasks 

Author                 Version      Date
-----------------------------------------------
S M Crawford (SAAO)    0.1          16 Mar 2010

"""

import os, shutil, time, ftplib, glob, pyfits
import numpy as np
import scipy as sp

#pysalt imports
from pyraf import iraf
from pyraf.iraf import pysalt
from pyraf.iraf import saltred
from pyraf.iraf import saltspec
#from pyraf.iraf import pipetools

import saltsafekey as saltkey
import saltsafeio as saltio

def cleandata(filename, iminfo=None, prodir='.', interp='linear', cleanup=True, clobber=False, logfile='saltclean.log', display_image=False, verbose=True):
   """Start the process to reduce the data and produce a single mosaicked image"""
   print filename
   status=0
   #create the input file name
   infile=os.path.basename(filename)
   rawpath=os.path.dirname(filename)
   outpath='./'
   outfile=outpath+'mbxp'+infile
   print infile, rawpath, outpath

   #check to see if the data have detmode
   if iminfo is not None:
      detmode=iminfo[headerList.index('DETMODE')].strip().upper()
      print 'DETMODE:' detmode

   #If it is a bin file, pre-process the data
   if filename.count('.bin'):
       print "I can't handle this yet"


   #check to see if it exists and return if clobber is no
   if os.path.isfile(outfile) and not clobber: return

   #set up the files needed
   if infile[0]=='P':
     gaindb = '/iraf/extern/pysalt/data/rss/RSSamps.dat'
     xtalkfile = '/iraf/extern/pysalt/data/rss/RSSxtalk.dat'
     geomfile = '/iraf/extern/pysalt/data/rss/RSSgeom.dat'
     usedb=True
   elif infile[0]=='S':
     gaindb = '/iraf/extern/pysalt/data/scam/SALTICAMamps.dat'
     xtalkfile = '/iraf/extern/pysalt/data/scam/SALTICAMxtalk.dat'
     geomfile = '/iraf/extern/pysalt/data/scam/SALTICAMgeom.dat'
 
   #verify the file
   hdu=saltio.openfits(rawpath+'/'+infile)
   hdu.verify('exception')
   

   #reduce the file
   saltred.saltprepare(images=infile,rawpath=rawpath,outimages='',outpref=outpath+'p',  \
                    clobber=clobber,logfile=logfile,verbose=verbose,status=status)
   pinfile=outpath+'p'+infile
   saltred.saltslot(images=pinfile,outimages='',outpref=outpath+'bx',gaindb=gaindb,
                 xtalkfile=xtalkfile,usedb=True, clobber=clobber,logfile=logfile,verbose=verbose,
                 status=status)
   biasfile=outpath+'bxp'+infile
   saltred.saltmosaic(images=biasfile,
                   outimages='',outpref=outpath+'m',geomfile=geomfile,
                   interp=interp,cleanup=cleanup,fill=True, clobber=clobber,logfile=logfile,
                   verbose=verbose, status=status)
   profile=outpath+'mbxp'+infile

   #remove intermediate steps
   if cleanup:
      if os.path.isfile(pinfile): os.remove(pinfile)
      if os.path.isfile(biasfile): os.remove(biasfile)

       i=headerList.index('CCDSUM')
       ccdbin=int(iminfo[i].split()[0])
       pix_scale=0.14
       r_ap=1.5/(pix_scale*ccdbin)
       print pix_scale, ccdbin, r_ap

       profile=outpath+'mbxp'+infile
       outcat=profile.split('.fits')[0]+'.cat'
       sexfile='/home/ccd/tools/qred.sex'
       backfile = profile.strip().strip('.fits')+'_back.fits'
       cmd='sex %s -c %s -PIXEL_SCALE %f -CATALOG_NAME %s -PHOT_APERTURES %f ' % (profile.strip(),sexfile, pix_scale,outcat,r_ap)
       print "SALTFIRST--Performing photometry on %s" % profile
       if os.path.isfile(sexfile): os.system(cmd)


   #If the images are spectral images, run specreduce on them
   if obsmode=='SPECTROSCOPY' and not(target in ['FLAT', 'BIAS']):
       try:
           profile=outpath+'mbxp'+infile
           badpixelimage=None
           caltype='rss'
           function='polynomial'
           order=3
           skysub=True
           if lampid: skysub=False
           skysection='1400:1500'
           findobj=False
           objsection='900:1100'
           thresh=3.0
           logfile='saltspec.log'
           saltspec.specreduce(images=profile,badpixelimage=badpixelimage,caltype=caltype,
                function=function, order=order, skysub=skysub, skysection=skysection,
                findobj=findobj, objsection=objsection, thresh=thresh, 
                clobber=clobber, logfile=logfile, verbose=verbose)
       except Exception,e:
           message="SALTFIRST--ERROR:  Could not wavelength calibrate %s because %s" % (infile, e)
           fout=open(logfile, 'a')
           fout.write(message)
           print message


   return
