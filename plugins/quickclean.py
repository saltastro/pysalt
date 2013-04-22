############################### LICENSE ##################################
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
import os 

#pysalt imports
from pyraf import iraf
from pyraf.iraf import pysalt
from pyraf.iraf import saltred

from saltcrclean import saltcrclean

import saltsafeio as saltio
import saltsafekey as saltkey

def quickclean(filename, interp='linear', cleanup=True, clobber=False, logfile='saltclean.log', verbose=True):
   """Start the process to reduce the data and produce a single mosaicked image"""
   print filename

   #create the input file name
   status=0
   infile=os.path.basename(filename)
   rawpath=os.path.dirname(filename)
   outpath='./'
   outfile=outpath+'mbxp'+infile
   print infile, rawpath, outpath

   #check to see if it exists and return if clobber is no
   if os.path.isfile(outfile) and not clobber: return

   #set up the files needed
   if infile[0]=='P':
     gaindb = iraf.osfn('pysalt$data/rss/RSSamps.dat')
     xtalkfile = iraf.osfn('pysalt$data/rss/RSSxtalk.dat')
     geomfile = iraf.osfn('pysalt$data/rss/RSSgeom.dat')
   elif infile[0]=='S':
     gaindb = iraf.osfn('pysalt$data/scam/SALTICAMamps.dat')
     xtalkfile = iraf.osfn('pysalt$data/scam/SALTICAMxtalk.dat')
     geomfile = iraf.osfn('pysalt$data/scam/SALTICAMgeom.dat')
 
   #verify the file
   hdu=saltio.openfits(rawpath+'/'+infile)
   hdu.verify('exception')
   
   #check to see if detmode is there
   if not saltkey.found('DETMODE', hdu[0]): 
      return 
 
   #reduce the file
   saltred.saltprepare(images=filename,outimages='',outpref=outpath+'p',  \
                    createvar=False, badpixelimage=None, clobber=clobber,logfile=logfile,verbose=verbose)
   pinfile=outpath+'p'+infile
   saltred.saltgain(pinfile, outimages=pinfile, outpref='', gaindb=gaindb,usedb=False, 
                    mult=True,clobber=True, logfile=logfile, verbose=verbose)
   saltred.saltxtalk(pinfile,outimages='',outpref='x',xtalkfile=xtalkfile,clobber=clobber,
                     logfile=logfile,verbose=verbose)
   #saltred.saltslot(images=pinfile,outimages='',outpref=outpath+'bx',gaindb=gaindb,
   #              xtalkfile=xtalkfile,clobber=clobber,logfile=logfile,verbose=verbose,
   #              status=0)
   xinfile=outpath+'xp'+infile
   saltred.saltbias(images=xinfile,outimages='',outpref='b',subover=True,trim=True,subbias=False, 
                    masterbias='', median=False,function='polynomial',order=5,rej_lo=3,rej_hi=3,niter=10,
                    plotover=False,turbo=False,logfile=logfile, clobber=clobber, verbose=verbose)
   biasfile=outpath+'bxp'+infile

   if hdu[0].header['CCDTYPE']=='OBJECT' and hdu[0].header['EXPTIME']>90:
       saltcrclean(images=biasfile, outimages=biasfile, outpref='', crtype='median',thresh=5,mbox=5,         \
                bthresh=3, flux_ratio=0.2, bbox=25, gain=1, rdnoise=5, fthresh=5,\
                bfactor=2, gbox=0, maxiter=5, multithread=True, clobber=True,          \
                logfile='salt.log', verbose=True)

   saltred.saltmosaic(images=biasfile,
                   outimages='',outpref=outpath+'m',geomfile=geomfile,
                   interp=interp,cleanup=cleanup,clobber=clobber,logfile=logfile,
                   verbose=verbose)
   profile=outpath+'mbxp'+infile

   #remove intermediate steps
   if cleanup:
      if os.path.isfile(pinfile): os.remove(pinfile)
      if os.path.isfile(xinfile): os.remove(xinfile)
      if os.path.isfile(biasfile): os.remove(biasfile)

   return
