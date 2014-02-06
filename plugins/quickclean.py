############################### LICENSE ##################################
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved. See License file for more details                   #
#                                                                          #
############################################################################


#!/usr/bin/env python

"""
QUICKCLEAN -- QUICKCLEAN is a plugin for saltfirst that provides quick 
reductions for normal imaging (ie not slotmode).  

Author                 Version      Date
-----------------------------------------------
S M Crawford (SAAO)    0.1          16 Mar 2010

"""
import os 

#pysalt imports
from pyraf import iraf
from iraf import pysalt

from pyraf.iraf import saltred

from saltprepare import prepare
from saltgain import gain
from saltxtalk  import xtalk
from saltbias import bias
from saltflat import flat
from saltcrclean import multicrclean


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
   struct=saltio.openfits(rawpath+'/'+infile)
   struct.verify('exception')
   
   #check to see if detmode is there
   if not saltkey.found('DETMODE', struct[0]): 
      return 
 
   #reduce the file
   struct=prepare(struct, createvar=False, badpixelstruct=None)
 
      #reset the names in the structures
   for i in range(1,len(struct)):
       struct[i].name=struct[i].header['EXTNAME']


   #gain correct the files
   usedb=True
   dblist= saltio.readgaindb(gaindb)
   log=open(logfile, 'a')
   ampccd = struct[0].header['NAMPS'] / struct[0].header['NCCDS']
   struct=gain(struct, mult=True,usedb=usedb, dblist=dblist, ampccd=ampccd, log=None, verbose=verbose)

   struct=bias(struct, subover=True,trim=True,subbias=False, 
                    median=False,function='polynomial',order=5,rej_lo=3,rej_hi=3,niter=10,
                    plotover=False,log=None, verbose=verbose)

   if struct[0].header['CCDTYPE']=='OBJECT' and struct[0].header['EXPTIME']>90:
      struct = multicrclean(struct, crtype='median', thresh=5, mbox=5, bbox=25, bthresh=5, flux_ratio=0.2, \
                          gain=1, rdnoise=5, bfactor=2, fthresh=5, gbox=0, maxiter=5, log=None, verbose=verbose)

   pinfile=outpath+'bxp'+infile
   saltio.writefits(struct, pinfile, clobber)

   saltred.saltmosaic(images=pinfile,
                   outimages='',outpref=outpath+'m',geomfile=geomfile,
                   interp=interp,cleanup=cleanup,clobber=clobber,logfile=logfile,
                   verbose=verbose)
   profile=outpath+'mbxp'+infile

   #remove intermediate steps
   if cleanup:
      if os.path.isfile(pinfile): os.remove(pinfile)

   return
