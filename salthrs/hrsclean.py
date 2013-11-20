################################# LICENSE ##################################
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.  See LICENSE for more details                       #
############################################################################

#!/usr/bin/env python

# Author                 Version      Date         Comment
# -----------------------------------------------------------------------
# S M Crawford (SAAO)    0.3          11 Oct 2013  

# hrsclean converts the MEF HRS files and cleans them together 
# into a single image 


import os, glob, time
import pyfits
import numpy as np
from pyraf import iraf
from pyraf.iraf import pysalt

import saltsafestring as saltstring
import saltsafeio as saltio
import saltsafekey as saltkey
from saltsafelog import logging, history

from saltobslog import obslog,createobslogfits
from saltprepare import prepare
from saltgain import gain
from saltxtalk  import xtalk
from saltbias import bias
from saltflat import flat
from saltcombine import saltcombine
from saltclean import compareimages

from hrsprepare import prepare as hrsprepare
from hrsstack import stack

from salterror import SaltError

debug=True

# -----------------------------------------------------------
# core routine
hrsbiasheader_list=['INSTRUME', 'DETNAM', 'DETMODE', 'CCDSUM', 'GAINSET', 'ROSPEED'] 
hrsflatheader_list=['INSTRUME', 'DETNAM', 'DETMODE', 'CCDSUM', 'GAINSET', 'ROSPEED']

def hrsclean(images, outpath, obslogfile=None, subover=True, trim=True, masbias=None, 
             subbias=True, median=False, function='polynomial', order=5, rej_lo=3, rej_hi=3,
             niter=5, interp='linear', clobber=False,  logfile='salt.log',verbose=True):
    """Convert MEF HRS data into a single image.  If variance frames and BPMs, then 
       convert them to the same format as well.   Returns an MEF image but that is
       combined into a single frame

    """
 
    with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outpath=saltio.abspath(outpath)

       if saltio.checkfornone(obslogfile) is None:
          raise SaltError('Obslog file is required')

       # Delete the obslog file if it already exists 
          

       if (os.path.isfile(obslogfile) and clobber) or not os.path.isfile(obslogfile): 
          if os.path.isfile(obslogfile): saltio.delete(obslogfile)
          #read in the obsveration log or create it
          headerDict=obslog(infiles, log)
          obsstruct=createobslogfits(headerDict)
          saltio.writefits(obsstruct, obslogfile)
       else:
          obsstruct=saltio.openfits(obslogfile)

       #create the list of bias frames and process them
       filename=obsstruct.data.field('FILENAME')
       detmode=obsstruct.data.field('DETMODE')
       ccdtype=obsstruct.data.field('OBJECT')

       biaslist=filename[ccdtype=='Bias']
       masterbias_dict={}
       if log: log.message('Processing Bias Frames')
       for img in infiles:
           if os.path.basename(img) in biaslist:
               #open the image
               struct=pyfits.open(img)
               bimg=outpath+'bgph'+os.path.basename(img)
               #print the message
               if log:
                   message='Processing Zero frame %s' % img
                   log.message(message, with_stdout=verbose, with_header=False)

               #process the image
               #struct=clean(struct, createvar=False, badpixelstruct=None, mult=True,
               #             subover=subover, trim=trim, subbias=False, imstack=False,
               #             bstruct=None, median=median, function=function, order=order,
               #             rej_lo=rej_lo, rej_hi=rej_hi, niter=niter, log=log,
               #             verbose=verbose)

               #write the file out
               # housekeeping keywords
               fname, hist=history(level=1, wrap=False, exclude=['images', 'outimages', 'outpref'])
               saltkey.housekeeping(struct[0],'HPREPARE', 'Images have been prepared', hist)
               saltkey.new('HGAIN',time.asctime(time.localtime()),'Images have been gain corrected',struct[0])
               #saltkey.new('HXTALK',time.asctime(time.localtime()),'Images have been xtalk corrected',struct[0])
               saltkey.new('HBIAS',time.asctime(time.localtime()),'Images have been de-biased',struct[0])

               # write FITS file
              # saltio.writefits(struct,bimg, clobber=clobber)
               saltio.closefits(struct)
 
               #add files to the master bias list
               masterbias_dict=compareimages(struct, bimg, masterbias_dict, keylist=hrsbiasheader_list)

       #create the master bias frame
       for i in masterbias_dict.keys():
           bkeys=masterbias_dict[i][0]
           blist=masterbias_dict[i][1:]
           mbiasname=outpath+createmasterbiasname(blist, bkeys, x1=5, x2=13)
           bfiles=','.join(blist)
           #saltcombine(bfiles, mbiasname, method='median', reject='sigclip', mask=False,
           #            weight=False, blank=0, scale=None, statsec=None, lthresh=3,    \
           #            hthresh=3, clobber=False, logfile=logfile,verbose=verbose)

       #apply full reductions to the science data
       for img in infiles:
           nimg=os.path.basename(img)
           if not nimg in biaslist:
               #open the image
               struct=pyfits.open(img)
               simg=outpath+'mbgph'+os.path.basename(img)

               #print the message
               if log:
                   message='Processing science frame %s' % img
                   log.message(message, with_stdout=verbose)

               #get master bias frame
               masterbias=get_masterbias(struct, masterbias_dict, keylist=hrsbiasheader_list) 
               if masterbias:
                  subbias=True
                  bstruct=saltio.openfits(masterbias)
               else:
                  subbias=False
                  bstruct=None

               #process the image
               struct=clean(struct, createvar=False, badpixelstruct=None, mult=True,
                            subover=subover, trim=trim, subbias=subbias, imstack=True,
                            bstruct=bstruct, median=median, function=function, order=order,
                            rej_lo=rej_lo, rej_hi=rej_hi, niter=niter, log=log,
                            verbose=verbose)

               #write the file out
               # housekeeping keywords
               fname, hist=history(level=1, wrap=False, exclude=['images', 'outimages', 'outpref'])
               saltkey.housekeeping(struct[0],'HPREPARE', 'Images have been prepared', hist)
               saltkey.new('HGAIN',time.asctime(time.localtime()),'Images have been gain corrected',struct[0])
               #saltkey.new('HXTALK',time.asctime(time.localtime()),'Images have been xtalk corrected',struct[0])
               saltkey.new('HBIAS',time.asctime(time.localtime()),'Images have been de-biased',struct[0])

               # write FITS file
               saltio.writefits(struct,simg, clobber=clobber)
               saltio.closefits(struct)


 
    return 

def get_masterbias(struct, imdict, keylist):
   """Return the name of the master bias frame if one is available
   """
   #create the list of header parameters
   klist=[]
   for k in keylist:
       try:
           value=str(struct[0].header[k]).strip()
       except:
           value=''
       klist.append(value)
   biasname=createmasterbiasname([os.path.basename(struct[0]._file.name)], klist, x1=1, x2=9)

   if os.path.isfile(biasname): return biasname
   return None


def clean(struct, createvar=False, badpixelstruct=None, mult=True,            \
          subover=True, trim=True, subbias=False, bstruct=None, imstack=False,  
          median=False, \
          function='polynomial', order=5, rej_lo=3, rej_hi=3, niter=5,        \
          log=None, verbose=True):
    """Clean HRS data and files.  This includes the following steps currently:
             * overscan correction
             * bias subtraction
             * gain correction
             * mosaic correction 
    """
   
    infile=struct

    tfile=struct[0]._file
    #prepare the HRS files
    struct=hrsprepare(struct)
    struct[0]._file=tfile
    ampccd=len(struct)-1

    #prepare the files
    struct=prepare(struct, createvar=createvar, badpixelstruct=badpixelstruct)

    #gain correct the files
    usedb=False
    struct=gain(struct, mult=mult,usedb=usedb, dblist=None, ampccd=ampccd, log=log, verbose=verbose)

    #overscan and bias correction 
    struct=bias(struct,subover=subover, trim=trim, subbias=subbias,
               bstruct=bstruct, median=median, function=function,
               order=order, rej_lo=rej_lo, rej_hi=rej_hi, niter=niter,
               plotover=False, log=log, verbose=verbose)

    #stack the data if requested
    if imstack:
       struct=stack(struct)
       struct=salt2iraf(struct)

    return struct

def salt2iraf(hdu, ext=1):
   """Convert a SALT MEF file into a single extension fits file"""

   #create the new output image
   odu = pyfits.PrimaryHDU(data=hdu[ext].data, header=hdu[0].header)

   #combine the headers from the primary and first exention
   for c in hdu[ext].header.ascardlist(): odu.header.update(c.key, c.value, c.comment)

   return hdu

 
def createmasterbiasname(infiles, biaskeys, x1=5, x2=13):
    """Create the name for the master bias file based on its parameters.  The format for 
       hte name is 
       [S/P][YYYYMMDD]Bias[DETNAM][BINNING][GAINSET][ROSPEED].fits
    
       where the following abbreviations are used:
      
       [H]--Scam or RSS
       [YYYYMMDD]--obsdate of the data or most common obsdate if multiple dates
       [DETNAM]--Detector Name
               Blue: B
               Red: R
       [BINNING]--CCD binning in XBINxYBIN
       [AMPS]--Number of Amps
                  Bright: BR
                  Faint: FA

       Parameters: 
       x1: Place to start to extract obsdate
       x2: Place to stop to extract obsdate
    """
    #setup the in the instrument
   
    #setup the obsdate--assumes fixed naming scheme
    obsdate=saltstring.makeobsdatestr(infiles, x1=x1, x2=x2)

    #if len(obsdate)<4: obsdate=''

    #set the mode string
    #mdstr=saltstring.makedetmodestr(biaskeys[2])
    detstr='B'
    instr='H'
    ampstr='A2'
    if biaskeys[1]=='08443-03-01': 
       detstr='R'
       instr='R'
       ampstr='A4'

    #set binning
    binstr=saltstring.makebinstr(biaskeys[3])

    #set the AMP string


    biasname='%s%sBias%s%s%s.fits' % (instr, obsdate, detstr, binstr, ampstr)
    return biasname

   



if not iraf.deftask('hrsclean'):
    parfile = iraf.osfn("salthrs$hrsclean.par")
    t = iraf.IrafTaskFactory(taskname="hrsclean",value=parfile,function=hrsclean, pkgname='salthrs')

