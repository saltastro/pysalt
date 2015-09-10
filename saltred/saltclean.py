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
SALTCLEAN performs the standard chain of primary reduction tasks on
lists of RSS or SALTICAM images.

 Author                 Version      Date
 -----------------------------------------------
 Martin Still (SAAO)    0.2          31 Jul 2007
 S M Crawford (SAAO)    0.3          22 Feb 2008
 S M Crawford (SAAO)    0.4          11 Sep 2011

Updates
------------------------------------------------
20110911 -  Major re-write to include new error handling
            and interaction with other programs
20120926 -  Updated so that slotmode always defaults to 
            order 1

Todo
------------------------------------------------

--Allow it to read in obslogfile

"""

from __future__ import with_statement

import sys,glob, os, shutil, time
import numpy as np
from astropy.io import fits

from pyraf import iraf
from pyraf.iraf import pysalt

from saltobslog import obslog, createobslogfits
from saltprepare import prepare
from saltgain import gain
from saltxtalk  import xtalk
from saltbias import bias
from saltflat import flat
from saltcombine import saltcombine 
from saltmosaic import saltmosaic

import saltsafestring as saltstring
import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging, history

from salterror import SaltError

debug=True

biasheader_list=['INSTRUME', 'DETMODE', 'CCDSUM', 'GAINSET', 'ROSPEED', 'NWINDOW']
flatheader_list=['INSTRUME', 'DETMODE', 'CCDSUM', 'GAINSET', 'ROSPEED', 'FILTER', 'GRATING', 'GR-ANGLE', 'AR-ANGLE', 'NWINDOW']


# -----------------------------------------------------------
# core routine

def saltclean(images, outpath, obslogfile=None, gaindb=None,xtalkfile=None, 
	geomfile=None,subover=True,trim=True,masbias=None, 
        subbias=False, median=False, function='polynomial', order=5,rej_lo=3,
        rej_hi=3,niter=5,interp='linear', clobber=False, logfile='salt.log', 
        verbose=True):
   """SALTCLEAN will provide basic CCD reductions for a set of data.  It will 
      sort the data, and first process the biases, flats, and then the science 
      frames.  It will record basic quality control information about each of 
      the steps.
   """
   plotover=False

   #start logging
   with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outpath=saltio.abspath(outpath)


       #does the gain database file exist
       if gaindb:
           dblist= saltio.readgaindb(gaindb)
       else:
           dblist=[]

       # does crosstalk coefficient data exist
       if xtalkfile:
           xtalkfile = xtalkfile.strip()
           xdict = saltio.readxtalkcoeff(xtalkfile)
       else:
           xdict=None
       #does the mosaic file exist--raise error if no
       saltio.fileexists(geomfile)


       # Delete the obslog file if it already exists
       if os.path.isfile(obslogfile) and clobber: saltio.delete(obslogfile)

       #read in the obsveration log or create it
       if os.path.isfile(obslogfile):
           msg='The observing log already exists.  Please either delete it or run saltclean with clobber=yes'
           raise SaltError(msg)
       else:
           headerDict=obslog(infiles, log)
           obsstruct=createobslogfits(headerDict)
           saltio.writefits(obsstruct, obslogfile)

       #create the list of bias frames and process them
       filename=obsstruct.data.field('FILENAME')
       detmode=obsstruct.data.field('DETMODE')
       ccdtype=obsstruct.data.field('CCDTYPE')

       #set the bias list of objects
       biaslist=filename[ccdtype=='ZERO']
       masterbias_dict={}
       for img in infiles:
           if os.path.basename(img) in biaslist:
               #open the image
               struct=fits.open(img)
               bimg=outpath+'bxgp'+os.path.basename(img)

               #print the message
               if log:
                   message='Processing Zero frame %s' % img
                   log.message(message, with_stdout=verbose)

               #process the image
               struct=clean(struct, createvar=False, badpixelstruct=None, mult=True, 
                            dblist=dblist, xdict=xdict, subover=subover, trim=trim, subbias=False,
                            bstruct=None, median=median, function=function, order=order,
                            rej_lo=rej_lo, rej_hi=rej_hi, niter=niter, plotover=plotover, log=log,
                            verbose=verbose)

               #write the file out
               # housekeeping keywords
               fname, hist=history(level=1, wrap=False, exclude=['images', 'outimages', 'outpref'])
               saltkey.housekeeping(struct[0],'SPREPARE', 'Images have been prepared', hist)
               saltkey.new('SGAIN',time.asctime(time.localtime()),'Images have been gain corrected',struct[0])
               saltkey.new('SXTALK',time.asctime(time.localtime()),'Images have been xtalk corrected',struct[0])
               saltkey.new('SBIAS',time.asctime(time.localtime()),'Images have been de-biased',struct[0])

               # write FITS file
               saltio.writefits(struct,bimg, clobber=clobber)
               saltio.closefits(struct)

               #add files to the master bias list
               masterbias_dict=compareimages(struct, bimg, masterbias_dict, keylist=biasheader_list)

       #create the master bias frame
       for i in masterbias_dict.keys():
           bkeys=masterbias_dict[i][0]
           blist=masterbias_dict[i][1:]
           mbiasname=outpath+createmasterbiasname(blist, bkeys)
           bfiles=','.join(blist)
           saltcombine(bfiles, mbiasname, method='median', reject='sigclip', mask=False, 
                       weight=False, blank=0, scale=None, statsec=None, lthresh=3,    \
                       hthresh=3, clobber=False, logfile=logfile,verbose=verbose)

           

       #create the list of flatfields and process them
       flatlist=filename[ccdtype=='FLAT']
       masterflat_dict={}
       for img in infiles:
           if os.path.basename(img) in flatlist:
               #open the image
               struct=fits.open(img)
               fimg=outpath+'bxgp'+os.path.basename(img)

               #print the message
               if log:
                   message='Processing Flat frame %s' % img
                   log.message(message, with_stdout=verbose)

               #process the image
               struct=clean(struct, createvar=False, badpixelstruct=None, mult=True, 
                            dblist=dblist, xdict=xdict, subover=subover, trim=trim, subbias=False,
                            bstruct=None, median=median, function=function, order=order,
                            rej_lo=rej_lo, rej_hi=rej_hi, niter=niter, plotover=plotover, log=log,
                            verbose=verbose)

               #write the file out
               # housekeeping keywords
               fname, hist=history(level=1, wrap=False, exclude=['images', 'outimages', 'outpref'])
               saltkey.housekeeping(struct[0],'SPREPARE', 'Images have been prepared', hist)
               saltkey.new('SGAIN',time.asctime(time.localtime()),'Images have been gain corrected',struct[0])
               saltkey.new('SXTALK',time.asctime(time.localtime()),'Images have been xtalk corrected',struct[0])
               saltkey.new('SBIAS',time.asctime(time.localtime()),'Images have been de-biased',struct[0])

               # write FITS file
               saltio.writefits(struct,fimg, clobber=clobber)
               saltio.closefits(struct)

               #add files to the master bias list
               masterflat_dict=compareimages(struct, fimg, masterflat_dict,  keylist=flatheader_list)

       #create the master flat frame
       for i in masterflat_dict.keys():
           fkeys=masterflat_dict[i][0]
           flist=masterflat_dict[i][1:]
           mflatname=outpath+createmasterflatname(flist, fkeys)
           ffiles=','.join(flist)
           saltcombine(ffiles, mflatname, method='median', reject='sigclip', mask=False, 
                       weight=False, blank=0, scale=None, statsec=None, lthresh=3,    \
                       hthresh=3, clobber=False, logfile=logfile,verbose=verbose)

       #process the science data
       for img in infiles:
           nimg=os.path.basename(img)
           #print nimg, nimg in flatlist, nimg in biaslist
           if not (nimg in biaslist):
               #open the image
               struct=fits.open(img)
               simg=outpath+'bxgp'+os.path.basename(img)

               #print the message
               if log:
                   message='Processing science frame %s' % img
                   log.message(message, with_stdout=verbose)

               #process the image
               struct=clean(struct, createvar=False, badpixelstruct=None, mult=True, 
                            dblist=dblist, xdict=xdict, subover=subover, trim=trim, subbias=False,
                            bstruct=None, median=median, function=function, order=order,
                            rej_lo=rej_lo, rej_hi=rej_hi, niter=niter, plotover=plotover, log=log,
                            verbose=verbose)

               #write the file out
               # housekeeping keywords
               fname, hist=history(level=1, wrap=False, exclude=['images', 'outimages', 'outpref'])
               saltkey.housekeeping(struct[0],'SPREPARE', 'Images have been prepared', hist)
               saltkey.new('SGAIN',time.asctime(time.localtime()),'Images have been gain corrected',struct[0])
               saltkey.new('SXTALK',time.asctime(time.localtime()),'Images have been xtalk corrected',struct[0])
               saltkey.new('SBIAS',time.asctime(time.localtime()),'Images have been de-biased',struct[0])

               # write FITS file
               saltio.writefits(struct,simg, clobber=clobber)
               saltio.closefits(struct)

               #mosaic the files--currently not in the proper format--will update when it is
               if not saltkey.fastmode(saltkey.get('DETMODE', struct[0])):
                   mimg=outpath+'mbxgp'+os.path.basename(img)
                   saltmosaic(images=simg, outimages=mimg,outpref='',geomfile=geomfile,
                        interp=interp,cleanup=True,clobber=clobber,logfile=logfile,
                        verbose=verbose)

                   #remove the intermediate steps
                   saltio.delete(simg)




def compareimages(struct, oimg, imdict, keylist):
   """See if the current structure is held in the dictionary of images.  If it is,
       then add it to the list.  If it isn't then create a new entry
   """
   #create the list of header parameters
   klist=[]
   for k in keylist:
       try:
           value=str(struct[0].header[k]).strip()
       except:
           value=''
       klist.append(value)

   if len(imdict)==0: 
       imdict[oimg]=[klist, oimg]
       return imdict

   #compare each value of imdict to the structure
   for i in imdict.keys():
       if klist==imdict[i][0]:
          imdict[i].append(oimg)
          return imdict

   #create a new one if it isn't found
   imdict[oimg]=[klist, oimg]
       
   return imdict
 

def clean(struct, createvar=False, badpixelstruct=None, mult=True, dblist=None, ampccd=2,
          xdict=[], subover=True,trim=True, subbias=False, bstruct=None,
         median=False, function='polynomial',order=3,rej_lo=3,rej_hi=3,niter=10,
         plotover=False, log=None, verbose=True):

   infile=struct

   #prepare the files
   struct=prepare(struct, createvar=createvar, badpixelstruct=badpixelstruct)

   #reset the names in the structures
   for i in range(1,len(struct)):
       struct[i].name=struct[i].header['EXTNAME']


   #gain correct the files
   usedb=False
   if dblist:  usedb=True
   struct=gain(struct, mult=mult,usedb=usedb, dblist=dblist, ampccd=ampccd, log=log, verbose=verbose)

   #xtalk correct the files
   usedb=False
   if xdict:  
       obsdate=saltkey.get('DATE-OBS', struct[0])
       try:
           obsdate=int('%s%s%s' % (obsdate[0:4],obsdate[5:7], obsdate[8:]))
           xkey=np.array(xdict.keys())
           date=xkey[abs(xkey-obsdate).argmin()]
           xcoeff=xdict[date]
       except Exception,e : 
           msg='WARNING--Can not find xtalk coefficient for %s because %s' % (e, infile)
           if log: log.warning(msg)
           xcoeff=xdict[xdict.keys()[-1]]
   else:
       xcoeff=[]
   struct = xtalk(struct, xcoeff, log=log, verbose=verbose)

   #bias correct the files
   if saltkey.fastmode(saltkey.get('DETMODE', struct[0])): order=1

   struct=bias(struct,subover=subover, trim=trim, subbias=subbias,
               bstruct=bstruct, median=median, function=function,
               order=order, rej_lo=rej_lo, rej_hi=rej_hi, niter=niter,
               plotover=plotover, log=log, verbose=verbose)


   #mosaic correct the files


   return struct

def createmasterbiasname(infiles, biaskeys):
    """Create the name for the master bias file based on its parameters.  The format for 
       hte name is 
       [S/P][YYYYMMDD]Bias[MODE][BINNING][GAINSET][ROSPEED].fits
    
       where the following abbreviations are used:
      
       [S/P]--Scam or RSS
       [YYYYMMDD]--obsdate of the data or most common obsdate if multiple dates
       [MODE]--Mode of the observations:
               Normal: NM
               Framte Transfer: FT
               Slot Mode:   SL
               Drift Scanning: DS
       [BINNING]--CCD binning in XBINxYBIN
       [GAINSET]--Gain setting
                  Bright: BR
                  Faint: FA
       [ROSPEED]--Read out speed
                  FAST: FA
                  SLOW: SL
       
    """
    
    #setup the in the instrument
    instr=saltstring.makeinstrumentstr(biaskeys[0])
    
    #setup the obsdate--assumes fixed naming scheme
    obsdate=saltstring.makeobsdatestr(infiles)
        
    #if len(obsdate)<4: obsdate=''
    print obsdate

    #set the mode string
    mdstr=saltstring.makedetmodestr(biaskeys[1])
 
    #set binning
    binstr=saltstring.makebinstr(biaskeys[2])

    #set gain
    gnstr=saltstring.makegainstr(biaskeys[3])

    #set readout
    rostr=saltstring.makereadoutstr(biaskeys[4])

    biasname='%s%sBias%s%s%s%s.fits' % (instr, obsdate, mdstr, binstr, gnstr, rostr)
    return biasname


def createmasterflatname(infiles, flatkeys):
    """Create the name for the master flat file based on its parameters.  The format for 
       hte name is 
       [S/P][YYYYMMDD]Flat[MODE][BINNING][GAINSET][ROSPEED][FILTER].fits
    
       where the following abbreviations are used:
      
       [S/P]--Scam or RSS
       [YYYYMMDD]--obsdate of the data or most common obsdate if multiple dates
       [MODE]--Mode of the observations:
               Normal: NM
               Framte Transfer: FT
               Slot Mode:   SL
               Drift Scanning: DS
       [BINNING]--CCD binning in XBINxYBIN
       [GAINSET]--Gain setting
                  Bright: BR
                  Faint: FA
       [ROSPEED]--Read out speed
                  FAST: FA
                  SLOW: SL
       [FILTER]--Filter used
       
    """
    #setup the in the instrument
    instr=saltstring.makeinstrumentstr(flatkeys[0])
    
    #setup the obsdate--assumes fixed naming scheme
    obsdate=saltstring.makeobsdatestr(infiles)
        
    #if len(obsdate)<4: obsdate=''
    print obsdate

    #set the mode string
    mdstr=saltstring.makedetmodestr(flatkeys[1])
 
    #set binning
    binstr=saltstring.makebinstr(flatkeys[2])

    #set gain
    gnstr=saltstring.makegainstr(flatkeys[3])

    #set readout
    rostr=saltstring.makereadoutstr(flatkeys[4])
    
    fltstr=flatkeys[5].strip()

    if flatkeys[6].count('SKY'): 
       skystr='Sky'
    else:
       skystr=''

    flatname='%s%s%sFlat%s%s%s%s%s.fits' % (instr, obsdate, skystr, mdstr, binstr, gnstr, rostr, fltstr)
    return flatname


# -----------------------------------------------------------
# main code

if not iraf.deftask('saltclean'):
   parfile = iraf.osfn("saltred$saltclean.par")
   t = iraf.IrafTaskFactory(taskname="saltclean",value=parfile,function=saltclean, pkgname='saltred')
