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
SALTXTALK--Correct cross talk between the amplifiers for the different CCDS

SALTXTALK corrects data listed by the images argument for amplifier
crosstalk. Each SALTICAM and RSS CCD has two readout
ammplifiers. There is crosstalk between them at the level of ~ 0.1%
which results in faint ghost sources across the image. Ghosts appear
as faint mirror images across amplifier boundaries of bright
sources. Provided images are not saturated or non-linear, crosstalk
can be mostly removed by simple subtraction of a scaled image of one
amplifier from it's neighbour. The scaling factors are supplied as an
ascii table through the xtalkfile argument or in the header keywords


Author                 Version      Date
-----------------------------------------------
Martin Still (SAAO)    0.1          16 Oct 2006
SM Crawford (SAA)      0.2          08 Aug 2011

Updates
-----------------------------------------------
20110808 --Added the new error handling
         --Added ability to read in the values from the header keywords
20130608 --Added variance frames

Todo
-----------------------------------------------

"""
from __future__ import with_statement

import numpy as np
from pyraf import iraf

import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging, history

from salterror import SaltError    

debug=True

# -----------------------------------------------------------
# core routine

def saltxtalk(images,outimages,outpref,xtalkfile=None, usedb=False,
              clobber=True, logfile='salt.log',verbose=True):

   #start logging
   with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       # are input and output lists the same length?
       saltio.comparelists(infiles,outfiles,'Input','output')

       # does crosstalk coefficient data exist
       if usedb:
           xtalkfile = xtalkfile.strip()
           xdict = saltio.readxtalkcoeff(xtalkfile)
       else:
           xdict=None

       for img, oimg in zip(infiles, outfiles):

           #open the fits file
           struct=saltio.openfits(img)

           #find the best xcoeff for the image if using the db
           if usedb:
              obsdate=saltkey.get('DATE-OBS', struct[0])
              obsdate=int('%s%s%s' % (obsdate[0:4],obsdate[5:7], obsdate[8:]))
              xkey=np.array(xdict.keys())
              date=xkey[abs(xkey-obsdate).argmin()]
              xcoeff=xdict[date]
           else:
              xcoeff=[]

           # identify instrument
           instrume,keyprep,keygain,keybias,keyxtalk,keyslot = saltkey.instrumid(struct)

           # has file been prepared already?
           if saltkey.found(keyxtalk, struct[0]):
               message='%s has already been xtalk corrected' % img
               raise SaltError(message)


           #apply the cross-talk correction
           struct = xtalk(struct, xcoeff, log=log, verbose=verbose)

           # housekeeping keywords
           fname, hist=history(level=1, wrap=False, exclude=['images', 'outimages', 'outpref'])
           saltkey.housekeeping(struct[0], 'SXTALK', 'Images have been xtalk corrected', hist)

           # write FITS file
           saltio.writefits(struct,oimg, clobber=clobber)
           saltio.closefits(struct)

# -----------------------------------------------------------
# crosstalk correction

def xtalk(struct,xcoeff,namps=2, log=None,verbose=False):
    """xtalk cross-talk corrects the amplifies.  It takes 

    """
    infile=saltkey.getimagename(struct[0])

    # how many extensions?
    nsciext = saltkey.get('NSCIEXT',struct[0])
    nextend = saltkey.get('NEXTEND',struct[0])
    nccd = saltkey.get('NCCDS',struct[0])

    # how many amplifiers?--this is hard wired
    amplifiers = namps * nccd

    #setup the log
    if log:
       message='%28s %23s' % ('HDU','Correction')
       message += '\n      ----------------------------------------------'
       log.message(message, with_header=False, with_stdout=verbose)

    #Loop over the image extensions and subtract one
    #set from the other--still hardwired at 2
    for i in range(1,nsciext, 2):
      if struct[i].name=='SCI' and struct[i+1].name=='SCI':
        #set up the first amplifier
        dat1=struct[i].data.copy()
        ext1=saltkey.get('EXTVER', struct[i])
        if xcoeff: 
           j=(ext1-1)%amplifiers
           xc1=float(xcoeff[j])
        else:
           xc1=1.0/saltkey.get('XTALK', struct[i])

        #set up the second amplifier
        dat2=struct[i+1].data.copy()
        ext2=saltkey.get('EXTVER', struct[i+1])
        if xcoeff: 
           j=(ext2-1)%amplifiers
           xc2=float(xcoeff[j])
        else:
           xc2=1.0/saltkey.get('XTALK', struct[i+1])

        #subtract one from the other
        struct[i].data=struct[i].data-xc2*dat2[:,::-1]
        struct[i+1].data=struct[i+1].data-xc1*dat1[:,::-1]

        #correct the variance frame
        if saltkey.found('VAREXT', struct[i]):
           vhdu1=saltkey.get('VAREXT', struct[i])
           vhdu2=saltkey.get('VAREXT', struct[i+1])
           try:
               struct[vhdu1].data+=xc2*struct[vhdu2].data
               struct[vhdu2].data+=xc1*struct[vhdu1].data
           except Exception, e:
               msg='Cannot update the variance frame in %s[%i] because %s' % (infile, vhdu1, e)
               raise SaltError(msg)


        #print the message
        if log:
           message = '%25s[%1d]  Amp%1d - Amp%1d * %8.6f' % \
                    (infile, i, ext1, ext2, xc2)
           log.message(message, with_header=False, with_stdout=verbose)
           message = '%25s[%1d]  Amp%1d - Amp%1d * %8.6f' % \
                    (infile, i+1, ext2, ext1, xc1)
           log.message(message, with_header=False, with_stdout=verbose)


    return struct

# -----------------------------------------------------------
# main code
if not iraf.deftask('saltxtalk'):
   parfile = iraf.osfn("saltred$saltxtalk.par")
   t = iraf.IrafTaskFactory(taskname="saltxtalk",value=parfile,function=saltxtalk, pkgname='saltred')
