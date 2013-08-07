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
saltflat is a flatfielding tool for salt images.  It reads
in the data, applies a flatfield, and then produces
flat images

Author                 Version      Date
-----------------------------------------------
SM Crawford (SAAO)    1.0       02 Feb 2009

TODO
-----------------------------------------------


Updates
-----------------------------------------------



"""

from __future__ import with_statement

import time, numpy
from pyraf import iraf

import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging, history

from salterror import SaltError

debug=True



# -----------------------------------------------------------
# core routine

def saltflat(images,outimages,outpref, flatimage,minflat=1, allext=False, clobber=False,     \
             logfile='salt.log',verbose=True):


   with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       #verify that the input and output lists are the same length
       saltio.comparelists(infiles,outfiles,'Input','output')


       # check flatfield image exists
       flaimage= flatimage.strip()
       saltio.fileexists(flatimage)

       flatstruct= saltio.openfits(flatimage)

       # Normalize the flat field image
       # This requires to go through each science extension and divide it by
       # mean of the image.  Things that have to be checked:
       # that data exist, that it is a science extension
       
       #determine the global mean
       fmean=0
       fcount=0
       for fext in flatstruct:
           if fext.shape:
              fmean += fext.data.sum()
              fcount += fext.data.size
       if fcount>0: fmean=fmean/fcount

       for fext in flatstruct:
           if fext.name=='PRIMARY':
                #is it a flat--if not throw a warning
                try:
                    key_ccdtype=saltkey.get('CCDTYPE', fext)
                except:
                    key_ccdtype=None
                if key_ccdtype!='FLAT':
                    message = '\nWARNING -- SALTFLAT: %s does not have CCDTYPE=FLAT' % flatimage
                    log.warning(message,with_stdout=verbose)

                #if there are data, normalize it
                if fext.shape:
                    fext.data=flatnormalize(fext.data, minflat)

           #Noramlize the science extensions
           if fext.name=='SCI':
                if fext.size()>0:
                    if allext is False: fmean=fext.data.mean()
                    fext.data=flatnormalize(fext.data, minflat, fmean)
                    

           #Apply to the variance frames
           if saltkey.found('VAREXT', fext):
               varext=saltkey.get('VAREXT',fext)
               flatstruct[varext].data=flatstruct[varext].data/fmean**2



       # open each raw image file
       for infile, outfile in zip(infiles,outfiles):
           struct = saltio.openfits(infile)

           # flat field correct the image
           outstruct = flat(struct, flatstruct)
           try:
               pass 
           except Exception,e:
               msg='Unable to flatten %s because %s' % (infile, e)
               raise SaltError(msg)

           #add any header keywords like history 
           fname, hist=history(level=1, wrap=False)
           saltkey.housekeeping(struct[0],'SFLAT', 'File flatfield corrected', hist)

           #write it out and close it
           saltio.writefits(outstruct,outfile,clobber=clobber)
           saltio.closefits(struct)

           #output the information
           log.message('Flatfields image %s using %s' % (infile, flatimage), with_header=False, with_stdout=verbose)

       #clost the flatfield image
       saltio.closefits(flatstruct)

# -----------------------------------------------------------
# Normalize the data

def flatnormalize(data, minflat, fmean=1):
    """Calculate the mean of the data and return a normalized data array
    """
    try:
        min_ind=numpy.where(data<minflat)
        data[min_ind]=minflat
        data=data/fmean
    except Exception, e:
        message='ERROR--FLATNORMALIZE:  Data could not be normalized because %s' % e
        raise SaltError(msg)

    return data

    



# -----------------------------------------------------------
# Flat field the data

def flat(struct,fstruct):
    """flat applies a flatfield to salt CCD data

       return  struct
    """
    # Determine the number of extensions
    nextend=len(struct)

    #flat field the data
    for i in range(nextend):
       if struct[i].name=='SCI' or len(struct)==1:
            #Account for variance frames 
            if saltkey.found('VAREXT', struct[i]):
               varext=saltkey.get('VAREXT', struct[i])
               if saltkey.found('VAREXT', fstruct[i]):
                   fvarext=saltkey.get('VAREXT', fstruct[i])
                   fvar=fstruct[fvarext].data
               else:
                   fvar=0
               struct[varext].data=(struct[i].data/fstruct[i].data)*(struct[varext].data/struct[i].data**2+fvar/fstruct[i].data**2)

            #flatten the data
            struct[i].data=struct[i].data/fstruct[i].data

    return struct

# -----------------------------------------------------------
# main code

parfile = iraf.osfn("saltred$saltflat.par")
t = iraf.IrafTaskFactory(taskname="saltflat",value=parfile,function=saltflat, pkgname='saltred')
