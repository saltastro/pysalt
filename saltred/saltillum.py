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
saltillum applies an illumination corection to a data set.   It will median
smooth the data and then divide the data through by the median smoothing

Author                 Version      Date
-----------------------------------------------
SM Crawford (SAAO)    1.0       02 Feb 2012

TODO
-----------------------------------------------


Updates
-----------------------------------------------



"""

from __future__ import with_statement

import time, numpy
from pyraf import iraf

from scipy.ndimage.filters import median_filter

import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging, history

from salterror import SaltError

debug=True



# -----------------------------------------------------------
# core routine

def saltillum(images,outimages,outpref, mbox=11, clobber=False,     \
             logfile='salt.log',verbose=True):


   with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       #verify that the input and output lists are the same length
       saltio.comparelists(infiles,outfiles,'Input','output')

       # open each raw image file
       for infile, outfile in zip(infiles,outfiles):
           struct = saltio.openfits(infile)

           struct = illum_cor(struct, mbox)

           #add any header keywords like history 
           fname, hist=history(level=1, wrap=False)
           saltkey.housekeeping(struct[0],'SILLUM', 'File Illumination corrected', hist)

           #write it out and close it
           saltio.writefits(struct,outfile,clobber=clobber)
           saltio.closefits(struct)

           #output the information
           log.message('Illumination corrected image %s ' % (infile), with_header=False, with_stdout=verbose)


# -----------------------------------------------------------
# Flat field the data

def illum_cor(struct,mbox):
    """Apply an illumination correction to a set of images

       return  struct
    """
    # Determine the number of extensions
    nextend=len(struct)

    #flat field the data
    for i in range(nextend):
       if struct[i].name=='SCI' or len(struct)==1:

            #create the median image
            mdata=median_filter(struct[i].data, size=(mbox, mbox))

            #create the data
            struct[i].data=struct[i].data/mdata

            #Account for variance frames 
            if saltkey.found('VAREXT', struct[i]):
               varext=saltkey.get('VAREXT', struct[i])
               struct[varext].data=struct[varext].data/mdata


    return struct

# -----------------------------------------------------------
# main code
if not iraf.deftask('saltillum'):
   parfile = iraf.osfn("saltred$saltillum.par")
   t = iraf.IrafTaskFactory(taskname="saltillum",value=parfile,function=saltillum, pkgname='saltred')
