################################# LICENSE ##################################
# Copyright (c) 2014, South African Astronomical Observatory (SAAO)        #
# All rights reserved.                                                     #
#                                                                          #
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
from saltfit import *

from salterror import SaltError

debug=True



# -----------------------------------------------------------
# core routine

def saltsurface(images,outimages,outpref, mask=True,  order=3, minlevel=0,
             clobber=False, logfile='salt.log',verbose=True):


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

           
           struct = surface_fit(struct, order=order, mask=mask, minlevel=minlevel)

           #add any header keywords like history 
           fname, hist=history(level=1, wrap=False)
           saltkey.housekeeping(struct[0],'SURFIT', 'File fit by a surface', hist)

           #write it out and close it
           saltio.writefits(struct,outfile,clobber=clobber)
           saltio.closefits(struct)

           #output the information
           log.message('Surface fitted image %s ' % (infile), with_header=False, with_stdout=verbose)



def surface_fit(struct,order=3, mask=True, minlevel=0):
    """Fit a surface to an image 

       return  struct
    """
    # Determine the number of extensions
    nextend=len(struct)

    #flat field the data
    for i in range(nextend):
       if struct[i].name=='SCI' or len(struct)==1:

            #create the median image
            if mask: 
               bpmext = saltkey.get('BPMEXT', struct[i])
               mask_arr = (struct[bpmext].data == 0) * (struct[i].data > minlevel)
            else:
               mask_arr = (struct[i].data < minlevel)

            #create the data
            coef=fit_surface(struct[i].data, mask=mask_arr, xorder=order, yorder=order, xyorder=0)
            y, x = np.indices(struct[i].data.shape)
            sf=surface(coef, xorder=order, yorder=order, xyorder=0)
            print coef
            struct[i].data=sf(x,y)

            #Account for variance frames 
            if mask:
               struct[bpmext].data = struct[bpmext].data*0


    return struct

# -----------------------------------------------------------
# main code
if not iraf.deftask('saltsurface'):
   parfile = iraf.osfn("saltred$saltsurface.par")
   t = iraf.IrafTaskFactory(taskname="saltsurface",value=parfile,function=saltsurface, pkgname='saltred')
