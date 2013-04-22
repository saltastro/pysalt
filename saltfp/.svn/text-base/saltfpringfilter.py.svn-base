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

"""
RINGFILTER determines the center coordinates of a ring, bins the ring radially and computes its power spectrum, and allows the user to select a smoothing filter for the ring. It uses T. Williams code.  The code assumes all the files are in the same directory. Also assumes that if there is a config file, it is also in the same directory as the data. Note that this config file is in the original FORTRAN code format so that the user does not have to write another file.

Updates:

20100706
    * First wrote the code
"""

# Ensure python 2.5 compatibility
from __future__ import with_statement

import os
import sys
import numpy as np
#import pyfits

from pyraf import iraf
from pyraf.iraf import pysalt
import saltsafekey
import saltsafeio
import fpsafeio
from saltsafelog import logging
from salterror import SaltIOError

# This reads the FORTRAN config file if it exists

from fortranfp import ringfilter_wrapper
from fortranfp.ringfilter_wrapper import getpfp
debug=True

def saltfpringfilter(axc,ayc,arad,rxc,ryc,filterfreq,filterwidth,itmax,conv, fitwidth,image,logfile,useconfig,configfile,verbose):  

    """ Determines the center coordinates of a ring, bins the ring radially and computes its power spectrum, and allows the user to select a smoothing filter for the ring.  """

# default parameter values are set up in the pyraf .par file. The values used are then changed if a FORTRAN config file exists and the user elects to override the pyraf .par file.

    # Is the input FORTRAN config file specified? 
    # If it is blank, then it will be ignored.        
    if useconfig:
        configfile = configfile.strip()
        if len(configfile) > 0:
            #check exists
            saltsafeio.fileexists(configfile) 

# read updated parameters from the file 
            array=getpfp(configfile,"axc")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                axc=float(array[0])

            array=getpfp(configfile,"ayc")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                ayc=float(array[0])                

            array=getpfp(configfile,"arad")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                arad=float(array[0])
   
            array=getpfp(configfile,"rxc")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                rxc=float(array[0])

            array=getpfp(configfile,"ryc")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                ryc=float(array[0])

            array=getpfp(configfile,"calring_filter_width")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                filterwidth=int(array[0])

            array=getpfp(configfile,"calring_filter_freq")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                filterfreq=int(array[0])

            array=getpfp(configfile,"calring_itmax")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                itmax=int(array[0])

            array=getpfp(configfile,"calring_conv")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                conv=float(array[0])

            array=getpfp(configfile,"calring_fitwidth")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                fitwidth=float(array[0])


# getting paths for filenames

    pathin = os.path.dirname(image)
    basein = os.path.basename(image)
    pathlog = os.path.dirname(logfile)
    baselog =  os.path.basename(logfile)

# forcing logfiles to be created in the same directory as the input data
# (we change to this directory once starting the fortran code)

    if len(pathin) > 0:
        logfile = baselog
    
# start log now that all parameter are set up          
    with logging(logfile, debug) as log:

# Some basic checks, many tests are done in the FORTRAN code itself
# is the input file specified?
        saltsafeio.filedefined('Input',image)

# if the input file is a file, does it exist?
        if basein[0] != '@':
            saltsafeio.fileexists(image)
            infile = image

# if the input file is a list, throw an error
        if basein[0] == '@':
            raise SaltIOError(basein + ' list input instead of a file' )

# optionally update the FORTRAN config file with new values  - not implemented currently

# If all looks OK, run the FORTRAN code
        if len(pathin) > 0:
            dir = pathin
        else:
            dir = './'
        infile = basein
        print dir, infile, 'input directory and input file'

        # Get current working directory as the Fortran code changes dir
        startdir = os.getcwd()
            
        ringfilter_wrapper.ringfilter(dir,axc, ayc,arad, rxc,ryc,filterfreq,filterwidth,itmax,conv,fitwidth,infile)

    # go back to starting directory
    os.chdir(startdir) 

# -----------------------------------------------------------
# main code

parfile = iraf.osfn("saltfp$saltfpringfilter.par")
t = iraf.IrafTaskFactory(taskname="saltfpringfilter",value=parfile,function=saltfpringfilter,pkgname='saltfp')
