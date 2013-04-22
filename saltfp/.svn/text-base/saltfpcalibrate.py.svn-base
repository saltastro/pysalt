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

"""SALTFPCALPROFILE fits a Voigt profile to one or more calibration rings,
displays the fit, and records the fit parameters for later use. It uses T. Williams code.  The code assumes all the files are in the same directory. Also assumes that if there is a config file, it is also in the same directory as the data. Note that this config file is in the original FORTRAN code format so that the user does not have to write another file.

Updates:

20100623
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

# This reads the FORTRAN config file if it exists

from fortranfp import calibrate_wrapper
debug=True

def saltfpcalibrate(plottype,infile,outfile,calibratelogfile,logfile,verbose):  

    """Determines the wavelength calibration of the Fabry-Perot    
   interactively, using the ring measurements from the CALRING program."""

    plottype = str(plottype)
    plottype = plottype.strip()
    if plottype == 'xwindow':
        plottype = '/xw'
    else:
        plottype = '/ps'
            
# start log now that all parameter are set up          
    with logging(logfile, debug) as log:

# Some basic checks, many tests are done in the FORTRAN code itself
# is the input file specified?
        saltsafeio.filedefined('Input',infile)

# The input file in this case is a file, does it exist?
        if infile[0] != '@':
            saltsafeio.fileexists(infile)

# is the proposed output file specified?

        saltsafeio.filedefined('Output',outfile)

# check whether the proposed output file exists

        path = os.path.dirname(infile)

        if len(path) > 0:
            dir = path + '/'
        else:
            dir = './'

        print dir, ' directory with data'
        outfile = outfile.strip()
        if os.path.isfile(outfile):
            print 'output file exists, appending'
#            saltsafeio.delete(outfile)

# check whether the calibrate logfile is defined

        saltsafeio.filedefined('Log',calibratelogfile)

        # Get current working directory as the Fortran code changes dir
        startdir = os.getcwd()
    
# If all looks OK, run the FORTRAN code
            
        calibrate_wrapper.calibrate(plottype,infile, outfile, calibratelogfile)

        # go back to starting directory
        os.chdir(startdir)   

# -----------------------------------------------------------
# main code

parfile = iraf.osfn("saltfp$saltfpcalibrate.par")
t = iraf.IrafTaskFactory(taskname="saltfpcalibrate",value=parfile,function=saltfpcalibrate,pkgname='saltfp')
