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

"""SALTFPPREP is a tool to prepare Fabry-Perot observations for the specific Fabry-Perot reductions using 
T. Williams code.  The code assumes all the files are in the same directory.

Updates:

20100611
    * First wrote the code

20110428
    *Updated the code to work with the output from the salt pipeline
    *Added a number of new features to include in the code
"""

# Ensure python 2.5 compatibility
from __future__ import with_statement

import os
import sys
import numpy as np
#import pyfits

from pyraf import iraf
from pyraf.iraf import pysalt
import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging

from salt2iraf import convertsalt
from salterror import SaltError

from fortranfp import RSSPrep_wrapper
debug=True

def saltfpprep(images, outimages, outpref, slit=False, plottype='postscript', 
               full_reduce=False, clobber=True, logfile='saltfp.log', verbose=True):  
    """Prepares data for Fabry-Perot reductions"""
    if full_reduce:
       rssprep_fortran(images, outpref, slit, logfile, plottype, verbose)
       return 

    with logging(logfile, debug) as log:
       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       #verify that the input and output lists are the same length
       saltio.comparelists(infiles,outfiles,'Input','output')

       for img, oimg in zip(infiles, outfiles):
           message="SALTFPPREP--Preparing image %s" % (img)
           log.message(message, with_header=False)
           try:
               convertsalt(img, oimg, clobber=clobber)
           except SaltError, e:
               log.message('%s' %e)

def rssprep_fortran(inlist, outpref, slit, logfile, plottype, verbose):
    plottype = str(plottype)
    plottype = plottype.strip()
    if plottype == 'xwindow':
        plottype = '/xw'
    else:
        plottype = '/ps'

    pathin = os.path.dirname(inlist)
    basein = os.path.basename(inlist)
    pathlog = os.path.dirname(logfile)
    baselog =  os.path.basename(logfile)

# forcing logfiles to be created in the same directory as the input data
# (we change to this directory once starting the fortran code)

    if len(pathin) > 0:
        logfile = baselog

# start log now that all parameter are set up         

    with logging(logfile, debug) as log:

# set up the variables
        slitstring = " "
        if not slit:
            slitstring = 'n'
        else:
            slitstring = 'y'

# Some basic checks, many tests are done in the FORTRAN code itself
# is the input file specified?
        saltio.filedefined('Input',inlist)


# check to see if input filelist exists, give error if it doesnt.
#        testlistb = indir + infiles
#        testlist = testlistb.strip()
        saltio.fileexists(inlist)

# read the file and check all the files inside exist, give and error if they do not.
        infile= open(inlist, "r")
        line = "dummy start"

        if len(pathin) > 0:
            indir = pathin
        else:
            indir = './'

        while line != "":
            line = infile.readline()
            if line != "":
                if len(pathin) > 0:
                    testfileb = pathin + "/" + line
                else:
                    testfileb = line
                testfile = testfileb.strip()
                saltio.fileexists(testfile)
# check to see if proposed output files exist, if they do remove them.
                if len(pathin) > 0:
                    outfileb = pathin + "/" + outpref + line
                else:
                    outfileb = outpref + line 
                outfile = outfileb.strip()
                if os.path.exists(outfile):
                    os.remove(outfile)
           
        print inlist, indir, slitstring, outpref

        # Get current working directory as the Fortran code changes dir
        startdir = os.getcwd()

        RSSPrep_wrapper.rss_prep(slitstring, indir, inlist, outpref, plottype)

        # go back to starting directory
        os.chdir(startdir)


# -----------------------------------------------------------
# main code

parfile = iraf.osfn("saltfp$saltfpprep.par")
t = iraf.IrafTaskFactory(taskname="saltfpprep",value=parfile,function=saltfpprep,pkgname='saltfp')
