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

from fortranfp import calprofile_wrapper
from fortranfp.calprofile_wrapper import getpfp
debug=True

def saltfpcalprofile(axc,ayc,arad,rxc,ryc,filter, filterfreq,filterwidth,plottype,itmax,conv, fitwidth,rlo,rhi,rfixed,cenfixed,images,outfile,comment,cala, calb, calc, cald, calf,calprofilelogfile,logfile,useconfig,configfile,verbose):  

    """Fits a Voigt profile to a ring in one or more Fabry-Perot calibration ring images"""

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

            array=getpfp(configfile,"calring_filter")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                filter=str(array[0])


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

            array=getpfp(configfile,"calring_rlo")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                rlo=float(array[0])

            array=getpfp(configfile,"calring_rhi")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                rhi=float(array[0])

            array=getpfp(configfile,"calring_fixed")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                cenfixed=str(array[0])
      
            
            array=getpfp(configfile,"calibration_a")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                cala=float(array[0])


            array=getpfp(configfile,"calibration_b")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                calb=float(array[0])

            array=getpfp(configfile,"calibration_c")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                calc=float(array[0])

            array=getpfp(configfile,"calibration_d")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                cald=float(array[0])

            array=getpfp(configfile,"calibration_f")
            s=len(array)
            flag = array[s-1]
            if flag == 1: 
                calf=float(array[0])


    cenfixed=str(cenfixed)
    if cenfixed[0] == 'y' or cenfixed[0] == 'Y':
        cenfixed = 1
    else:
        cenfixed = 0
    cenfixed = bool(cenfixed)
        
    filter = str(filter)
    if filter[0] == 'y' or filter[0] == 'Y':
        filter = 1
    else:
        filter = 0
    filter = bool(filter)

    plottype = str(plottype)
    plottype = plottype.strip()
    if plottype == 'xwindow':
        plottype = '/xw'
    else:
        plottype = '/ps'
        
# getting paths for filenames

    pathin = os.path.dirname(images)
    basein = os.path.basename(images)
    pathout = os.path.dirname(outfile)
    baseout = os.path.basename(outfile)
    pathlog = os.path.dirname(logfile)
    baselog =  os.path.basename(logfile)
    pathcalprofilelog = os.path.dirname(calprofilelogfile)
    basecalprofilelog =  os.path.basename(calprofilelogfile)

# forcing logfiles to be created in the same directory as the input data
# (we change to this directory once starting the fortran code)

    if len(pathin) > 0:
        logfile = baselog
        calprofilelogfile = basecalprofilelog
    
# start log now that all parameter are set up          
    with logging(logfile, debug) as log:

# Some basic checks, many tests are done in the FORTRAN code itself
# is the input file specified?
        saltsafeio.filedefined('Input',images)

# if the input file is a file, does it exist?
        if basein[0] != '@':
            saltsafeio.fileexists(images)
            infile = images

# if the input file is a list, does it exist?
        if basein[0] == '@':
            infile2 = basein.lstrip('@')
            if (len(pathin) > 0):
                infile = pathin + '/' + infile2
            else:
                infile = infile2
            saltsafeio.listexists('Input',infile)

# parse list of input files
            infiles=fpsafeio.fplistparse('Raw image',images,'','','')
# check input files exist
            saltsafeio.filesexist(infiles,'','r')

# is the proposed output file specified?

        saltsafeio.filedefined('Output',outfile)

# check whether the proposed output file exists

        path = os.path.dirname(images)

        if len(path) > 0:
            dir = path + '/'
        else:
            dir = './'

#        print dir, ' directory with data'
        outfile = outfile.strip()
        if os.path.isfile(outfile):
            print 'output file exists, appending'
#            saltsafeio.delete(outfile)


# optionally update the FORTRAN config file with new values  - not implemented currently

# If all looks OK, run the FORTRAN code

        infile = basein
        if len(pathout) > 0 :
            if (pathin == pathout):
                outfile = baseout
        # absolute path then leave alone
            elif pathout[0] == '/':
                outfile = outfile
            else:
   # as will move to directory of input data, need to strip this off if its in the path
                a = pathout.lstrip(pathin)
                b = a.lstrip('/')
                outfile = b + '/' + baseout
            
        
    
        print dir, infile, outfile, 'input directory, input and output files'
        if basein[0] == '@':
            print infiles,'infiles'

        # Get current working directory as the Fortran code changes dir
        startdir = os.getcwd()
            
        calprofile_wrapper.calprofile(dir,calprofilelogfile,outfile,comment, axc, ayc,arad, rxc,ryc,filter,filterfreq,filterwidth,plottype,itmax,conv,fitwidth,rlo,rhi,rfixed, cenfixed,infile, cala, calb, calc, cald, calf)

        # go back to starting directory
        os.chdir(startdir)   

# -----------------------------------------------------------
# main code

parfile = iraf.osfn("saltfp$saltfpcalprofile.par")
t = iraf.IrafTaskFactory(taskname="saltfpcalprofile",value=parfile,function=saltfpcalprofile,pkgname='saltfp')
