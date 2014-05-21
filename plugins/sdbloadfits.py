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

# Author                 Version      Date         Comment
# -----------------------------------------------------------------------
# S M Crawford (SAAO)    0.3          10 Dec 2008

# sdbloadfits adds or updates a fitsdata record in the science database
#
# Limitations
# --Has to be in the directory to work properly

# Ensure python 2.5 compatibility
from __future__ import with_statement


import os, time, glob, string
from pyraf import iraf
import saltsafeio as saltio
import saltsafemysql as saltmysql
from saltsafelog import logging
from salterror import SaltError


debug=True


def sdbloadfits(infile, sdb, logfile, verbose):
    """Add a fits file to the science database
    """
    #determine the base file name
    FileName = findrawfilename(infile)

    #determine the reduced file name
    if FileName!=infile:
            PipelineFileName=infile
    else:
            PipelineFileName=''

    #open the image and extract the header
    infits=saltio.openfits(infile)

    #Extract the primary header
    try:
        ext=0
        ImageHeader = infits[ext].header
    except:
        message='SALTSDBLOADFITS -- ERROR:  Can not access header for %s' % infile
        raise SaltError(message)


    #get the FileData_Id if it already exists
    FileData_Id = checksdbforfits(FileName, sdb, logfile, verbose)

    #create a new entry or update
    if FileData_Id == -1:
        #add to database and get its new id number
        FileData_Id=saltmysql.createnewFileData(sdb,ImageHeader,FileName, PipelineFileName)
    else:
        #update the FileData information
        saltmysql.updateFileData(sdb, ImageHeader, FileData_Id, FileName, PipelineFileName)

    #Update all of the fits header tables
    saltmysql.updateFitsHeaders(sdb, infits, FileData_Id)


    saltio.closefits(infits)

    return 

def checksdbforfits(inname, sdb, log,verbose):
    """
    Check to see if the image is in the database

    returns FileData_Id, status
    """
    FileData_Id=-1

    #check to see if the name is the FileName

    logic="FileName='%s'" % inname
    records=saltmysql.select(sdb,'FileData_Id','FileData',logic)


    if not records:
       FileData_Id=-1
       message='SALTSDBLOADFITS: Adding %s to database' % inname
    elif len(records[0])==1 :
       FileData_Id=records[0][0]
       message='SALTSDBLOADFITS: Updating %s in database' % inname
    elif len(records[0])>1:
       message='SALTLOADSDBFITS -- ERROR: Ambigous record for %s' % inname
       raise SaltError( message)

    if verbose:
       log.message(message, with_header=False)

    return FileData_Id


def findrawfilename(infits):
    i=-1
    try:
        i=infits.index("S")
    except:
        message = infits+' is not a SALT File \n'
    try:
        i=infits.index("P")
    except:
        message = infits+' is not a SALT File \n'

    try:
        i=infits.index("H")
    except:
        message = infits+' is not a SALT File \n'
    try:
        i=infits.index("R")
    except:
        message = infits+' is not a SALT File \n'

    if i>=0:
        inname=infits[i:]
    else:
        inname=''
        message = 'SALTSDBLOADFITS -- ERROR : ' + message
        raise SaltError(message)

    return inname


def findimagenumber (filename):
    """find the number for each image file"""
    #split the file so that name is a string equal to OBSDATE+number
    name=filename.split('/')[-1].split('.')[0]
    return int(name[9:])




