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

"""Module for handling IO in SALT software.

Some calls to FP IO functions should go through this module.

"""

from salterror import SaltIOError
from time import strftime
import pyfits
import numpy as np
import os
import glob

def stringdefined(nameofstring,stringvalue):
    """has a string value been defined?"""

    stringvalue = stringvalue.strip()
    if (len(stringvalue) == 0 or stringvalue.count(' ') > 0):
        raise SaltIOError(nameofstring+' is not specified')

def fplistparse(listtype,inlist,pref,altlist,path):
    """create a list from a file or parameter, allowing user to do this when the list is not in the cwd"""
    outlist = []

    pathin = os.path.dirname(inlist)
    basein = os.path.basename(inlist)
#    print pathin, basein, 'path and base in'
    # open the file and read in the arguments
    if (len(inlist) > 0 and basein[0] == '@' and len(pref) == 0):
        line = 'dummy'
        infile2 = basein.lstrip('@')
        if len(pathin) > 0:
            infile2 = pathin + "/" + infile2
        infile = open(infile2)
        while (len(line.strip()) != 0):
            line = infile.readline()
            if (len(line.strip()) > 0):
                pathtest = os.path.dirname(line)
                basetest = os.path.basename(line)
                if (len(pathtest) == 0 and len(basetest) > 0):
                    pathtest = pathin
                    if (len(pathin) > 0):
                        line = pathtest + "/" + basetest
#                print line.strip(), len(line.strip())
                if (len(line.strip()) > 0):
                    outlist.append(line.rstrip('\r\n'))

    # Include a single entry or a comma separated list of entries.
    elif (len(inlist) > 0 and inlist[0] != '@' and inlist.count('*') == 0 and len(pref)  == 0):
        if (inlist.count(',') == 0):
            outlist.append(inlist)
        else:
            list = inlist.split(',')
            for listitem in list:
                outlist.append(listitem)

    # Include entries with a wildcard
    elif (len(inlist) > 0 and inlist[0] != '@' and inlist.count('*') > 0 and len(pref) == 0):
        globfiles = glob.glob(path+inlist)
        for globitem in globfiles:
            outlist.append(globitem.lstrip(path))
            outlist.sort()

    # Have an alternate or default list to include
    elif (len(pref) > 0):
        for file in altlist:
            filepath = file.split('/')
            outlist.append(pref+filepath[len(filepath)-1])

    # If nothing is found, throw an error
    if (len(outlist) == 0):
        raise SaltIOError(listtype + ' list is empty')

    return outlist


