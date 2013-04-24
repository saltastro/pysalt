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

import saltprint

# perform string functions on a list
# -----------------------------------------------------------

def listfunc(oldlist,proc,logfile):

    status = 0
    newlist = []
    if (proc == 'upper'):
        for value in oldlist:
            newlist.append(value.upper())
    elif (proc == 'lower'):
        for value in oldlist:
            newlist.append(value.lower())
    elif (proc == 'lstrip'):
        for value in oldlist:
            newlist.append(value.lstrip())
    elif (proc == 'rstrip'):
        for value in oldlist:
            newlist.append(value.rstrip())
    elif (proc == 'clean'):
        for value in oldlist:
            newlist.append(value.lstrip().rstrip().upper())
    else:
        message = 'ERROR: SALTSTRING.LISTFUNC: unknown string function ' + proc
        status = saltprint.err(logfile,message)

    return newlist, status

# extract file number from file name
# -----------------------------------------------------------

def filenumber(filename):

    status = 0
    try:
        fileno = int(filename.rstrip()[9:-5])
    except:
        message = 'ERROR: SALTSTRING.FILENUMBER -- could not extract file number from filename' + filename
        status = saltprint.err(logfile,message)

    return fileno, status

# construct file name from file number
# -----------------------------------------------------------

def filename(prefix,obsdate,filenum):

    status = 0
    try:
        if (filenum < 10):
            name = '000' + str(filenum)
        elif (filenum > 9 and filenum < 100):
            name = '00' + str(filenum)
        elif (filenum > 99 and filenum < 1000):
            name = '0' + str(filenum)
        else:
            name = str(filenum)
        name = prefix + obsdate + name + '.fits'
    except:
        message = 'ERROR: SALTSTRING.FILENAME -- could not construct file name from file number' + filenum
        status = saltprint.err(logfile,message)

    return name, status

# find common string in two separate strings
# -----------------------------------------------------------

def extract(a,b):

    m = min(len(a),len(b))
    for i in range(m):
        if a[i] != b[i]:
            return a[:i]

    return a[:m]

# extract x and y ranges from SEC keyword value
# -----------------------------------------------------------

def secsplit(sec,file,logfile):

    status = 0
    x = [None] * 2
    y = [None] * 2
    try:
        sec = sec.strip('[').strip(']')
        ranges = sec.split(',')
        x[0] = int(ranges[0].split(':')[0])
        x[1] = int(ranges[0].split(':')[1])
        y[0] = int(ranges[1].split(':')[0])
        y[1] = int(ranges[1].split(':')[1])
    except:
        message = 'ERROR: SALTSTRING.SECSPLIT -- failed to split SEC keyword from ' + file
        status = saltprint.err(logfile,status)

    return x, y, status
