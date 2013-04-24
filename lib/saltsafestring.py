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

# perform string functions on a list
# -----------------------------------------------------------

import os, string
from salterror import SaltError

def listfunc(oldlist,proc):
    """Apply a string function to all entries in a list 

       oldlist: original list

       proc: Process to apply.  Options include upper, lower, lstrip, rstrip, 
             strip, clean (apply all).

    """

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
        message = 'Unknown string function ' + proc
        raise SaltError(message)

    return newlist


def filenumber(filename, x1=9, x2=-5):
    """Extract the file number from a raw filename
 
       *TODO* This will not work if more than 9999 images are taken on a night
    """

    try:
        return int(filename.rstrip()[x1:x2])
    except:
        message = 'Could not extract file number from filename' + filename
        raise SaltError(message)


# construct file name from file number
# -----------------------------------------------------------

def filename(prefix,obsdate,filenum):
   """Construct the file name from a file number"""

   try:
       name = string.zfill(filenum,4)
       name = prefix + obsdate + name + '.fits'
   except:
       message = 'ERROR: SALTSTRING.FILENAME -- could not construct file name from file number' + filenum
       raise SaltError(message)
   return name

def filedate(filename):
   """Extract the date from the filename"""
   return filename[1:9]


def extract(a,b):
    """Find common string in two separate strings"""
    m = min(len(a),len(b))
    for i in range(m):
        if a[i] != b[i]:
            return a[:i]

    return a[:m]

def makeinstrumentstr(instrument):
   """Return a shortened string for the instrument"""
   if instrument=='SALTICAM':
      instr='S'
   elif instrument=='RSS':
      instr='P'
   else:
      instr=''
   return instr

def makeobsdatestr(infiles, x1=5, x2=13):
   """Determine a common obsdate for all infiles"""
   obsdate=os.path.basename(infiles[0])[x1:x2]
   for img in infiles:
        obsdate=extract(obsdate, os.path.basename(img)[x1:x2])
   if len(obsdate)<4: obsdate=''
   return obsdate



def makebinstr(binning):
   """Return a string for the binning"""
   binning=binning.strip().split()
   return '%sx%s' % (binning[0], binning[1])

def makedetmodestr(detmode):
   """Return a shortened string for the obsmode"""
   if detmode=='Normal':
      mdstr='NM'
   elif detmode=='Frame Transfer':
      mdstr='FT'
   elif detmode.upper().count('SLOT'):
      mdstr='SL'
   elif detmode.upper().count('DRIFT'):
      mdstr='DS'
   else:
      mdstr=''
   return mdstr

def makegainstr(gainset):
   """Return a shorted string for the gainsetting"""
   if gainset.upper()=='FAINT':
       gnstr='FA'
   elif gainset.upper()=='BRIGHT':
       gnstr='BR'
   else:
       gnstr=''
   return gnstr

def makereadoutstr(rospeed):
   """Return a shorted string for the read out setting"""
   if rospeed.upper()=='FAST':
       rostr='FA'
   elif rospeed.upper()=='SLOW':
       rostr='SL'
   else:
       rostr=''
   return rostr



def secsplit(sec):
    """Extract x and y ranges from SEC keyword value
       
       --depreciated for getSection in saltsafeio
    """
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
        message = 'Failed to split SEC keyword from ' + sec
        raise SaltError(message)

    return x, y
