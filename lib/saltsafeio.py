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

All calls to IO functions should go through this module.

**Note** this is a replacement for the `saltio` module.
"""

import os, glob, shutil, re, tempfile
from time import strftime
import smtplib 
from email.mime.text import MIMEText


import pyfits
import numpy as np

from salterror import SaltIOError

def readlist(param):
    """Returns list from epar parameter.

    It accepts the following for *param*
    A filename prefixed by @, for example::

            readlist('@images')

    will read names from the file images.
    A comma separated list of filenames
    """
    
    # Ensure parameter is a string without leading or trailing whitespace
    try:
        param=str(param).strip()
    except:
        raise SaltIOError('Cannot convert input argument to string.')

    if param[0]=='@':
        try:
            f=open(param[1:],'r')
            output=[name.strip() for name in f]
        except:
            raise SaltIOError('Could not read from input file '+param[1:])
    else:
        output=[name.strip() for name in param.split(',')]

    return output

def get_exposure(files, number=1):
    """Given a list of fits files returns exposure data as numpy array.

    By default `get_exposure` will return the first exposure with data.
    If *number* parameter is given it will browse the file list until
    the requested exposure number is found.

    Because of this browing `get_exposure` is only fast for small exposure
    *number*. If the exposure numbers are always big or if you need to
    access a lot of exposures you should use `build_exposure_index` and
    `get_indexed_exposure` instead.
    """

    try:
        f=pyfits.open(files[0])

        # Check if primary HDU contains data
        if f[0].data is None:
            N=len(f)-1
            offset=1
        else:
            N=len(f)
            offset=0

        # Check if requested exposure number is in current file
        if number<=N:
            print 'Exposure found in file ',files[0]
            output=np.asarray(f[number-1+offset].data)
        else:
            output=get_exposure(files[1:],number-N)

        f.close()
    except IOError:
        raise SaltIOError('Cannot open FITS file '+str(files[0]))

    return output

def build_exposure_index(files):
    """Given a list of fits files returns an index of exposures for use
    with `get_indexed_exposure`.
    """

    index=[]

    for name in files:
        try:
            f=pyfits.open(name)

            index+=[(name,number) for number in range(len(f)) if f[number].data is not None]

            f.close()
        except IOError:
            raise SaltIOError('Cannot open FITS file '+str(name))

    return index

def get_indexed_exposure(files,index,number=1):
    """Given a list of fits files and a index generated by
    `build_exposure_index` it returns the requested exposure as a numpy array.

    By default `get_exposure` will return the first exposure with data.
    If *number* parameter is given it will browse the file list until
    the requested exposure number is found.

    The index consumes a lot of memory. Only use `get_indexed_exposure`
    when you need to access a large number of exposures fast, or when the
    exposure number is always big.
    Otherwise you should use `get_exposure` instead.
    """

    try:
        f=pyfits.open(index[number-1][0])

        output=np.asarray(f[index[number-1][1]].data)

        f.close()
    except IOError:
        raise SaltIOError('Cannot open FITS file '+str(index[number-1]))

    return output

def openbinary(file, type):
    """Open binary file."""

    content=[]
    try:
        content = open(file,type)
    except:
        raise SaltIOError('Cannot open ASCII file ' + file)

    return content

def readbinary(content, size, format):
    """read data from a binary file."""

    import struct
    value=''
    try:
        value=content.read(size)
        value=struct.unpack(format,value)[0]
    except:
        raise SaltIOError('Cannot read value')

    return value

def openascii(file,type):
    """open ASCII file"""

    content=[]
    try:
        content = open(file,type)
    except:
        raise SaltIOError('Cannot open ASCII file '+file)

    return content

def closeascii(file):
    """close ASCII file"""

    try:
        file.close()
    except:
        raise SaltIOError('Cannot close ASCII file '+file)

def tmpfile(path):
    """create a temporary file name"""
    try:
        tempfile.tempdir = path
        infile = tempfile.mktemp()
    except Exception, e:
        infile = ''
        raise SaltIOError('Cannot create temporary file name because %s' % e)
    return infile

def openfits(infile, mode='copyonwrite', memmap=False):
    """open FITS file"""
    try:
        struct = pyfits.open(infile, mode=mode, memmap=memmap)
    except Exception, e:
        msg='Cannot open %s as a FITS file because %s'  % (infile, e)
        raise SaltIOError(msg)
        struct = None
    return struct

def openupdatefits(infile):
    """open FITS file for updating"""
    try:
        struct = pyfits.open(infile,mode='update')
    except:
        raise SaltIOError('Cannot open '+infile+' as a FITS file')
        struct = None
    return struct

def updatefits(struct):
    """update existing FITS file"""
    try:
        struct.flush()
    except:
        raise SaltIOError('Cannot update FITS file')

def writefits(struct,outfile, clobber=True):
    """write to FITS file"""
    if (os.path.isfile(outfile) and clobber): 
        delete(outfile)
    try:
        struct.writeto(outfile)
    except Exception,e :
        raise SaltIOError('Cannot write %s because %s' % (outfile, e))

def readimage(struct,hdu):
    """read image from HDU structure"""

    imagedata=[]
    try:
        imagedata = struct[hdu].data
    except:
        raise SaltIOError('Cannot read image data from HDU '+str(hdu))
    return imagedata

def readheader(struct,hdu):
    """read image from HDU structure"""

    headerdata=[]
    try:
        headerdata = struct[hdu].header
    except:
        raise SaltIOError('Cannot read header data from HDU '+str(hdu))
    return headerdata

def writeimage(struct,hdu,imagedata):
    """write image from HDU structure"""
    try:
        struct[hdu].data = imagedata
    except Exception, e:
        raise SaltIOError('Cannot write image data to HDU ' + str(hdu))
    return struct

def readtab(hdu,infile=''):
    """read FITS table HDU"""
    table=''
    
    try:
        table = hdu.data
    except:
        raise SaltIOError('could not extract table from '+infile)

    return table

def fitscolumns(columns):
    """construct FITS table columns"""
    table=''
    try:
        table = pyfits.ColDefs(columns)
    except:
        raise SaltIOError('Cannot define table columns')
    return table

def newfitstable(table,infile=None):
    """write FITS table"""
    struct=''
    try:
        struct = pyfits.new_table(table)
    except Exception, e:
        raise SaltIOError('Cannot create new table because %s' % e)
    return struct

def closefits(struct):
    """close HDU structure"""
    try:
        struct.close()
    except:
        raise SaltIOError('Cannot close HDU structure')

def logname(file):
    """test the log file exists"""
    import string
    newlog = string.join(file.split(),"")

    if (len(newlog) == 0):
        newlog = 'salt.log'

    return newlog

def overwrite(infile,clobber):
    """clobber if file is to be overwritten"""
    if (os.path.isfile(infile) and clobber):
        delete(infile)
    elif (os.path.isfile(infile) and not clobber):
        raise SaltIOError('file '+infile+' exists. use clobber=y')

def fileexists(infile):
    """check that a file exists"""
    if not os.path.isfile(infile):
        raise SaltIOError('File '+infile+' does not exist')

def filedoesnotexist(infile):
    """check that a file does not exist"""
    if os.path.isfile(infile):
        raise SaltIOError('File '+infile+' already exists')

def delete(infile):
    """delete a file"""
    try:
        os.remove(infile)
    except Exception, e:
        raise SaltIOError('Could not delete '+infile+' because '+str(e))

def deletedir(path):
    """delete a directory"""
    try:
        shutil.rmtree(path)
    except:
        raise SaltIOError('Could not delete directory '+path)

def pathexists(path):
    """check that a path exists and name ends with a /"""
    path = path.strip()
    if (path[-1] != '/'): path += '/'
    if (not os.path.exists(path)):
        raise SaltIOError('Path '+path[:-1]+' does not exist')

    return path

def abspath(path):
    """convert relative path to absolute path"""
    try:
        path=pathexists(path)
        curpath = os.getcwd()
        changedir(path)
        path = os.getcwd() + '/'
        changedir(curpath)
    except:
        raise SaltIOError('Could not determine absolute path to '+path)

    return path

def createdir(path):
    """create a directory"""

    path = path.strip()
    if (path[-1] != '/'): path += '/'
    if (not os.path.exists(path)):
        try:
            os.mkdir(path)
        except:
            raise SaltIOError('Could not create directory '+path)

def changedir(path):
    """change working directory"""
    path = path.strip()
    try:
        os.chdir(path)
    except:
        raise SaltIOError('Could not move to directory '+path)

def copy(file1,file2):
    """copy file"""
    try:
        shutil.copy2(file1,file2)
    except Exception, e:
        raise SaltIOError('Could not copy %s to %s due to %s' % (file1, file2, e))

def copydir(file1,file2):
    """copy direcotry recursively"""
    try:
        shutil.copytree(file1,file2)
    except:
        raise SaltIOError('Could not copy ' + file1 + ' to ' + file2)

def move(file1,file2):
    """move file"""
    try:
        shutil.move(file1,file2)
    except Exception,e :
        raise SaltIOError('Could not move %s to %s due to %s' % (file1, file2, e))

def symlink(infile,linkfile,clobber):
    """create symbolic link"""

    # delete file if one of the same name already exists
    if (os.path.exists(linkfile) and not clobber):
        raise SaltIOError('file ' + linkfile + ' exists, use clobber=y')
    if clobber:
        try:
            os.remove(linkfile)
        except:
            pass

    # create symbolic link
    try:
        os.symlink(infile,linkfile)
    except:
        raise SaltIOError('could not create symbolic link from '+infile+' to '+linkfile)

def filedefined(filetype,file):
    """has a file been defined?"""

    file = file.strip()
    if (len(file) == 0 or file.count(' ') > 0):
        raise SaltIOError('filetype '+filetype+'file(s) '+file+' not specified')

def argunpack(argument, value):
    """For arguments that might be a file or list, unpack to make a single list"""
    try:
        argdefined(argument, value)

        if value[0] == '@':
            listexists(argument, value)

        return listparse(argument,value,'','','')
    except:
        raise SaltIOError('Unable to unpack ' + argument)


def argdefined(argument,value):
    """has an argument been defined?"""
    value = value.strip()
    if (len(value) == 0):
        raise SaltIOError(argument + ' argument not defined')

def listexists(filetype,file):
    """does a list file exist, i.e. a parameter that begins with the '@' character"""
    file = file.lstrip('@')
    if not os.path.isfile(file):
        raise SaltIOError(filetype + ' list '+file+' does not exist')

def readimages(filetype, images):
    """Read in and parse a list of input images
    """

    infiles=[]

    # check to see if the file exists
    filedefined(filetype,images)

    listexists(filetype,images)

    infiles=listparse(filetype,images,'',infiles,'')

    filesexist(infiles,'','r')

    return infiles

def listparse(listtype,inlist,pref,altlist,path):
    """create a list from a file or parameter"""
    outlist = []

    # open the file and read in the arguements
    if (len(inlist) > 0 and inlist[0] == '@' and len(pref) == 0):
        line = ' '
        infile = open(inlist.lstrip('@'))
        while line:
            line = infile.readline()
            if (len(line.strip()) > 0) and not line.strip().startswith('#'):
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
        for infile in altlist:
            basefile = os.path.basename(infile)
            outlist.append(pref+basefile)

    # If nothing is found, throw an error
    if (len(outlist) == 0):
        raise SaltIOError(listtype + ' list is empty')

    return outlist

def filesexist(infiles,path,mode):
    """check files in list exist"""

    if (path != ''):
        if (path[len(path)-1] != '/'): path += '/'
    for fileitem in infiles:
        if (mode == 'r'):
            if (not os.path.isfile(path+fileitem)):
                raise SaltIOError('file '+path+fileitem+' does not exist')
        elif (mode == 'w'):
            if (os.path.isfile(path+fileitem)):
                raise SaltIOError('file '+path+fileitem+' already exists')

def comparelists(list1,list2,name1,name2):
    """are two lists the same length?"""
    if (len(list1) != len(list2)):
        raise SaltIOError(name1+' and '+name2+' lists are of unequal length')

def cleanpropcode(pids, propids):
    """Create the list of appropriate proprosal codes"""
    props = []

    # Make a list of all the Propcodes that are observed
    for pid in propids:
        props.extend(pid.upper().split(','))

    # Check to make sure that the requested propocodes
    # are in that night's observation or set the propocodes
    # to all of that nights observatoins
    if (pids[0].upper() != 'ALL'):
        for pid in pids:
            for pid in pid.split(','):
                if (pid.upper().strip() not in set(props)):
                    raise SaltIOError('Propcode ' + pid.upper()+' is not recorded in the observation log ')
    else:
        pids = set(props)

    pids=removebadpids(pids)

    if not pids:
        raise SaltIOError('Propcode list is empty')

    return pids

def removebadpids(pids):
    """Remove propcodes that you do not want --namely junk, bias, test"""
    badnames=('BIAS','COMMON','JUNK','NONE','TEST','UNKNOWN')
    original_pids=set(pids)
    for  pid in original_pids:
        for bn in badnames:
            if pid.upper().count(bn): pids.remove(pid)

    return pids

def removeengineeringpids(pids):
   """Removing propcodes that are associated with engineering and calibration proposals"""
   new_pids=[]
   for pid in pids:
       if not pid.count('ENG_') and not pid.count('CAL_'):
          new_pids.append(pid)
   return new_pids


def readgaindb(gaindb):
    """read gain database file"""
    dbspeed = []
    dbrate  = []
    dbgain  = []
    dbnoise = []
    dbbias  = []
    dbamp   = []

    try:
        gainfile = open(gaindb,'r')
        for line in gainfile:
            if (len(line.strip()) > 0 and line[0] != '#'):
                line = line.rstrip('\r\n')
                line = re.sub("\s+",",",line)
                line.rstrip(',')
                entries = line.split(',')
                dbspeed.append(entries[0])
                dbrate.append(entries[1])
                dbgain.append(entries[2])
                dbnoise.append(entries[3])
                dbbias.append(entries[4])
                dbamp.append(entries[5].strip('amp'))
    except:
        raise SaltIOError('Cannot read gain database file '+gaindb)

    return dbspeed, dbrate, dbgain, dbnoise, dbbias, dbamp

def readxtalkcoeff(xtalkfile):
    """read crosstalk coefficent file"""
    xdict = {}
    try:
        xfile = open(xtalkfile,'r')
        for line in xfile:
            if (len(line.strip()) > 0 and line[0] != '#'):
                line = line.rstrip('\r\n')
                line = line.rstrip()
                line = re.sub("\s+",",",line)
                line.rstrip(',')
                line = line.split(',')
                xdict[int(line[0])]=line[1:]
    except:
        raise SaltIOError('Cannot read crosstalk coefficient file '+xtalkfile)

    return xdict

def readccdgeom(geomfile):
    """read CCD geometry definition file"""

    gap = 0.
    xshift = [0., 0.]
    yshift = [0., 0.]
    rot = [0., 0.]
    try:
        gfile = open(geomfile,'r')
        for line in gfile:
            if (len(line.strip()) > 0 and line[0] != '#'):
                line = line.rstrip('\r\n')
                line = line.rstrip().lstrip()
                line = re.sub("\s+",",",line)
                pars = line.split(',')
                gap = float(pars[1])
                xshift[0] = float(pars[2])
                yshift[0] = float(pars[3])
                rot[0] = float(pars[4])
                if (len(pars) == 8):
                    xshift[1] = float(pars[5])
                    yshift[1] = float(pars[6])
                    rot[1] = float(pars[7])
    except Exception, e :
        raise SaltIOError('Cannot read geometry definition parameters in file %s because %s'%(geomfile, e))

    return gap, xshift, yshift, rot

def checkfornone(inval):
    if inval is None:  return None
    try:
       inval=inval.strip().upper()
    except:
       pass
    if inval in ['NONE','']: return None 
    try:
       if not inval.strip(): return None
    except:
       pass
    return inval

def getSection(section, iraf_format=True):
  """Given an input string for a section in an image, it will
     return the input as a list.  

   section: An input string given a section in the image
            Set to None to return the whole image

   iraf_format: It will invert the x and y values

  """
  #return None if section is None
  if section is None:
     return None

  #remove the brackets
  section=section.replace('[','')
  section=section.replace(']','')

  #loop through the axis
  sect_list=[]
  for s in section.split(','):
      for t in s.split(':'):
          sect_list.append(int(t))

  #flip things around for use with python
  if iraf_format and len(sect_list)==4:
     return [sect_list[2]-1, sect_list[3], sect_list[0]-1, sect_list[1]]

  return sect_list


def ask (msg):
    """Ask for a user response

       returns the response
    """
    resp=''
    try:
        resp=raw_input(msg)
    except Exception, e:
        msg='Could not get response because %s' % e
        raise SaltIOError(msg)

    return resp

def yn_ask(msg):
    """Ask for a [y/n] user response

       returns the response
    """
    resp=ask(msg)

    while not (resp=='y' or resp=='n'):
        pmsg="Please respond with a 'y' or 'n':"
        resp=ask(pmsg)

    if resp=='y':
        resp=True
    else:
        resp=False

    return resp

def email(server,username,password,sender,recipient,bcc, subject,message):
    """send email"""

    # connect to email server
    try:
        smtp = smtplib.SMTP()
        smtp.connect(server)
        smtp.ehlo()
        smtp.starttls()
        smtp.ehlo()
        smtp.login(username,password)
    except Exception, e:
        message = 'Cannot connect to %s because %s' % (server, e)
        raise SaltIOError(message)

    #set up to send to all recipients
    recip = []
    recip.append(recipient)
    for bccobj in bcc.split(','):
        recip.append(bccobj)

    # send emails
    msg = MIMEText(message)
    msg['Subject'] = subject
    msg['From'] = sender
    msg['To'] = recipient
    msg['bcc'] = bcc
    try:
        smtp.sendmail(sender,recip,msg.as_string())
        #print msg
    except Exception, e:
        raise SaltIOError('Failed to send email to %s because %s'% (recipient, e))

    # disconnect from email server
    try:
        smtp.quit()
    except:
        raise SaltIOError('Cannot disconnect from email server '+server)
