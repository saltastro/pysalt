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

# saltio is a general library to handle the input
# and output of different types of files

import saltprint
import pyfits, os, sys, string, shutil, glob, re, tempfile, smtplib
from email.mime.text import MIMEText

# -----------------------------------------------------------
# open binary file

def openbinary(file,type,logfile):

    content=[]
    status = 0
    try:
        content = open(file,type)
    except:
        message = 'ERROR: SALTIO.OPENBINARY -- Cannot open ASCII file ' + file
        status = saltprint.err(logfile,message)

    return content,status

#-----------------------------------------------------------
# read data from a binary file

def readbinary(content, size, format,logfile):
    import struct
    value=''
    status=0
    try:
        value=content.read(size)
        value=struct.unpack(format,value)[0]
    except:
        message = 'ERROR: SALTIO.READBINARY -- Cannot read value '
        status = saltprint.err(logfile,message)

    return value, status


# -----------------------------------------------------------
# open ASCII file

def openascii(file,type,logfile):

    content=[]
    status = 0
    try:
        content = open(file,type)
    except:
        message = 'ERROR: SALTIO.OPENASCII -- Cannot open ASCII file ' + file
        status = saltprint.err(logfile,message)

    return content, status

# -----------------------------------------------------------
# close ASCII file

def closeascii(infile,logfile):

    status = 0
    try:
        infile.close()
    except:
        message = 'ERROR: SALTIO.CLOSEASCII -- Cannot close ASCII file ' + infile
        status = saltprint.err(logfile,message)

    return status

# -----------------------------------------------------------
# create a temporary file name

def tmpfile(path,verbose,logfile):

    status = 0
    message = 'SALTIO.TMPFILE -- Created temporary file name '
    try:
        tempfile.tempdir = path
        file = tempfile.mktemp()
        saltprint.log(logfile,message + file,verbose)
    except:
        file = ''
        message = 'ERROR -- SALTPRINT.TMPFILE: Cannot create temporary file name'
        status = saltprint.err(logfile,message)
    return file, status

# -----------------------------------------------------------
# open FITS file

def openfits(file,logfile):

    status = 0
    try:
        struct = pyfits.open(file)
    except:
        message = 'ERROR -- SALTIO.OPENFITS: cannot open ' + file
        message += ' as a FITS file'
        status = saltprint.err(logfile,message)
        struct = None
    return struct, status

# -----------------------------------------------------------
# open FITS file for updating

def openupdatefits(file,logfile):

    status = 0
    try:
        struct = pyfits.open(file,mode='update')
    except:
        message = 'ERROR -- SALTIO.OPENUPDATEFITS: cannot open ' + file
        message += ' as a FITS file'
        struct = None
        status = saltprint.err(logfile,message)
    return struct, status

# -----------------------------------------------------------
# update existing FITS file

def updatefits(struct,logfile):

    status = 0
    try:
        struct.flush()
    except:
        message = 'ERROR -- SALTIO.UPDATEFITS: cannot update FITS file'
        status = saltprint.err(logfile,message)
    return status

# -----------------------------------------------------------
# write to FITS file

def writefits(struct,file,logfile):

    status = 0
    try:
        struct.writeto(file)
    except Exception, e:
        message = 'ERROR -- SALTIO.WRITE: cannot write %s because %s' % (file, e)
        status = saltprint.err(logfile,message)
    return status


# -----------------------------------------------------------
# read image from HDU structure

def readimage(struct,hdu,logfile):

    imagedata=[]
    status = 0
    try:
        imagedata = struct[hdu].data
    except:
        message = 'ERROR -- SALTIO.READIMAGE: cannot read image data from HDU ' + str(hdu)
        status = saltprint.err(logfile,message)
    return imagedata, status

# -----------------------------------------------------------
# write image from HDU structure

def writeimage(struct,hdu,imagedata,logfile):

    status = 0
    try:
        struct[hdu].data = imagedata
    except Exception, e:
        message = 'ERROR -- SALTIO.WRITEIMAGE: Cannot write image data to HDU ' + str(hdu)
        message += ' because %s ' % e
        status = saltprint.err(logfile,message)
    return struct, status

# -----------------------------------------------------------
# read FITS table HDU

def readtab(hdu,file,logfile):

    table=''
    status = 0
    message = 'ERROR -- SALTIO.READTAB: could not extract table from ' + file
    try:
        table = hdu.data
    except:
        status = saltprint.err(logfile,message)
    return table, status

# -----------------------------------------------------------
# construct FITS table columns

def fitscolumns(columns,file,logfile):

    table=''
    status = 0
    try:
        table = pyfits.ColDefs(columns)
    except:
        message = 'ERROR -- SALTIO.FITSCOLUMNS: cannot define table columns in ' + file
        status = saltprint.err(logfile,message)
    return table, status

# -----------------------------------------------------------
# write FITS table

def newfitstable(table,file,logfile):

    struct=''
    status = 0
    try:
        struct = pyfits.new_table(table)
    except:
        message = 'ERROR -- SALTIO.NEWFITSTABLE: cannot create new table in ' + file
        status = saltprint.err(logfile,message)
    return struct, status

# -----------------------------------------------------------
# close HDU structure

def closefits(struct,logfile):

    status = 0
    try:
        struct.close()
    except:
        message = 'ERROR -- SALTIO.CLOSE: Cannot close HDU structure'
        status = saltprint.err(logfile,message)
    return status

# -----------------------------------------------------------
# test the log file exists

def logname(file):

    newlog = string.join(file.split(),"")
    if (len(newlog) == 0): newlog = 'salt.log'
    if (newlog != file):
        message = 'WARNING -- SALTIO.LOGNAME: logfile renamed to ' + newlog + '\n\n'
        log(newlog,message,chatter)
    return newlog

# -----------------------------------------------------------
# clobber if file is to be overwritten

def overwrite(file,clobber,logfile):

    status = 0
    if (os.path.isfile(file) and clobber):
        status = delete(file,False,logfile)
    elif (os.path.isfile(file) and not clobber):
        message = 'ERROR: SALTIO.OVERWRITE -- file ' + file + ' exists. use clobber=y'
        status = saltprint.err(logfile,message)

    return status

# -----------------------------------------------------------
# check that a file exists

def fileexists(file,logfile):

    status = 0
    if not os.path.isfile(file):
        message = 'ERROR -- SALTIO.FILEEXISTS: File ' + file
        message += ' does not exist'
        status = saltprint.err(logfile,message)
    return status

# -----------------------------------------------------------
# check that a file does not exist

def filedoesnotexist(file,verbose,logfile):

    status = 0
    if os.path.isfile(file):
        if (verbose == 'no'):
            status = 1
        else:
            message = 'ERROR -- SALTIO.FILEDOESNOTEXIST: File ' + file
            message += ' already exists'
            status = saltprint.err(logfile,message)
    return status

# -----------------------------------------------------------
# delete a file

def delete(file,verbose,logfile):

    status = 0
    message = 'SALTIO.DELETE -- deleted file ' + file
    try:
        os.remove(file)
        saltprint.log(logfile,message,verbose)
    except Exception, e:
        message = 'ERROR -- SALTIO.DELETE: Could not delete %s because %s' % (file ,e)
        status = saltprint.err(logfile,message)
    return status

# -----------------------------------------------------------
# delete a directory

def deletedir(path,logfile):

    status = 0
    try:
        os.rmdir(path)
    except:
        message = 'ERROR -- SALTIO.DELETEDIR: Could not delete directory' + file
        status = saltprint.err(logfile,message)
    return status

# -----------------------------------------------------------
# check that a path exists and name ends with a "/"

def pathexists(path,logfile):

    status = 0
    path = path.strip()
    if (path[-1] != '/'): path += '/'
    if (not os.path.exists(path)):
        message = 'ERROR -- SALTIO.PATHEXISTS: Path ' + path[:-1] + ' does not exist'
        status = saltprint.err(logfile,message)

    return path, status

# -----------------------------------------------------------
# convert relative path to absolute path

def abspath(path,logfile):

    status = 0
    path, status = pathexists(path,logfile)
    if (status == 0):
        curpath = os.getcwd()
        status = changedir(path,False,logfile)
    if (status == 0):
        path = os.getcwd() + '/'
    else:
        message = 'ERROR: SALTIO.ABSPATH -- cannot determine absolute path of ' + path
        status = saltprint.err(logfile,message)
    if status == 0:
        status = changedir(curpath,False,logfile)

    return path, status

# -----------------------------------------------------------
# create a directory

def createdir(path,verbose,logfile):

    status = 0
    path = path.strip()
    message = 'SALTIO.CREATEDIR -- Created directory ' + path
    if (path[-1] != '/'): path += '/'
    if (not os.path.exists(path)):
        try:
            os.mkdir(path)
            saltprint.log(logfile,message,verbose)
        except:
            message = 'ERROR -- SALTIO.CREATEDIR: Could not create directory ' + path
            status = saltprint.err(logfile,message)
    else:
        message = 'SALTIO.CREATEDIR -- ' + path + ' directory exists'
        saltprint.log(logfile,message,verbose)

    return status

# -----------------------------------------------------------
# change working directory

def changedir(path,verbose,logfile):

    status = 0
    path = path.strip()
    message = 'SALTIO.CHANGEDIR -- Moved to directory ' + path
    try:
        os.chdir(path)
        saltprint.log(logfile,message,verbose)
    except:
        message = 'ERROR -- SALTIO.CREATEDIR: Could not move to directory ' + path
        status = saltprint.err(logfile,message)

    return status

# -----------------------------------------------------------
# copy file

def copy(file1,file2,verbose,logfile):

    status = 0
    message = 'SALTIO.COPY -- copied ' + file1 + ' to ' + file2
    try:
        shutil.copy2(file1,file2)
        saltprint.log(logfile,message,verbose)
    except:
        message = 'ERROR -- SALTIO.COPY: could not copy ' + file1 + ' to ' + file2
        status = saltprint.err(logfile,message)

    return status

# -----------------------------------------------------------
# copy direcotry recursively

def copydir(file1,file2,verbose,logfile):

    status = 0
    message = 'SALTIO.COPYDIR -- copied ' + file1 + ' to ' + file2
    try:
        shutil.copytree(file1,file2)
    except Exception, e:
        message = 'ERROR -- SALTIO.COPYDIR: could not copy %s to %s because %s' % (file1,file2,e)
        status = saltprint.err(logfile,message)

    return status

# -----------------------------------------------------------
# move file

def move(file1,file2,verbose,logfile):

    status = 0
    message = 'SALTIO.MOVE -- moved ' + file1 + ' to ' + file2
    try:
        shutil.move(file1,file2)
        saltprint.log(logfile,message,verbose)
    except:
        message = 'ERROR -- SALTIO.MOVE: Could not move ' + file1 + ' to ' + file2
        status = saltprint.err(logfile,message)

    return status

# -----------------------------------------------------------
# create symbolic link

def symlink(infile,linkfile,clobber,verbose,logfile):

# delete file if one of the same name already exists

    status = 0
    message = 'SALTIO.SYMLINK -- created symbolic link from ' + infile + ' to ' + linkfile
    if (os.path.exists(linkfile) and not clobber):
        message = 'ERROR: SALTIO.SYMLINK -- file ' + linkfile + ' exists, use clobber=y'
        status = saltprint.err(logfile,message)
    if (status == 0 and clobber):
        try:
            os.remove(linkfile)
        except:
            status = 0

# create symbolic link

    if (status == 0):
        try:
            os.symlink(infile,linkfile)
        except:
            message  = 'ERROR: SALTIO.SYMLINK -- could not create symbolic link from '
            message += infile + ' to ' + linkfile
            status = saltprint.err(logfile,message)
    if (status == 0): saltprint.log(logfile,message,verbose)

    return status

# -----------------------------------------------------------
# has a file been defined?

def filedefined(filetype,file,logfile):

    status = 0
    file = file.strip()
    if (len(file) == 0 or file.count(' ') > 0):
        message = 'ERROR -- SALTIO.FILEDEFINED: ' + filetype + ' file(s) not specified'
        status = saltprint.err(logfile,message)
    return status

#-----------------------------------------------------------
#For arguments that might be a file or list, unpack to make
#a single list

def argunpack(argument, value, logfile):
    import saltstring
    arg=''

    #check that it is defined
    status= argdefined(argument, value, logfile)

    #check to see if it is a list in a file, if it exists
    if (status == 0 and value[0] == '@'):
        status = listexists(argument, value,logfile)

    #read into the list 
    if (status == 0 ) :
       arg,status=listparse(argument,value,'','','',logfile)
       #arg,status=saltstring.listfunc(arg,'clean',logfile)

    if status==1:
        message = 'ERROR -- SALTIO.ARGUNPACK: Unable to unpack ' + argument + ' '
        status = saltprint.err(logfile,message)
    return arg,status

# -----------------------------------------------------------
# has an argument been defined?

def argdefined(argument,value,logfile):

    status = 0
    value = value.strip()
    if (len(value) == 0): # or value.count(' ') > 0):
        message = 'ERROR -- SALTIO.ARGDEFINED: ' + argument + ' argument not defined'
        status = saltprint.err(logfile,message)
    return status

# -----------------------------------------------------------
# does a list file exist, i.e. a parameter that begins with the '@' character

def listexists(filetype,file,logfile):

    status = 0
    file = file.lstrip('@')
    if not os.path.isfile(file):
        message = 'ERROR -- SALTIO.LISTEXIST: ' + filetype + ' list '+ file
        message += ' does not exist'
        status = saltprint.err(logfile,message)
    return status

# ----------------------------------------------------------
# Read in input image files

def readimages(filetype, images, logfile):
    """Read in and parse a list of input images
    """

    infiles=[]
    #check to see if the file exists
    status = filedefined(filetype,images,logfile)

    if status==0:
        status = listexists(filetype,images,logfile)

    if (status == 0):
        infiles, status = listparse(filetype,images,'',infiles,'',logfile)

    if (status == 0): status = filesexist(infiles,'','r',logfile)

    return infiles, status

# -----------------------------------------------------------
# create a list from a file or parameter

def listparse(listtype,inlist,pref,altlist,path,logfile):

    status = 0
    outlist = []
    #open the file and read in the arguements
    if (len(inlist) > 0 and inlist[0] == '@' and len(pref) == 0):
        line = ' '
        infile = open(inlist.lstrip('@'))
        while line:
            line = infile.readline()
            if (len(line.strip()) > 0):
                outlist.append(line.rstrip('\r\n'))
    #Include a single entry or a comma separated list of entries.
    elif (len(inlist) > 0 and inlist[0] != '@' and inlist.count('*') == 0 and len(pref)  == 0):
        if (inlist.count(',') == 0):
            outlist.append(inlist)
        else:
            list = inlist.split(',')
            for listitem in list:
                outlist.append(listitem)

    #Include entries with a wildcard
    elif (len(inlist) > 0 and inlist[0] != '@' and inlist.count('*') > 0 and len(pref) == 0):
        globfiles = glob.glob(path+inlist)
        for globitem in globfiles:
            outlist.append(globitem.lstrip(path))
            outlist.sort()

    #have an alternate or default list to include
    elif (len(pref) > 0):
        for file in altlist:
            filepath = file.split('/')
            outlist.append(pref+filepath[len(filepath)-1])
    #if nothing is found, throw an error
    if (len(outlist) == 0):
        message = 'ERROR -- SALTIO.LISTPARSE: ' + listtype + ' list is empty'
        status = saltprint.err(logfile,message)
    #outlist.sort()
    return outlist, status

# -----------------------------------------------------------
# check files in list exist

def filesexist(infiles,path,mode,logfile):

    status = 0
    if (path != ''):
        if (path[len(path)-1] != '/'): path += '/'
    for fileitem in infiles:
        if (mode == 'r'):
            if (not os.path.isfile(path+fileitem)):
                message = 'ERROR -- SALTIO.FILESEXIST file ' + path + fileitem
                message += ' does not exist'
                status = saltprint.err(logfile,message)
        elif (mode == 'w'):
            if (os.path.isfile(path+fileitem)):
                message = 'ERROR -- SALTIO.FILESEXIST file ' + path + fileitem
                message += ' already exists'
                status = saltprint.err(logfile,message)
    return status

# -----------------------------------------------------------
# are two lists the same length?

def comparelists(list1,list2,name1,name2,logfile):

    status = 0
    if (len(list1) != len(list2)):
        message = 'ERROR -- SALTIO.COMPARELISTS: ' + name1 + ' and ' + name2
        message += ' lists are of unequal length'
        status = saltprint.err(logfile,message)
    return status

#-----------------------------------------------------------
#Create the list of appropriate proprosal codes

def cleanpropcode(pids, propids, logfile):
    status=0
    props = []

    #Make a list of all the Propcodes that are observed
    for pid in propids:
        props.extend(pid.upper().split(','))

    #Check to make sure that the requested propocodes
    # are in that night's observation or set the propocodes
    # to all of that nights observatoins
    if (pids[0].upper() != 'ALL'):
        for pid in pids:
            for pid in pid.split(','):
                if (pid.upper().strip() not in set(props)):
                    message  = 'ERROR: CLEANPROPCODE -- Propcode ' + pid.upper()
                    message += ' is not recorded in the observation log '
                    status = saltprint.err(logfile,message)
    else:
        pids = set(props)

    pids=removebadpids(pids)

    if not pids:
        message  = 'ERROR: CLEANPROPCODE -- Propcode list is empty'
        status = saltprint.err(logfile,message)
    return pids, status

# -----------------------------------------------------------
# Remove propcodes that you do not want --namely junk, bias, test

def removebadpids(pids):
    badnames=('BIAS','COMMON','JUNK','NONE','TEST','UNKNOWN')
    original_pids=set(pids)
    for  pid in original_pids:
        for bn in badnames:
            if pid.upper().count(bn): pids.remove(pid)

    return pids

# -----------------------------------------------------------
# read gain database file

def readgaindb(gaindb,logfile):

    dbspeed = []
    dbrate  = []
    dbgain  = []
    dbnoise = []
    dbbias  = []
    dbamp   = []
    status  = 0
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
        message = 'Cannot read gain database file ' + gaindb
        status = saltprint.err(logfile,message)

    return dbspeed, dbrate, dbgain, dbnoise, dbbias, dbamp, status

# -----------------------------------------------------------
# read crosstalk coefficent file

def readxtalkcoeff(xtalkfile,logfile,status):

    xcoeff = []
    status=0
    try:
        xfile = open(xtalkfile,'r')
        for line in xfile:
            if (len(line.strip()) > 0 and line[0] != '#'):
                line = line.rstrip('\r\n')
                line = line.rstrip()
                line = re.sub("\s+",",",line)
                line.rstrip(',')
                xcoeff = line.split(',')
    except:
        message = 'ERROR -- SALTIO.READXTALKCOEFF: Cannot read crosstalk coefficient '
        message += 'file ' + xtalkfile
        status = saltprint.err(logfile,message)

    return xcoeff, status

# -----------------------------------------------------------
# read CCD geometry definition file

def readccdgeom(geomfile,logfile,status):

    status = 0
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
    except:
        message = 'ERROR -- SALTIO.READCCDGEOM: Cannot read geometry definition parameters '
        message += 'in file ' + geomfile
        status = saltprint.err(logfile,message)

    return gap, xshift, yshift, rot, status

# -----------------------------------------------------------
# ask for a response from the user

def ask (msg, logfile):
    """Ask for a user response

       returns the response
    """
    resp=''
    status=0
    try:
        resp=raw_input(msg)
    except Exception, e:
        msg='ERROR--SALTIO.ASK: Could not get response because %s' % e
        status = saltprint.err(logfile,message)

    return resp,status

# -----------------------------------------------------------
# ask for a [y/n] response from the user

def yn_ask (msg, logfile):
    """Ask for a [y/n] user response

       returns the response
    """
    resp, status =ask(msg, logfile)

    if (status==0):
        while not (resp=='y' or resp=='n') and status==0:
            pmsg="Please respond with a 'y' or 'n':"
            resp, status = ask (pmsg, logfile)

    if (status==0):
        if resp=='y':
            resp=True
        else:
            resp=False

    return resp, status





# -----------------------------------------------------------
# send email

def email(server,username,password,sender,recipient,subject,message,logfile):

    status = 0

# connect to email server

    try:
        smtp = smtplib.SMTP()
        smtp.connect(server)
    except:
        message = 'ERROR: SALTIO.EMAIL -- cannot connect to email server ' + server
        status = saltprint.err(logfile,message)

# login to email server

    if (status == 0):
        try:
            smtp.login(username,password)
        except:
            message = 'ERROR: SALTEMAIL -- cannot login to email server ' + server + ' as user ' + username
            status = saltprint.err(logfile,message)

# send emails

    if (status == 0):
        msg = MIMEText(message)
        msg['Subject'] = subject
        msg['From'] = sender
        msg['To'] = recipient
        try:
            smtp.sendmail(sender,recipient,msg.as_string())
        except:
            message = 'ERROR: SALTEMAIL -- failed to send email to ' + recipient
            status = saltprint.err(logfile,message)

# disconnect from email server

    if (status == 0):
        try:
            smtp.quit()
        except:
            message = '\nERROR: SALTEMAIL -- cannot disconnect from email server ' + server
            status = saltprint.err(logfile,message)

    return status
