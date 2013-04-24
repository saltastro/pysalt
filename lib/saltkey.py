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

"""Working with header keys.

.. note::
    Perhaps this module needs to be depricated.
    Most of it's code seems to be dedicated to hiding nice error handling behind ugly status returns.
    Please convert all code that depens on it to use :mod:`salterror` and :mod:`saltsafeio` instead.

.. todo::
    Clean up module to use proper error handling and remove unnecesary code.
"""

from pyraf import iraf
import saltprint, saltio, pyfits, string

def get(keyword,hdu,file,logfile):
    """get keyword value"""

    status = 0
    try:
        value = hdu.header[keyword]
    except:
        message = 'ERROR -- SALTKEY.GET: Cannot read keyword ' + keyword
        message += ' in file ' + file
        status = saltprint.err(logfile,message)
        value = None

    return value, status

def exist(keyword,hdu,file,logfile):
    """does keyword exist"""

    status = 0
    try:
        value = hdu.header[keyword]
    except:
        message = 'ERROR -- SALTKEY.EXIST: ' + keyword + ' does not exist'
        message += ' in file ' + file
        status = saltprint.err(logfile,message)

    return status

def found(keyword,hdu):
    """does keyword exist"""

    status = True
    try:
        value = hdu.header[keyword]
    except:
        status = False

    return status

def match(keyword,value1,hdu,file,logfile):
    """does keyword match prediction"""

    status = 0
    value2, status = get(keyword,hdu,file,logfile)
    if (status == 0 and value1 != value2):
        message  = 'ERROR -- SALTKEY.MATCH -- ' + file
        message += '[' + keyword + '] .not. = ' + value1 + ' : ERROR'
        status = saltprint.err(logfile,message)

    return status

def put(keyword,value,hdu,file,logfile):
    """change existing keyword value"""

    status = 0
    try:
        hdu.header.update(keyword,value)
    except:
        message = 'ERROR -- SALTKEY.PUT: Cannot update keyword ' + keyword
        message += ' in ' + file
        status = saltprint.err(logfile,message)

    return status

def keypar(file,hdu,key,logfile):
    """read keyword woth keypar IRAF tool i.e. without opening the whole file"""

    try:
        iraf.tables.keypar(file+'['+str(hdu)+']',key)
        value = iraf.tables.keypar.value
    except:
        message = 'ERROR -- SALTKEY.KEYPAR: Cannot read keyword ' + key
        message += ' in ' + file + '[' + str(hdu) + ']'
        status = saltprint.err(logfile,message)
        value = None

    return value

def mkheader(file,keyword,value,comment,verbose,logfile):
    """create keyword with mkheader IRAF tool i.e. without opening the whole file"""

    status = 0
    message = 'SALTKEY.MKHEADER: Created keyword ' + keyword + ' in ' + file
    try:
        tmpfile, status = saltio.tmpfile('.',False,logfile)
        tmp, status = saltio.openascii(tmpfile,'w',logfile)
        tmp.write('%-8s= \'%-18s\' / %-s\n' % (keyword,value,comment))
        status = saltio.closeascii(tmp,logfile)
        iraf.noao.artdata.mkheader(file,tmpfile,append='y',verbose='n')
        saltio.delete(tmpfile,False,logfile)
        saltprint.log(logfile,message,verbose)
    except:
        message = 'ERROR -- SALTKEY.MKHEADER: Cannot edit keyword ' + keyword
        message += ' in ' + file
        status = saltprint.err(logfile,message)

    return status

def new(keyword,value,comment,hdu,file,logfile):
    """add new keyword"""

    status = 0
    try:
        hdu.header.update(keyword,value,comment)
    except Exception, e:
        message = 'ERROR -- SALTKEY.NEW: Cannot create keyword %s in %s because %s ' % (keyword, file, e)
        status = saltprint.err(logfile,message)

    return status

def rem(keyword,hdu,file,logfile):
    """delete keyword"""

    status = 0
    try:
        del hdu.header[keyword]
    except:
        message = 'ERROR -- SALTKEY.DEL: Cannot delete keyword '+keyword
        message += ' in '+file
        status = saltprint.err(logfile,message)

    return status

def dateobs(struct,file,logfile):
    """observation date keyword"""

    status = 0
    year = ''
    month = ''
    day = ''
    try:
        date = struct.header['DATE-OBS']
        year = int(date.split('-')[0])
        month = int(date.split('-')[1])
        day = int(date.split('-')[2])
    except:
        message = 'ERROR -- SALTKEY.DATEOBS: DATE-OBS keyword not '
        message += 'recognized in file ' + file
        status = saltprint.err(logfile,message)

    return year, month, day, status

def timeobs(struct,file,logfile):
    """observation time keyword"""

    status = 0
    hour = ''
    minute = ''
    second = ''
    try:
        utcobs = struct.header['UTC-OBS']
        hour = int(utcobs.split(':')[0])
        minute = int(utcobs.split(':')[1])
        second = float(utcobs.split(':')[2])
    except:
        message = 'ERROR -- SALTKEY.TIMEOBS: UTC-OBS keyword not '
        message += 'recognized in file ' + file
        status = saltprint.err(logfile,message)

    return hour, minute, second, status

def ccdbin(struct,file,logfile):
    """CCD binning keyword"""

    status = 0
    xbin = 0
    ybin = 0
    try:
        ccdsum = struct.header['CCDSUM']
        xbin = int(ccdsum.split(' ')[0])
        ybin = int(ccdsum.split(' ')[1])
    except:
        message = 'ERROR -- SALTKEY.CCDBIN: CCDSUM keyword not recognized in file ' + file
        status = saltprint.err(logfile,message)

    return xbin, ybin, status

def instrumid(struct,file,logfile):
    """identify instrument in keywords"""

    instrume = ''
    keyprep = ''
    keygain = ''
    keybias = ''
    keyxtalk = ''
    keyslot = ''
    status = 0
    try:
        instrume = struct[0].header['INSTRUME']
        if (string.join(instrume.split(),"") == 'RSS' or string.join(instrume.split(),"") == 'PFIS'):
            keyprep = 'PPREPARE'
            keygain = 'PGAIN'
            keybias = 'PBIAS'
            keyxtalk = 'PXTALK'
            keyslot = 'PSLOT'
        elif (string.join(instrume.split(),"") == 'SALTICAM'):
            keyprep = 'SPREPARE'
            keygain = 'SGAIN'
            keybias = 'SBIAS'
            keyxtalk = 'SXTALK'
            keyslot = 'SSLOT'
        else:
            message = 'ERROR -- SALTKEY.INSTRUMID: INSTRUME keyword not '
            message += 'recognized in file ' + file
            status = saltprint.err(logfile,message)
    except:
        message = 'ERROR -- SALTKEY.INSTRUMID: INSTRUME keyword not found '
        message += 'in file ' + file
        status = saltprint.err(logfile,message)

    return instrume,keyprep,keygain,keybias,keyxtalk,keyslot,status

def prepare(struct,file,keyprep,logfile):
    """has file been prepared?"""

    status = 0
    try:
        prep = struct[0].header[keyprep]
    except:
        message = 'ERROR -- SALTIO.PREPARE: File ' + file + ' has not been prepared for '
        message += 'SALT IRAF/PyRAF tools'
        status = saltprint.err(logfile,message)

    return status

def clean(struct,file,keyslot,keygain,logfile):
    """has file been cleaned?"""

    status = 0
    try:
        slot = struct[0].header[keyslot]
    except:
        try:
            gain = struct[0].header[keygain]
        except:
            message = 'WARNING -- SALTIO.CLEAN File ' + file + ' is probably in '
            message += 'a raw state'
            status = saltprint.err(logfile,message)

    return status

def history(struct,message,file,logfile):
    """add history keyword"""

    import pyfits

    status = 0
    try:
        struct.header.add_history(message)
    except:
        message = 'ERROR -- SALTKEY.HISTORY: Cannot write HISTORY keyword to ' + file
        status = saltprint.err(logfile,message)

    return status

def copy(new,old,key,logfile):
    """copy keyword"""

    status = 0

    if found(key,old):
        try:
            oldcard=old.header.ascardlist()
            new.header.update(key,old.header[key],oldcard[key].comment)
        except:
            message  = 'ERROR -- SALTKEY.COPY: Cannot COPY KEYWORD ' + key
            status = saltprint.err(logfile,message)


    return status

def housekeeping(hdu, keytask, keycomment, hist, file,logfile):
    """house cleaning keywords"""

    import time
    status = 0
    try:
        status = put('SAL-TLM',time.asctime(time.localtime()),hdu, file,logfile)
        status = new(keytask,time.asctime(time.localtime()),keycomment,hdu,file,logfile)
        status = history(hdu,hist,file,logfile)
    except Exception, e:
        message = 'ERROR--SALTKEY.HOUSECLEANING: Unable to append keywords because %s ' % e
        status = saltprint.err(logfile,message)
    return status

def compare(ahdu, bhdu, keyword, afile, bfile, logfile, verbose):
    """Verify keywords are the same between two different structures
    Compare a keyword between two headers and return a boolean"""

    #get the value in the first header
    avalue,status = get(keyword,ahdu,afile,logfile)

    #get the value in the second header
    if status==0:
        bvalue,status = get(keyword,bhdu,bfile,logfile)

    if status==0 and avalue==bvalue:
        return True

    if status==1:
        message = 'ERROR--SALTKEY.COMPARE: Was not able to compare files %s and %s' % (afile, bfile)
        status = saltprint.err(logfile,message)

    return False
