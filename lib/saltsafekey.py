################################# LICENSE ##################################
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.                                                     #
#                                                                          #
############################################################################

"""Working with header keys.


20091209
   * Changed the behaviour of found to return true, false, or to throw an error

.. note::
    Perhaps this module needs to be depricated.
    Most of it's code seems to be dedicated to hiding nice error handling behind ugly status returns.
    Please convert all code that depens on it to use :mod:`salterror` and :mod:`saltsafeio` instead.

    Do exist and found need to both be present?  Can we clean one up?

.. todo::
    Clean up module to use proper error handling and remove unnecesary code.
"""

import os
from pyraf import iraf
import string
import saltsafeio
from salterror import SaltError, SaltIOError

def getimagename(hdu, base=True):
   """Return the file name of an hdu"""
   try:
       name=hdu._file.name
   except:
       raise SaltIOError('Cannot determine file name')
   if base:  return os.path.basename(name)
   return name


def get(keyword,hdu,file=None):
   """get keyword value"""
       
   try:
        value = hdu.header[keyword]
   except:
       if file is None:
          file=getimagename(hdu, base=False)
       raise SaltIOError('Cannot read keyword '+keyword+' in file '+file)
       value = None

   return value

def exist(keyword,hdu,file=None):
   """does keyword exist"""
   if file is None:
       file=getimagename(hdu, base=False)

   try:
        value = hdu.header[keyword]
   except:
        raise SaltIOError(keyword + ' does not exist in file '+file)

def found(keyword,hdu):
   """does keyword exist.  True if it does, false if a keyerror is thrown, and 
   
   returns bool
   """
   found=False

   #try to access the keyword
   try: 
       hdu.header[keyword]
       found=True
   except KeyError:
       found=False
   except Exception, e:
       raise SaltIOError(e)

   return found
    


def match(keyword,value1,hdu,file=None):
    """does keyword match prediction"""
    if file is None:
       file=getimagename(hdu, base=False)


    value2=get(keyword,hdu,file)
    if value1 != value2:
        raise SaltIOError(file+'['+keyword+'] .not. = '+value1)

def put(keyword,value,hdu,infile=None):
    """change existing keyword value"""


    try:
        hdu.header[keyword] = value
    except:
        if infile is None: infile=getimagename(hdu, base=False)
        raise SaltIOError('Cannot update keyword '+keyword+' in '+infile)

def keypar(file,hdu,key):
    """read keyword woth keypar IRAF tool i.e. without opening the whole file"""

    try:
        iraf.tables.keypar(file+'['+str(hdu)+']',key)
        value = iraf.tables.keypar.value
    except:
        raise SaltIOError('Cannot read keyword '+key+' in '+file+'['+str(hdu)+']')
        value = None

    return value

def mkheader(file,keyword,value,comment):
    """create keyword with mkheader IRAF tool i.e. without opening the whole file"""

    try:
        tmpfile=saltsafeio.tmpfile('.',False)
        tmp=saltsafeio.openascii(tmpfile,'w')
        tmp.write('%-8s= \'%-18s\' / %-s\n' % (keyword,value,comment))
        saltsafeio.closeascii(tmp)
        iraf.noao.artdata.mkheader(file,tmpfile,append='y',verbose='n')
        saltsafeio.delete(tmpfile,False)
    except:
        raise SaltIOError('Cannot edit keyword '+keyword+' in '+file)

def new(keyword,value,comment,hdu,infile=None):
    """add new keyword"""

    try:
       hdu.header[keyword] = (value,comment)
    except Exception, e:
       if infile is None: infile=getimagename(hdu, base=False)
       msg='Cannot create keyword %s in %s because %s ' % (keyword, infile, e)
       raise SaltIOError(msg)

def rem(keyword,hdu,file):
    """delete keyword"""

    try:
        del hdu.header[keyword]
    except:
        raise SaltIOError('Cannot delete keyword '+keyword+' in '+file)

def dateobs(struct,file=None):
    """observation date keyword"""

    year = ''
    month = ''
    day = ''
    try:
        date = struct.header['DATE-OBS']
        year = int(date.split('-')[0])
        month = int(date.split('-')[1])
        day = int(date.split('-')[2])
    except:
       if file is None:
          file=getimagename(hdu, base=False)
       raise SaltIOError('DATE-OBS keyword not recognized in file '+file)

    return year, month, day

def timeobs(struct,file=None):
    """observation time keyword"""

    hour = ''
    minute = ''
    second = ''
    try:
        utcobs = struct.header['UTC-OBS']
        hour = int(utcobs.split(':')[0])
        minute = int(utcobs.split(':')[1])
        second = float(utcobs.split(':')[2])
    except:
        if file is None:
          file=getimagename(hdu, base=False)
        raise SaltIOError('UTC-OBS keyword not recognized in file '+file)

    return hour, minute, second

def ccdbin(struct,file=None):
    """CCD binning keyword"""

    xbin = 0
    ybin = 0
    try:
        ccdsum = struct.header['CCDSUM']
        xbin = int(ccdsum.split(' ')[0])
        ybin = int(ccdsum.split(' ')[1])
    except: 
        if file is None:
            file=getimagename(hdu, base=False)
        raise SaltIOError('CCDSUM keyword not recognized in file '+file)

    return xbin, ybin

def instrumid(struct,file=''):
    """identify instrument in keywords"""
    if not file: struct.filename

    instrume = ''
    keyprep = ''
    keygain = ''
    keybias = ''
    keyxtalk = ''
    keyslot = ''
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
        elif (string.join(instrume.split(),"") == 'HRS'):
            keyprep = 'HPREPARE'
            keygain = 'HGAIN'
            keybias = 'HBIAS'
            keyxtalk = 'HXTALK'
            keyslot = 'HSLOT'
        else:
            raise SaltIOError('INSTRUME keyword not recognized in file '+file)
    except:
        raise SaltIOError('INSTRUME keyword not found in file '+file)

    return instrume,keyprep,keygain,keybias,keyxtalk,keyslot

def prepare(struct,file,keyprep):
    """has file been prepared?"""

    try:
        prep = struct[0].header[keyprep]
    except:
        raise SaltIOError('File '+file+' has not been prepared for SALT IRAF/PyRAF tools')

def clean(struct,file,keyslot,keygain):
    """has file been cleaned?"""

    try:
        slot = struct[0].header[keyslot]
    except:
        try:
            gain = struct[0].header[keygain]
        except:
            raise SaltIOError('File '+file+' is probably in a raw state')

def history(struct,message,infile=None):
    """add history keyword"""
    struct.header.add_history(str(message))
    try:
      pass
    except Exception, e:
        if infile is None: infile=getimagename(struct, base=False)
        raise SaltIOError('Cannot write HISTORY keyword to %s because %s'%(infile, e))

def copy(new,old,key):
    """copy keyword"""

    if found(key,old):
        try:
            oldcard=old.header.cards
            new.header[key] = (oldcard[key].value,oldcard[key].comment)
        except:
            raise SaltIOError('Cannot COPY KEYWORD '+key)

def housekeeping(hdu, keytask, keycomment, hist, infile=None):
    """house cleaning keywords"""

    import time
    if 1:
        put('SAL-TLM',time.asctime(time.localtime()),hdu)
        new(keytask,time.asctime(time.localtime()),keycomment,hdu)
        history(hdu,hist)
    try:
        pass
    except Exception, e:
        msg='Unable to append keywords because %s ' % e
        raise SaltIOError(msg)

def compare(ahdu, bhdu, keyword, afile, bfile):
    """Verify keywords are the same between two different structures
    Compare a keyword between two headers and return a boolean"""

    # get the value in the first header
    avalue=get(keyword,ahdu,afile)

    # get the value in the second header
    bvalue=get(keyword,bhdu,bfile)

    if avalue==bvalue:
        return True
    else:
        return False

def fastmode(mode):
   """Returns True if mode for the detector is one of the high speed readout
       modes
                   
       mode--detectormode
                                 
       returns boolean
   """
   mode=mode.strip().upper()
   if mode in ['VIDEO MODE', 'SLOT MODE', 'SLOTMODE', 'SLOT', 'VIDEO']: return True
   return False
