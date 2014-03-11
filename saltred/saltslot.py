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

"""SALTSLOT--Fast data reduction for slot mode data

Author                 Version      Date
-----------------------------------------------
Martin Still (SAAO)    0.1          13 Nov 2006
SM Crawford  (SAAO)    0.1          08 Aug 2011

Updates
-----------------------------------------------
20110808   Re-wrote to use the new error handling

"""

from pyraf import iraf
from pyraf.iraf import pysalt
import os, time

import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging, history

from salterror import SaltError

debug=True

# -----------------------------------------------------------
# core routine

def saltslot(images,outimages,outpref,gaindb='',xtalkfile='',usedb=False, clobber=False,logfile='salt.log',verbose=True):


   #start logging
   with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       # are input and output lists the same length?
       saltio.comparelists(infiles,outfiles,'Input','output')

       # does crosstalk coefficient data exist
       if usedb:
           dblist= saltio.readgaindb(gaindb)
           xtalkfile = xtalkfile.strip()
           xdict = saltio.readxtalkcoeff(xtalkfile)
       else:
           dblist=[]
           xdict=None

       for img, oimg in zip(infiles, outfiles):
           #open the fits file
           struct=saltio.openfits(img)

           # identify instrument
           instrume,keyprep,keygain,keybias,keyxtalk,keyslot = saltkey.instrumid(struct)

           # has file been prepared already?
           if saltkey.found(keygain, struct[0]):
               message='%s has already been reduced' % img
               raise SaltError(message)



           # housekeeping keywords
           fname, hist=history(level=1, wrap=False, exclude=['images', 'outimages', 'outpref'])
           saltkey.housekeeping(struct[0],keyslot, 'Images have been slotmode reduced', hist)

           # write FITS file
           saltio.writefits(struct,oimg, clobber=clobber)
           saltio.closefits(struct)


# -----------------------------------------------------------
# fast SLOT mode reduction

def slot(struct,infile,dbspeed,dbrate,dbgain,dbnoise,dbbias,dbamp,xcoeff,gaindb,xtalkfile,
         logfile,verbose):

    import saltprint, saltkey, saltio, saltstat, time

# identify instrument

    instrume,keyprep,keygain,keybias,keyxtalk,keyslot,status = saltkey.instrumid(struct,infile,logfile)

# number of image HDU

    nextend = 0
    while (status == 0):
        try:
            struct[nextend+1].header['XTENSION']
            nextend += 1
        except:
            break
    nccds,status = saltkey.get('NCCDS',struct[0],infile,logfile)
    amplifiers = nccds * 2
    if (nextend%(amplifiers) != 0):
        message = '\nERROR -- SALTSLOT: Number of image extensions and'
        message += 'number of amplifiers are not consistent'
        status = saltprint.err(saltlog,message)
    status = saltkey.new('NSCIEXT',nextend,'Number of science extensions',struct[0],infile,logfile)
    status = saltkey.new('NEXTEND',nextend,'Number of data extensions',struct[0],infile,logfile)

# check image file and gain database are compatible

    if (status == 0):
        ngains = len(dbgain)
        if (int(max(dbamp)) != amplifiers):
            message  = '\nERROR -- SALTGSLOT: ' + infile + ' contains ' + str(amplifiers) + ' amplifiers'
            message += ', the gaindb file ' + gaindb + ' contains ' + str(max(dbamp)) + ' amplifiers'
            status = saltprint.err(logfile,message)

# check image file and cross talk database are compatible

    if (status == 0):
        if (len(xcoeff)-1 != amplifiers):
            message  = '\nERROR -- SALTSLOT: ' + infile + ' contains ' + str(amplifiers) + ' amplifiers'
            message += ', the cross talk file ' + xtalkfile + ' contains ' + str(len(xcoeff)-1) + ' amplifiers'
            status = saltprint.err(logfile,message)

# housekeeping keywords

    if (status == 0):
        status = saltkey.put('SAL-TLM',time.asctime(time.localtime()),struct[0],infile,logfile)
        status = saltkey.new(keyslot,time.asctime(time.localtime()),
                             'Data have been cleaned by SALTSLOT',struct[0],infile,logfile)

# keywords for image extensions

    for i in range(nextend):
        hdu = i + 1
        status = saltkey.new('EXTNAME','SCI','Extension name',struct[hdu],infile,logfile)
        status = saltkey.new('EXTVER',hdu,'Extension number',struct[hdu],infile,logfile)

# log coefficent table

    if (status == 0):
        message = '%30s %5s %4s %8s' % ('HDU','Gain','Bias','Xtalk')
        saltprint.log(logfile,'\n     ---------------------------------------------',verbose)
        saltprint.log(logfile,message,verbose)
        saltprint.log(logfile,'     ---------------------------------------------',verbose)

# loop over image extensions

    if (status == 0):
        for i in range(nextend/2):
            hdu = i * 2 + 1
            amplifier = hdu%amplifiers
            if (amplifier == 0): amplifier = amplifiers
            if (status == 0):
                value,status = saltkey.get('NAXIS1',struct[hdu],infile,logfile)
                naxis1 = int(value)
            if (status == 0):
                value,status = saltkey.get('NAXIS1',struct[hdu+1],infile,logfile)
                naxis2 = int(value)
            if (status == 0 and hdu == 1): biassec, status = saltkey.get('BIASSEC',struct[hdu],infile,logfile)
            if (status == 0 and hdu == 1):
                ranges = biassec.lstrip('[').rstrip(']').split(',')
                x1_1 = int(ranges[0].split(':')[0]) - 1
                x2_1 = int(ranges[0].split(':')[1]) - 1
                y1_1 = int(ranges[1].split(':')[0]) - 1
                y2_1 = int(ranges[1].split(':')[1]) - 1
            if (status == 0 and hdu == 1): biassec, status = saltkey.get('BIASSEC',struct[hdu+1],infile,logfile)
            if (status == 0 and hdu == 1):
                ranges = biassec.lstrip('[').rstrip(']').split(',')
                x1_2 = int(ranges[0].split(':')[0]) - 1
                x2_2 = int(ranges[0].split(':')[1]) - 1
                y1_2 = int(ranges[1].split(':')[0]) - 1
                y2_2 = int(ranges[1].split(':')[1]) - 1
            if (status == 0 and hdu == 1): datasec,status = saltkey.get('DATASEC',struct[hdu],infile,logfile)
            if (status == 0 and hdu == 1):
                ranges = datasec.lstrip('[').rstrip(']').split(',')
                dx1_1 = int(ranges[0].split(':')[0]) - 1
                dx2_1 = int(ranges[0].split(':')[1])
                dy1_1 = int(ranges[1].split(':')[0]) - 1
                dy2_1 = int(ranges[1].split(':')[1])
            if (status == 0 and hdu == 1): datasec,status = saltkey.get('DATASEC',struct[hdu+1],infile,logfile)
            if (status == 0 and hdu == 1):
                ranges = datasec.lstrip('[').rstrip(']').split(',')
                dx1_2 = int(ranges[0].split(':')[0]) - 1
                dx2_2 = int(ranges[0].split(':')[1])
                dy1_2 = int(ranges[1].split(':')[0]) - 1
                dy2_2 = int(ranges[1].split(':')[1])
            if (status == 0 and dx2_1 - dx1_1 != dx2_2 - dx1_2):
                message = 'ERROR -- SALTSLOT: HDUs '+infile
                message += '['+str(hdu)+'] and '+infile+'['+str(hdu+1)+']'
                message += ' have different dimensions'
                status = saltprint.err(logfile,message)

# read speed and gain of each exposure

            if (status == 0 and hdu == 1):
                gainset,status = saltkey.get('GAINSET',struct[0],infile,logfile)
                rospeed,status = saltkey.get('ROSPEED',struct[0],infile,logfile)
                if (rospeed == 'NONE'):
                    saltprint.log(logfile," ",verbose)
                    message = "ERROR -- SALTSLOT: Readout speed is 'NONE' in "
                    message += "primary keywords of " + infile
                    status = saltprint.err(logfile,message)

# read raw images

            if (status == 0):
                imagedata1,status = saltio.readimage(struct,hdu,logfile)
                imagedata2,status = saltio.readimage(struct,hdu+1,logfile)

# gain correction

            if (status == 0):
                for j in range(len(dbgain)):
                    if (gainset == dbrate[j] and rospeed == dbspeed[j] and amplifier == int(dbamp[j])):
                        try:
                            gain1 = float(dbgain[j])
                            imagedata1 *= gain1
                        except:
                            mesage = 'ERROR -- SALTSLOT: Cannot perform gain correction on image '
                            message += infile+'['+str(hdu)+']'
                            status = saltprint.err(logfile,message)
                    elif (gainset == dbrate[j] and rospeed == dbspeed[j] and amplifier + 1 == int(dbamp[j])):
                        try:
                            gain2 = float(dbgain[j])
                            imagedata2 *= gain2
                        except:
                            mesage = 'ERROR -- SALTSLOT: Cannot perform gain correction on image '
                            message += infile+'['+str(hdu+1)+']'
                            status = saltprint.err(logfile,message)

# crosstalk correction

            if (status == 0):
                revimage1 = imagedata1 * float(xcoeff[amplifier])
                revimage2 = imagedata2 * float(xcoeff[amplifier+1])
                for j in range(dx2_1-dx1_1+1):
                    imagedata1[:,j] -= revimage2[:,dx2_2-j-1]
                    imagedata2[:,j] -= revimage1[:,dx2_1-j-1]

# bias subtraction

            if (status == 0):
                overx_val_1 = []
                overx_val_2 = []
                for x in range(x1_1,x2_1+1):
                    list_1 = imagedata1[y1_1:y2_1,x] * 1.0
                    overx_val_1.append(saltstat.median(list_1,logfile))
                    overlevel_1 = saltstat.median(overx_val_1,logfile)
                for x in range(x1_2,x2_2+1):
                    list_2 = imagedata2[y1_2:y2_2,x] * 1.0
                    overx_val_2.append(saltstat.median(list_2,logfile))
                    overlevel_2 = saltstat.median(overx_val_2,logfile)
                imagedata1 -= overlevel_1
                imagedata2 -= overlevel_2

# trim overscan

            if (status == 0):
                imagedata1 = imagedata1[dy1_1:dy2_1,dx1_1:dx2_1]
                imagedata2 = imagedata2[dy1_2:dy2_2,dx1_2:dx2_2]
                datasec = '[1:'+str(dx2_1-dx1_1)+',1:'+str(dy2_1-dy1_1)+']'
                status = saltkey.put('DATASEC',datasec,struct[hdu],infile,logfile)
                status = saltkey.rem('BIASSEC',struct[hdu],infile,logfile)
                datasec = '[1:'+str(dx2_2-dx1_2)+',1:'+str(dy2_2-dy1_2)+']'
                status = saltkey.put('DATASEC',datasec,struct[hdu+1],infile,logfile)
                status = saltkey.rem('BIASSEC',struct[hdu+1],infile,logfile)

# log coefficient table

            if (status == 0):
                infilename = infile.split('/')
                infilename = infilename[len(infilename)-1]
                message = '%25s[%3d] %5.2f %4d %8.6f' % \
                    (infilename, hdu, gain1, overlevel_1, float(xcoeff[amplifier+1]))
                saltprint.log(logfile,message,verbose)
                message = '%25s[%3d] %5.2f %4d %8.6f' % \
                    (infilename, hdu+1, gain2, overlevel_2,float(xcoeff[amplifier]))
                saltprint.log(logfile,message,verbose)

# update image in HDU structure

            if (status == 0):
                struct,status = saltio.writeimage(struct,hdu,imagedata1,logfile)
                struct,status = saltio.writeimage(struct,hdu+1,imagedata2,logfile)

    return struct, status

# -----------------------------------------------------------
# main code
if not iraf.deftask('saltslot'):
   parfile = iraf.osfn("saltred$saltslot.par")
   t = iraf.IrafTaskFactory(taskname="saltslot",value=parfile,function=saltslot, pkgname='saltred')
