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

"""A tool to produce to produce a fits file which only contains the fits extensions of interest.
If selected, the tool will either create a background
subtracted image or the background image.
"""

# Ensure python 2.5 compatibility
from __future__ import with_statement

import os
import sys
from pyraf import iraf
from pyraf.iraf import pysalt
import saltsafekey as saltkey
import saltsafeio as saltio
from slotbackground import subbackground
from saltsafelog import logging
from salterror import SaltIOError
import pyfits

debug=True

def slotback(images,outfits,extension,imgtype='image',subbacktype='median',
             sigback=3,mbin=7,sorder=3,niter=5,ampperccd=2,ignorexp=6,
             clobber=False,logfile='salt.log',verbose=True):

    with logging(logfile,debug) as log:
        # set up the variables
        order=sorder
        plotfreq=0.01
        ftime=plotfreq
    
        # is the input file specified?
        infiles = saltio.argunpack ('Input',images)

        # is the output file specified?
        saltio.filedefined('Output',outfits)
    
        #open the first file and check its charactistics
        struct=saltio.openfits(infiles[0])

        # how many extensions?
        nextend=saltkey.get('NEXTEND',struct[0])
        if nextend < extension:
            msg='Insufficient number of extensions in %s' % (infile)
            raise SaltIOError(msg)
    
        # how many amplifiers?
        amplifiers=saltkey.get('NCCDS',struct[0])
        amplifiers = int(ampperccd*float(amplifiers))
        if ampperccd>0:
            nframes = nextend/amplifiers
            nstep=amplifiers
        else:
            nframes = nextend
            nstep=1
        ntotal=nframes*len(infiles)
    
        # image size
        naxis1=saltkey.get('NAXIS1',struct[extension])
        naxis2=saltkey.get('NAXIS2',struct[extension])
    
        # CCD binning
        ccdsum=saltkey.get('CCDSUM',struct[0])
        binx=int(ccdsum.split(' ')[0])
        biny=int(ccdsum.split(' ')[1])
    
        # If a total file is to written out, create it and update it
        hdu = 1

        # delete the file if clobber is on and the file exists
        if os.path.isfile(outfits):
           if clobber:
               saltio.delete(outfits)
           else:
               raise SaltIOError('Output fits file '+writenewfits+'already exists. Use Clobber=yes to overwrite file')
    
        # create the output file
        try:
           hdu = pyfits.PrimaryHDU()
           hdu.header=struct[0].header
           hdu.header['NCCDS']=1
           hdu.header['NSCIEXT']=ntotal-ignorexp
           hdu.header['NEXTEND']=ntotal-ignorexp
           hduList = pyfits.HDUList(hdu)
           hduList.verify()
           hduList.writeto(outfits)
        except:
           raise SaltIOError('Could not create new fits file named '+writenewfits)
    
        # read it back in for updating
        hduList = saltio.openfits(outfits,mode='update')
    
        # set up the completeness tracker
        if verbose:
            j=0
            x2=float(j)/float(ntotal)
            ctext='Percentage Complete: %3.2f' % x2
            sys.stdout.write(ctext)
            sys.stdout.flush()
    
        for infile in infiles:
           struct=saltio.openfits(infile)

           # Skip through the frames and process each frame individually
           for i in range(nframes):
                if not (infile==infiles[0] and i < ignorexp):
                    hdu=extension+i*nstep
                    try:
                        header=struct[hdu].header
                        array=struct[hdu].data
                        array=array*1.0
                    except Exception, e:
                        msg='Unable to open extension %i in image %s because %s' % (hdu, infile, e)
                        raise SaltIOError(msg)
    
                    # start the analysis of each frame
                    # gain and readout noise
                    try:
                        gain=float(header['GAIN'])
                    except:
                        gain=1
                        log.warning('Gain not specified in image header')
                    try:
                        rdnoise=float(header['RDNOISE'])
                    except:
                        rdnoise=0
                        log.warning('RDNOISE not specified in image header')

                    # background subtraction
                    if not subbacktype=='none':
                       try:
                           narray=subbackground(array, sigback, mbin, order, niter, subbacktype)
                       except:
                           log.warning('Image '+infile+' extention '+str(i)+' is blank, skipping')
                           continue


                    # create output array
                    try:
                       if imgtype=='background':
                           oarray=narray-array
                       else:
                           oarray=narray
                    except:
                       oarray=array
    
                    # print progress
                    if verbose:
                       x2=float(j)/float(ntotal)
                       if x2 > ftime:
                           ctext='\b\b\b\b\b %3.2f' % x2
                           sys.stdout.write(ctext)
                           sys.stdout.flush()
                           ftime += plotfreq
    

                    # update the header values with the image name and extension number
                    try:
                       hdue = pyfits.ImageHDU(oarray)
                       hdue.header=header
                       hdue.header.update('ONAME',infile,'Original image name')
                       hdue.header.update('OEXT',hdu,'Original extension number')
                    except Exception, e:
                       msg='SALTPHOT--WARNING:  Could not update image in newfits file for %s ext %i because %s' \
                           % (infile, hdu, e)
                       raise SaltIOError(msg)

                    hduList.append(hdue)
                    j +=1
    
           # close FITS file
           saltio.closefits(struct)
    
        # close the output fits file:
        try:
           # write out the file
           hduList.flush()
           hduList.close()
        except Exception, e:
           raise SaltIOError('Failed to write %s because %s' % (outfits, e))
    
# -----------------------------------------------------------
# main code

parfile = iraf.osfn("slottools$slotback.par")
t = iraf.IrafTaskFactory(taskname="slotback",value=parfile,function=slotback,pkgname='slottools')
