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

"""Slotphot is a tool to do photometry on a serious of data obtained with 
   slotmode configuration of salticam.  After reading in the arguments, the 
   algorithm is the following:

#. Open th fits file
#. Subtact a smooth background
#. Centroid the comparison star to acount for drift
#. Photometry on Comparison and target star
#. Write data
#. Plot data [optional]
#. Write fits [optional]

One of the basic assumptions of this program are that
the target and comparison star are on the same amplifier

Updates:

20080616
    * Fixed a number of small bugs
    * Updated the handling of multiple amplifiers so that the program can now handle having the newfits as as an input file

20080910
    * Corrected the background subtraction. Now a pixel is replaced with the median of its row if its two bright
    * Added option to allow row-by-row median subtraction

20081027
    * Added limits to the finddrift to limit it to within a set number of pixels
    * Gives the users on local background subraction and annullur definition

20090609
    * Added some comments
"""

# Ensure python 2.5 compatibility
from __future__ import with_statement

import os
import sys
import numpy as np
import pyfits

from pyraf import iraf
from pyraf.iraf import pysalt
import saltsafekey
import saltsafeio
import slottool
from slotbackground import subbackground
from saltsafelog import logging
from salterror import SaltError, SaltIOError

debug=True

def slotphot(images,outfile,srcfile,newfits=None,phottype='square', 
             subbacktype='median',sigback=3,mbin=7,sorder=3,niter=5,sigdet=5,
             contpix=10,ampperccd=2,ignorexp=6,driftlimit=10.,finddrift=True,
             outtype='ascii',reltime=True,clobber=True,logfile='salt.log',
             verbose=True):
    """Perform photometry on listed SALT slotmode *images*."""

    with logging(logfile,debug) as log:
        # set up the variables
        entries = []
        vig_lo = {}
        vig_hi = {}
        amp = {}
        x = {}
        y = {}
        x_o = {}
        y_o = {}
        r = {}
        br1 = {}
        br2 = {}
        hour = 0
        min = 0
        sec = 0.
        time0 = 0.
        nframes = 0
        bin=mbin
        order=sorder

        # is the input file specified?
        saltsafeio.filedefined('Input',images)

        # if the input file is a list, does it exist?
        if images[0] == '@':
            saltsafeio.listexists('Input',images)

        # parse list of input files
        infiles=saltsafeio.listparse('Raw image',images,'','','')

        # check input files exist
        saltsafeio.filesexist(infiles,'','r')

        # is the output file specified?
        saltsafeio.filedefined('Output',outfile)

        # check output file does not exist, optionally remove it if it does exist
        if os.path.exists(outfile) and clobber:
            os.remove(outfile)
        elif os.path.exists(outfile) and not clobber:
            raise SaltIOError('File '+outfile+' already exists, use clobber=y')

        # open output ascii file
        if outtype=='ascii':
            try:
                lc = open(outfile,'a')
            except:
                raise SaltIOError('Cannot open ouput file '+outfile)

        # is the extraction region defintion file specified?
        saltsafeio.filedefined('Extraction region defintion',srcfile)

        # check extraction region defintion file exists
        srcfile = srcfile.strip()
        saltsafeio.fileexists(srcfile)

        # read extraction region defintion file
        amp, x, y, x_o, y_o, r, br1, br2=slottool.readsrcfile(srcfile)

        # set the writenewfits parameter
        if not newfits or newfits=='none':
            writenewfits=False
        else:
            writenewfits=newfits

        # get time of first exposure and basic information about the observations
        infile=infiles[0]
        struct=saltsafeio.openfits(infile)

        # identify instrument
        instrume,keyprep,keygain,keybias,keyxtalk,keyslot=saltsafekey.instrumid(struct,infile)

        # how many extensions?
        nextend=saltsafekey.get('NEXTEND',struct[0],infile)
        if nextend < amp['comparison']:
            msg='Insufficient number of extensions in %s' % (infile)
            raise SaltIOError(msg)

        # how many amplifiers?
        amplifiers=saltsafekey.get('NCCDS',struct[0],infile)
        amplifiers = int(ampperccd*float(amplifiers))
        if ampperccd>0:
            nframes = int(nextend/amplifiers)
            nstep=amplifiers
        else:
            nframes = nextend
            nstep=1
        ntotal=nframes*len(infiles)

        # image size
        naxis1=saltsafekey.get('NAXIS1',struct[amp['comparison']],infile)
        naxis2=saltsafekey.get('NAXIS2',struct[amp['comparison']],infile)

        # CCD binning
        ccdsum=saltsafekey.get('CCDSUM',struct[0],infile)
        binx=int(ccdsum.split(' ')[0])
        biny=int(ccdsum.split(' ')[1])

        # Identify the time of the observations
        ext = 1
        try:
            time0=slottool.getobstime(struct[ext], infile+'['+str(ext)+']')
            dateobs=saltsafekey.get('DATE-OBS',struct[ext],infile)
            dateobs=dateobs.replace('-','/')
        except:
            raise SaltIOError('No time or obsdate in first image')

        # If a total file is to be written out, create it and update it
        if writenewfits:
            if os.path.isfile(writenewfits):
                if clobber:
                    saltsafeio.delete(writenewfits)
                else:
                    raise SaltIOError('Newfits file exists, use clobber')

            try:
                hdu=pyfits.PrimaryHDU()
                hdu.header=struct[0].header
                hdu.header['NCCDS']=1
                hdu.header['NSCIEXT']=ntotal-ignorexp
                hdu.header['NEXTEND']=ntotal-ignorexp
                hduList=pyfits.HDUList(hdu)
                hduList.verify()
                hduList.writeto(writenewfits)
            except:
                raise SaltIOError('Could not create newfits file, '+writenewfits)

        # Close image file
        saltsafeio.closefits(struct)

        # Read newfits file back in for updating
        if writenewfits:
            try:
                hduList=pyfits.open(writenewfits,mode='update')
            except:
                raise SaltIOError('Cannot open newfits file '+writenewfits+' for updating.')

        # set up the arrays
        j=0
        time=np.zeros(ntotal ,dtype='float')-1.0
        dx=np.zeros(ntotal ,dtype='float')
        dy=np.zeros(ntotal ,dtype='float')
        tflux=np.zeros(ntotal ,dtype='float')
        terr =np.zeros(ntotal ,dtype='float')
        cflux=np.zeros(ntotal ,dtype='float')
        cerr =np.zeros(ntotal ,dtype='float')
        ratio=np.zeros(ntotal ,dtype='float')
        rerr =np.zeros(ntotal ,dtype='float')
        tgt_x=np.zeros(ntotal ,dtype='float')
        tgt_y=np.zeros(ntotal ,dtype='float')
        cmp_x=np.zeros(ntotal ,dtype='float')
        cmp_y=np.zeros(ntotal ,dtype='float')

        p_one=100./ntotal            # One percent
        p_old=-1                     # Previous completed percentage
        p_new=0                      # New completed percentage
        p_n=1                        # Counter number
        for infile in infiles:
            # Log
            if verbose:
                log.message('Starting photometry on file '+infile, with_stdout=False)

            struct=pyfits.open(infile)

            # Skip through the frames and process each frame individually
            for i in range(nframes):
                # Show progress
                if verbose:
                    p_new=int(p_n*p_one)
                    p_n+=1
                    if p_new!=p_old:
                        ctext='Percentage Complete: %d\r' % p_new
                        sys.stdout.write(ctext)
                        sys.stdout.flush()
                        p_old=p_new

                if not (infile==infiles[0] and i < ignorexp):

                    ext=amp['comparison']+i*nstep
                    try:
                        header=struct[ext].header
                        array=struct[ext].data
                        array=array*1.0
                    except:
                        msg='Unable to open extension %i in image %s' % (ext, infile)
                        raise SaltIOError(msg)

                    # starti the analysis of each frame
                    # get the time
                    time[j]=slottool.getobstime(struct[ext],infile+'['+str(ext)+']')

                    # gain and readout noise
                    try:
                        gain=float(header['GAIN'])
                    except:
                        gain=1
                        raise SaltIOError('Gain not specified in image header')
                    try:
                        rdnoise=float(header['RDNOISE'])
                    except:
                        rdnoise=0
                        raise SaltIOError('RDNOISE not specified in image header')

                    # background subtraction
                    if not subbacktype=='none':
                        try:
                            array=subbackground(array, sigback, bin, order, niter, subbacktype)
                        except SaltError:
                            log.warning('Image '+infile+' extention '+str(ext)+' is blank, skipping')
                            continue

                    # x-y fit to the comparison star and update the x,y values
                    if finddrift:
                        carray, fx,fy=slottool.finddrift(array, x['comparison'], y['comparison'], r['comparison'], naxis1, naxis2, sigdet, contpix, sigback, driftlimit, niter)
                        if fx > -1  and fy > -1:
                            if fx < naxis1 and fy < naxis2:
                                dx[j]=x['comparison']-fx
                                dy[j]=y['comparison']-fy
                                x['comparison']=fx
                                y['comparison']=fy
                                x['target']=x['target']-dx[j]
                                y['target']=y['target']-dy[j]
                            else:
                                dx[j]=0
                                dy[j]=0
                                x['comparison']=x_o['comparison']
                                y['comparison']=y_o['comparison']
                                x['target']=x_o['target']
                                y['target']=y_o['target']
                        else:
                            msg='No comparison object found in image file ' + infile+' on extension %i skipping.' % ext
                            log.warning(msg)
                            pass

                    # do photometry
                    try:
                        tflux[j],terr[j],cflux[j],cerr[j],ratio[j],rerr[j]=slottool.dophot(phottype, array, x, y, r, br1, br2, gain, rdnoise, naxis1, naxis2)
                    except SaltError, e:
                        msg='Could not do photometry on extension %i in image %s because %s skipping.' % (ext, infile, e)
                        log.warning(msg)

                    tgt_x[j]=x['target']
                    tgt_y[j]=y['target']
                    cmp_x[j]=x['comparison']
                    cmp_y[j]=y['comparison']

                    # record results
                    # TODO! This should be removed in favor of the write all to disk in the end
                    if outtype=='ascii':
                        slottool.writedataout(lc, j+1, time[j], x, y, tflux[j], terr[j], cflux[j],cerr[j],ratio[j],rerr[j],time0, reltime)

                    # write newfits file
                    if writenewfits:
                        # add original name and extension number to header
                        try:
                            hdue=pyfits.ImageHDU(array)
                            hdue.header=header
                            hdue.header.update('ONAME',infile,'Original image name')
                            hdue.header.update('OEXT',ext,'Original extension number')
                            hduList.append(hdue)
                        except:
                            log.warning('Could not update image in newfits '+infile+' '+str(ext))

                    # increment counter
                    j+=1

            # close FITS file
            saltsafeio.closefits(struct)

        # close newfits file
        if writenewfits:
            try:
                hduList.flush()
                hduList.close()
            except:
                raise SaltIOError('Cannot close newfits file.')

        # write to output
        if outtype=='ascii':
            # close output ascii file
            try:
                lc.close()
            except:
                raise SaltIOError('Cannot close ouput file ' + outfile)
        elif outtype=='fits':
            print 'writing fits'
            try:
                c1=pyfits.Column(name='index',format='D',array=np.arange(ntotal))
                if reltime:
                    c2=pyfits.Column(name='time',format='D',array=time-time0)
                else:
                    c2=pyfits.Column(name='time',format='D',array=time)

                c3=pyfits.Column(name='tgt_x',format='D',array=tgt_x)
                c4=pyfits.Column(name='tgt_y',format='D',array=tgt_y)
                c5=pyfits.Column(name='tgt_flux',format='D',array=tflux)
                c6=pyfits.Column(name='tgt_err',format='D',array=terr)
                c7=pyfits.Column(name='cmp_x',format='D',array=cmp_x)
                c8=pyfits.Column(name='cmp_y',format='D',array=cmp_y)
                c9=pyfits.Column(name='cmp_flux',format='D',array=cflux)
                c10=pyfits.Column(name='cmp_err',format='D',array=cerr)
                c11=pyfits.Column(name='flux_ratio',format='D',array=ratio)
                c12=pyfits.Column(name='flux_ratio_err',format='D',array=rerr)
                tbhdu=pyfits.new_table([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12])
                # Add header information
                tbhdu.header.update('RELTIME',str(reltime),'Time relative to first datapoint or absolute.')

                tbhdu.writeto(outfile)
                print 'fits written to ',outfile
            except:
                raise SaltIOError('Could not write to fits table.')

# -----------------------------------------------------------
# main code

parfile = iraf.osfn("slottools$slotphot.par")
t = iraf.IrafTaskFactory(taskname="slotphot",value=parfile,function=slotphot,pkgname='slottools')
