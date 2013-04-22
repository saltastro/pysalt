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

"""Will combined images distributed over multiple
extensions into a single image.
It will be able to handle multiple images in a single FITS file and will also update the header information with all of the relevant details.
"""

# Ensure python 2.5 compatibility
from __future__ import with_statement

from pyraf import iraf
from pyraf.iraf import pysalt

import saltsafekey
import saltsafeio
import saltstring
import salttran
import os
import time
import numpy as np
import pyfits
from saltsafelog import logging, history
from salterror import SaltIOError, SaltError

debug=True

def slotmerge(images,outimages,outpref,geomfile,clobber,logfile,verbose):

    with logging(logfile,debug) as log:
        # are the arguments defined
        saltsafeio.argdefined('images',images)
        saltsafeio.argdefined('geomfile',geomfile)
        saltsafeio.argdefined('logfile',logfile)

        # if the input file is a list, does it exist?
        if images[0] == '@':
            saltsafeio.listexists('Input',images)

        # parse list of input files
        infiles=saltsafeio.listparse('Raw image',images,'','','')

        # check input files exist
        saltsafeio.filesexist(infiles,'','r')

        # load output name list: @list, * and comma separated
        outimages = outimages.strip()
        outpref = outpref.strip()
        if len(outpref) == 0 and len(outimages) == 0:
            raise SaltIOError('Output file(s) not specified')

        # test output @filelist exists
        if len(outimages) > 0 and outimages[0] == '@':
            saltsafeio.listexists('Output',outimages)

        # parse list of output files
        outfiles=saltsafeio.listparse('Output image',outimages,outpref,infiles,'')

        # are input and output lists the same length?
        saltsafeio.comparelists(infiles,outfiles,'Input','output')

        # do the output files already exist?
        if not clobber:
            saltsafeio.filesexist(outfiles,'','w')

        # does CCD geometry definition file exist
        geomfilefile = geomfile.strip()
        saltsafeio.fileexists(geomfile)

        # read geometry definition file
        gap = 0
        xshift = [0, 0]
        yshift = [0, 0]
        rotation = [0, 0]

        gap, xshift, yshift, rotation=saltsafeio.readccdgeom(geomfile)
        for ro in rotation:
            if ro!=0:
                log.warning('SLOTMERGE currently ignores CCD rotation')

        # Begin processes each file
        for infile, outfile in zip(infiles, outfiles):
            # determine the name for the output file
            outpath = outfile.rstrip(os.path.basename(outfile))
            if (len(outpath) == 0):
                outpath = '.'

            # open each raw image
            struct=saltsafeio.openfits(infile)

            # identify instrument
            instrume,keyprep,keygain,keybias,keyxtalk,keyslot=saltsafekey.instrumid(struct,infile)

            # how many amplifiers?
            nccds=saltsafekey.get('NCCDS',struct[0],infile)
            amplifiers = nccds * 2
            #if (nccds != 2):
            #    raise SaltError('Can not currently handle more than two CCDs')

            # CCD geometry coefficients
            if instrume == 'RSS' or instrume == 'PFIS':
                xsh = [xshift[0], 0., xshift[1]]
                ysh = [yshift[0], 0., yshift[1]]
                rot = [rotation[0], 0., rotation[1]]
                refid = 1
            if instrume == 'SALTICAM':
                xsh = [xshift[0], 0.]
                ysh = [yshift[0], 0.]
                rot = [rotation[0], 0]
                refid = 1

            # how many extensions?
            nextend=saltsafekey.get('NEXTEND',struct[0],infile)

            # how many exposures
            exposures = nextend/amplifiers

            # CCD on-chip binning
            xbin, ybin=saltsafekey.ccdbin(struct[0],infile)
            gp = int(gap / xbin)

            # create output hdu structure
            outstruct = [None] * int(exposures+1)
            outstruct[0]=struct[0]

            # iterate over exposures, stitch them to produce file of CCD images
            for i in range(exposures):
                # Determine the total size of the image
                xsize=0
                ysize=0
                for j in range(amplifiers):
                    hdu=i*amplifiers+j+1
                    try:
                        xsize += len(struct[hdu].data[0])
                        if ysize < len(struct[hdu].data):
                            ysize=len(struct[hdu].data)
                    except:
                        msg='Unable to access extension %i ' % hdu
                        raise SaltIOError(msg)
                xsize += gp* (nccds-1)
                maxxsh, minxsh = determineshifts(xsh)
                maxysh, minysh = determineshifts(ysh)
                xsize += (maxxsh-minxsh)
                ysize += (maxysh-minysh)

                # Determine the x and y origins for each frame
                xdist=0
                ydist=0
                shid=0
                x0=np.zeros(amplifiers)
                y0=np.zeros(amplifiers)
                for j in range(amplifiers):
                    x0[j]=xdist+xsh[shid]-minxsh
                    y0[j]=ysh[shid]-minysh
                    hdu=i*amplifiers+j+1
                    darr=struct[hdu].data
                    xdist += len(darr[0])
                    if j%2==1:
                        xdist += gp
                        shid += 1

                # make the out image
                outarr=np.zeros((ysize, xsize), np.float64)

                # Embed each frame into the output array
                for j in range(amplifiers):
                    hdu=i*amplifiers+j+1
                    darr=struct[hdu].data
                    outarr=salttran.embed(darr, x0[j], y0[j], outarr)

                # Add the outimage to the output structure
                hdu=i*amplifiers+1
                outhdu=i+1
                outstruct[outhdu] = pyfits.ImageHDU(outarr)
                outstruct[outhdu].header=struct[hdu].header

                # Fix the headers in each extension
                datasec='[1:%4i,1:%4i]' % (xsize, ysize)
                saltsafekey.put('DATASEC',datasec, outstruct[outhdu], outfile)
                saltsafekey.rem('DETSIZE',outstruct[outhdu],outfile)
                saltsafekey.rem('DETSEC',outstruct[outhdu],outfile)
                saltsafekey.rem('CCDSEC',outstruct[outhdu],outfile)
                saltsafekey.rem('AMPSEC',outstruct[outhdu],outfile)

                # add housekeeping key words
                outstruct[outhdu]=addhousekeeping(outstruct[outhdu], outhdu, outfile)

            # close input FITS file
            saltsafeio.closefits(struct)

            # housekeeping keywords
            keymosaic='SLOTMERG'
            fname, hist=history(level=1, wrap=False)
            saltsafekey.housekeeping(struct[0],keymosaic,'Amplifiers have been mosaiced', hist)
            #saltsafekey.history(outstruct[0],hist)

            # this is added for later use by
            saltsafekey.put('NCCDS', 0.5, outstruct[0])
            saltsafekey.put('NSCIEXT', exposures, outstruct[0])
            saltsafekey.put('NEXTEND', exposures, outstruct[0])

            # write FITS file of mosaiced image
            outstruct=pyfits.HDUList(outstruct)
            saltsafeio.writefits(outstruct, outfile, clobber=clobber)

def determineshifts(shift):
    """Determine the maximum shifts in the data."""

    maxsh=0
    minsh=0
    for sh in shift:
        if sh > maxsh: maxsh=sh
        if sh < minsh: minsh=sh
    return maxsh, minsh

def addWCS(struct, xbin, ybin, outfile):
    """WCS keywords"""

    saltsafekey.new('CRPIX1',1,'WCS: X reference pixel',struct,outfile)
    saltsafekey.new('CRPIX2',1,'WCS: Y reference pixel',struct,outfile)
    saltsafekey.new('CRVAL1',float(xbin), 'WCS: X reference coordinate value',struct,outfile)
    saltsafekey.new('CRVAL2',float(ybin), 'WCS: Y reference coordinate value',struct,outfile)
    saltsafekey.new('CDELT1',float(xbin),'WCS: X pixel size',struct,outfile)
    saltsafekey.new('CDELT2',float(ybin),'WCS: Y pixel size',struct,outfile)
    saltsafekey.new('CTYPE1','pixel','X type',struct,outfile)
    saltsafekey.new('CTYPE2','pixel','Y type',struct,outfile)

    return struct

def addhousekeeping(struct, outhdu, outfile):
    """housekeeping keywords"""

    saltsafekey.put('SAL-TLM',time.asctime(time.localtime()),struct,outfile)
    saltsafekey.new('EXTNAME','SCI','Extension name',struct,outfile)
    saltsafekey.new('EXTVER',outhdu,'Extension number',struct,outfile)

    return struct

# -----------------------------------------------------------
# main code

parfile = iraf.osfn("slottools$slotmerge.par")
t = iraf.IrafTaskFactory(taskname="slotmerge",value=parfile,function=slotmerge, pkgname='slottools')
