#!/usr/bin/env python

# LICENSE
# Copyright (c) 2014, South African Astronomical Observatory (SAAO)
# All rights reserved.  See License file for more details

"""
SALTMOSAIC is a task to apply the CCD geometric corrections to MEF style SALT
data.

Author                 Version      Date
-----------------------------------------------
Martin Still (SAAO)    0.1          16 Oct 2006
SM Crawford (SAAO)      0.2           19 Mar 2006

Updates
--------------------
20120319   - Update to new error handling
           - Changed the mosaic to use the whole frame and not trim some data
             off
20141111   - Added option to replace the masked regions
"""

import os
import time
import numpy
from scipy import ndimage as nd
from astropy.io import fits
from pyraf import iraf

from math import cos, sin, pi
from scipy.ndimage import geometric_transform

import saltsafekey as saltkey
import saltsafeio as saltio
import saltsafestring as saltstring
from saltsafelog import logging, history

from salterror import SaltError

debug = True


# -----------------------------------------------------------
# core routine

def saltmosaic(images, outimages, outpref, geomfile, interp='linear',
               geotran=True, fill=False, cleanup=True, clobber=False,
               logfile=None, verbose=True):

    # Start the logging
    with logging(logfile, debug) as log:

        # Check the input images
        infiles = saltio.argunpack('Input', images)

        # create list of output files
        outfiles = saltio.listparse('Outfile', outimages, outpref, infiles, '')

        # verify that the input and output lists are the same length
        saltio.comparelists(infiles, outfiles, 'Input', 'output')

        # does CCD geometry definition file exist
        geomfilefile = geomfile.strip()
        saltio.fileexists(geomfile)

        gap = 0
        xshift = [0, 0]
        yshift = [0, 0]
        rotation = [0, 0]
        gap, xshift, yshift, rotation = saltio.readccdgeom(geomfile)

        # open each raw image file and apply the transformation to it
        for img, oimg in zip(infiles, outfiles):

            # open the structure
            struct = saltio.openfits(img)

            # create the mosaic
            ostruct = make_mosaic(
                struct,
                gap,
                xshift,
                yshift,
                rotation,
                interp_type=interp,
                geotran=geotran,
                fill=fill,
                cleanup=cleanup,
                log=log,
                verbose=verbose)

            # update the header information
            # housekeeping keywords
            fname, hist = history(
                level=1, wrap=False, exclude=[
                    'images', 'outimages', 'outpref'])
            saltkey.housekeeping(
                ostruct[0],
                'SMOSAIC',
                'Images have been mosaicked',
                hist)

            # write the image out
            ostruct.writeto(oimg, clobber=clobber, output_verify='ignore')

            # close the files
            struct.close()
            ostruct.close()


def make_mosaic(struct, gap, xshift, yshift, rotation, interp_type='linear',
                boundary='constant', constant=0, geotran=True, fill=False,
                cleanup=True, log=None, verbose=False):
    """Given a SALT image struct, combine each of the individual amplifiers and
        apply the geometric CCD transformations to the image
    """

    # get the name of the file
    infile = saltkey.getimagename(struct[0], base=True)
    outpath = './'

    # identify instrument
    instrume, keyprep, keygain, keybias, keyxtalk, keyslot = \
        saltkey.instrumid(struct)

    # how many amplifiers?
    nsciext = saltkey.get('NSCIEXT', struct[0])
    nextend = saltkey.get('NEXTEND', struct[0])
    nccds = saltkey.get('NCCDS', struct[0])
    amplifiers = nccds * 2

    if nextend > nsciext:
        varframe = True
    else:
        varframe = False

    # CCD geometry coefficients
    if (instrume == 'RSS' or instrume == 'PFIS'):
        xsh = [0., xshift[0], 0., xshift[1]]
        ysh = [0., yshift[0], 0., yshift[1]]
        rot = [0., rotation[0], 0., rotation[1]]
    elif instrume == 'SALTICAM':
        xsh = [0., xshift[0], 0.]
        ysh = [0., yshift[0], 0.]
        rot = [0., rotation[0], 0]

    # how many extensions?
    nextend = saltkey.get('NEXTEND', struct[0])

    # CCD on-chip binning
    xbin, ybin = saltkey.ccdbin(struct[0])

    # create temporary primary extension
    outstruct = []
    outstruct.append(struct[0])
    # define temporary FITS file store tiled CCDs

    tilefile = saltio.tmpfile(outpath)
    tilefile += 'tile.fits'
    if varframe:
        tilehdu = [None] * (3 * int(nsciext / 2) + 1)
    else:
        tilehdu = [None] * int(nsciext / 2 + 1)
    tilehdu[0] = fits.PrimaryHDU()
    #tilehdu[0].header = struct[0].header

    if log:
        log.message('', with_stdout=verbose)

    # iterate over amplifiers, stich them to produce file of CCD images
    for i in range(int(nsciext / 2)):
        hdu = i * 2 + 1
        # amplifier = hdu%amplifiers
        # if (amplifier == 0): amplifier = amplifiers

        # read DATASEC keywords
        datasec1 = saltkey.get('DATASEC', struct[hdu])
        datasec2 = saltkey.get('DATASEC', struct[hdu + 1])
        xdsec1, ydsec1 = saltstring.secsplit(datasec1)
        xdsec2, ydsec2 = saltstring.secsplit(datasec2)

        # read images
        imdata1 = saltio.readimage(struct, hdu)
        imdata2 = saltio.readimage(struct, hdu + 1)

        # tile 2n amplifiers to yield n CCD images
        outdata = numpy.zeros((int(ydsec1[1] + abs(ysh[i + 1] / ybin)),
                               int(xdsec1[1] + xdsec2[1] +
                               abs(xsh[i + 1] / xbin))), 
                               numpy.float32)

        # set up the variance frame
        if varframe:
            vardata = outdata.copy()
            vdata1 = saltio.readimage(struct, struct[hdu].header['VAREXT'])
            vdata2 = saltio.readimage(struct, struct[hdu + 1].header['VAREXT'])

            bpmdata = outdata.copy()
            bdata1 = saltio.readimage(struct, struct[hdu].header['BPMEXT'])
            bdata2 = saltio.readimage(struct, struct[hdu + 1].header['BPMEXT'])

        x1 = xdsec1[0] - 1
        if x1 != 0:
            msg = 'The data in %s have not been trimmed prior to mosaicking.' \
                  % infile
            log.error(msg)
        if xsh[i + 1] < 0:
            x1 += int(abs(xsh[i + 1] / xbin))
        x2 = x1 + xdsec1[1]
        y1 = ydsec1[0] - 1
        if ysh[i + 1] < 0:
            y1 += int(abs(ysh[i + 1] / ybin))
        y2 = y1 + ydsec1[1]
        outdata[y1:y2, x1:x2] =\
            imdata1[ydsec1[0] - 1:ydsec1[1], xdsec1[0] - 1:xdsec1[1]]

        if varframe:
            vardata[y1:y2, x1:x2] =\
                vdata1[ydsec1[0] - 1:ydsec1[1], xdsec1[0] - 1:xdsec1[1]]
            bpmdata[y1:y2, x1:x2] =\
                bdata1[ydsec1[0] - 1:ydsec1[1], xdsec1[0] - 1:xdsec1[1]]

        x1 = x2
        x2 = x1 + xdsec2[1]
        y1 = ydsec2[0] - 1
        if ysh[i + 1] < 0:
            y1 += abs(ysh[i + 1] / ybin)
        y2 = y1 + ydsec2[1]
        outdata[y1:y2, x1:x2] =\
            imdata2[ydsec1[0] - 1:ydsec1[1], xdsec1[0] - 1:xdsec1[1]]

        if varframe:
            vardata[y1:y2, x1:x2] =\
                vdata2[ydsec1[0] - 1:ydsec1[1], xdsec1[0] - 1:xdsec1[1]]
            bpmdata[y1:y2, x1:x2] =\
                bdata2[ydsec1[0] - 1:ydsec1[1], xdsec1[0] - 1:xdsec1[1]]

        # size of new image
        naxis1 = str(xdsec1[1] + xdsec2[1])
        naxis2 = str(ydsec1[1])

        # add image and keywords to HDU list
        tilehdu[i + 1] = fits.ImageHDU(outdata)
        tilehdu[i + 1].header = struct[hdu].header
        #tilehdu[
        #    i + 1].header['DATASEC'] = '[1:' + naxis1 + ',1:' + naxis2 + ']'

        if varframe:
            vext = i + 1 + int(nsciext / 2.)
            tilehdu[vext] = fits.ImageHDU(vardata)
            #tilehdu[vext].header = struct[struct[hdu].header['VAREXT']].header
            #tilehdu[vext].header[
            #    'DATASEC'] = '[1:' + naxis1 + ',1:' + naxis2 + ']'

            bext = i + 1 + 2 * int(nsciext / 2.)
            tilehdu[bext] = fits.ImageHDU(bpmdata)
            #tilehdu[bext].header = struct[struct[hdu].header['BPMEXT']].header
            #tilehdu[bext].header[
            #    'DATASEC'] = '[1:' + naxis1 + ',1:' + naxis2 + ']'

        # image tile log message #1
        if log:
            message = os.path.basename(infile) + '[' + str(hdu) + ']['
            message += str(xdsec1[0]) + ':' + str(xdsec1[1]) + ','
            message += str(ydsec1[0]) + ':' + str(ydsec1[1]) + '] --> '
            message += os.path.basename(tilefile) + '[' + str(i + 1) + ']['
            message += str(xdsec1[0]) + ':' + str(xdsec1[1]) + ','
            message += str(ydsec1[0]) + ':' + str(ydsec1[1]) + ']'
            log.message(message, with_stdout=verbose, with_header=False)
            message = os.path.basename(infile) + '[' + str(hdu + 1) + ']['
            message += str(xdsec1[0]) + ':' + str(xdsec1[1]) + ','
            message += str(ydsec1[0]) + ':' + str(ydsec1[1]) + '] --> '
            message += os.path.basename(tilefile) + '[' + str(i + 1) + ']['
            message += str(xdsec1[1] + 1) + ':' + \
                str(xdsec1[1] + xdsec2[1]) + ','
            message += str(ydsec2[0]) + ':' + str(ydsec2[1]) + ']'
            log.message(message, with_stdout=verbose, with_header=False)

    # write temporary file of tiled CCDs
    hdulist = fits.HDUList(tilehdu)
    hdulist.writeto(tilefile)

    # iterate over CCDs, transform and rotate images
    yrot = [None] * 4
    xrot = [None] * 4

    tranfile = [' ']
    tranhdu = [0]
    if varframe:
        tranfile = [''] * (3 * int(nsciext / 2) + 1)
        tranhdu = [0] * (3 * int(nsciext / 2) + 1)
    else:
        tranfile = [''] * int(nsciext / 2 + 1)
        tranhdu = [0] * int(nsciext / 2 + 1)

    # this is hardwired for SALT where the second CCD is considered the
    # fiducial
    for hdu in range(1, int(nsciext / 2 + 1)):
        tranfile[hdu] = saltio.tmpfile(outpath)
        tranfile[hdu] += 'tran.fits'
        if varframe:
            tranfile[hdu + nccds] = saltio.tmpfile(outpath) + 'tran.fits'
            tranfile[hdu + 2 * nccds] = saltio.tmpfile(outpath) + 'tran.fits'

        ccd = hdu % nccds
        if (ccd == 0):
            ccd = nccds

        # correct rotation for CCD binning
        yrot[ccd] = rot[ccd] * ybin / xbin
        xrot[ccd] = rot[ccd] * xbin / ybin
        dxshift = xbin * int(float(int(gap) / xbin) + 0.5) - gap

        # transformation using geotran IRAF task
        # if (ccd == 1):
        if (ccd != 2):

            if geotran:
                message = '\nSALTMOSAIC -- geotran ' + tilefile + \
                    '[' + str(ccd) + '] ' + tranfile[hdu]
                message += ' \"\" \"\" xshift=' + \
                    str((xsh[ccd] + (2 - ccd) * dxshift) / xbin) + ' '
                message += 'yshift=' + \
                    str(ysh[ccd] / ybin) + ' xrotation=' + str(xrot[ccd]) + ' '
                message += 'yrotation=' + \
                    str(yrot[ccd]) + ' xmag=1 ymag=1 xmin=\'INDEF\''
                message += 'xmax=\'INDEF\' ymin=\'INDEF\' ymax=\'INDEF\' '
                message += 'ncols=\'INDEF\' '
                message += 'nlines=\'INDEF\' verbose=\'no\' '
                message += 'fluxconserve=\'yes\' nxblock=2048 '
                message += 'nyblock=2048 interpolant=\'' + \
                    interp_type + '\' boundary=\'constant\' constant=0'
                log.message(message, with_stdout=verbose)

                yd, xd = tilehdu[ccd].data.shape
                ncols = 'INDEF'  # ncols=xd+abs(xsh[ccd]/xbin)
                nlines = 'INDEF'  # nlines=yd+abs(ysh[ccd]/ybin)
                geo_xshift = xsh[ccd] + (2 - ccd) * dxshift / xbin
                geo_yshift = ysh[ccd] / ybin
                iraf.images.immatch.geotran(tilefile + "[" + str(ccd) + "]",
                                            tranfile[hdu],
                                            "",
                                            "",
                                            xshift=geo_xshift,
                                            yshift=geo_yshift,
                                            xrotation=xrot[ccd],
                                            yrotation=yrot[ccd],
                                            xmag=1, ymag=1, xmin='INDEF',
                                            xmax='INDEF', ymin='INDEF',
                                            ymax='INDEF', ncols=ncols,
                                            nlines=nlines, verbose='no',
                                            fluxconserve='yes', nxblock=2048,
                                            nyblock=2048, interpolant="linear",
                                            boundary="constant", constant=0)
                if varframe:
                    var_infile = tilefile + "[" + str(ccd + nccds) + "]"
                    iraf.images.immatch.geotran(var_infile,
                                                tranfile[hdu + nccds],
                                                "",
                                                "",
                                                xshift=geo_xshift,
                                                yshift=geo_yshift,
                                                xrotation=xrot[ccd],
                                                yrotation=yrot[ccd],
                                                xmag=1, ymag=1, xmin='INDEF',
                                                xmax='INDEF', ymin='INDEF',
                                                ymax='INDEF', ncols=ncols,
                                                nlines=nlines, verbose='no',
                                                fluxconserve='yes',
                                                nxblock=2048, nyblock=2048,
                                                interpolant="linear",
                                                boundary="constant",
                                                constant=0)
                    var2_infile = tilefile + "[" + str(ccd + 2 * nccds) + "]"
                    iraf.images.immatch.geotran(var2_infile,
                                                tranfile[hdu + 2 * nccds],
                                                "",
                                                "",
                                                xshift=geo_xshift,
                                                yshift=geo_yshift,
                                                xrotation=xrot[ccd],
                                                yrotation=yrot[ccd],
                                                xmag=1, ymag=1, xmin='INDEF',
                                                xmax='INDEF', ymin='INDEF',
                                                ymax='INDEF', ncols=ncols,
                                                nlines=nlines, verbose='no',
                                                fluxconserve='yes',
                                                nxblock=2048, nyblock=2048,
                                                interpolant="linear",
                                                boundary="constant",
                                                constant=0)

                # open the file and copy the data to tranhdu
                tstruct = fits.open(tranfile[hdu])
                tranhdu[hdu] = tstruct[0].data
                tstruct.close()
                if varframe:
                    tranhdu[
                        hdu +
                        nccds] = fits.open(
                        tranfile[
                            hdu +
                            nccds])[0].data
                    tranhdu[
                        hdu +
                        2 *
                        nccds] = fits.open(
                        tranfile[
                            hdu +
                            2 *
                            nccds])[0].data

            else:
                log.message(
                    "Transform CCD #%i using dx=%s, dy=%s, rot=%s" %
                    (ccd,
                     xsh[ccd] /
                        2.0,
                        ysh[ccd] /
                        2.0,
                        xrot[ccd]),
                    with_stdout=verbose,
                    with_header=False)
                tranhdu[hdu] = geometric_transform(
                    tilehdu[ccd].data,
                    tran_func,
                    prefilter=False,
                    order=1,
                    extra_arguments=(
                        xsh[ccd] / 2,
                        ysh[ccd] / 2,
                        1,
                        1,
                        xrot[ccd],
                        yrot[ccd]))
                tstruct = fits.PrimaryHDU(tranhdu[hdu])
                tstruct.writeto(tranfile[hdu])
                if varframe:
                    tranhdu[hdu + nccds] = geometric_transform(
                        tilehdu[hdu + 3].data,
                        tran_func,
                        prefilter=False,
                        order=1,
                        extra_arguments=(
                            xsh[ccd] / 2, ysh[ccd] / 2,
                            1, 1,
                            xrot[ccd], yrot[ccd]))
                    tranhdu[hdu + 2 * nccds] = geometric_transform(
                        tilehdu[hdu + 6].data,
                        tran_func,
                        prefilter=False,
                        order=1,
                        extra_arguments=(
                            xsh[ccd] / 2, ysh[ccd] / 2,
                            1, 1,
                            xrot[ccd], yrot[ccd]))

        else:
            log.message(
                "Transform CCD #%i using dx=%s, dy=%s, rot=%s" %
                (ccd, 0, 0, 0), with_stdout=verbose, with_header=False)
            tranhdu[hdu] = tilehdu[ccd].data
            if varframe:
                tranhdu[hdu + nccds] = tilehdu[ccd + nccds].data
                tranhdu[hdu + 2 * nccds] = tilehdu[ccd + 2 * nccds].data

    # open outfile
    if varframe:
        outlist = 4 * [None]
    else:
        outlist = 2 * [None]

    #outlist[0] = struct[0].copy()
    outlist[0] = fits.PrimaryHDU()
    outlist[0].header = struct[0].header

    naxis1 = int(gap / xbin * (nccds - 1))
    naxis2 = 0
    for i in range(1, nccds + 1):
        yw, xw = tranhdu[i].shape
        naxis1 += xw + int(abs(xsh[ccd] / xbin)) + 1
        naxis2 = max(naxis2, yw)
    outdata = numpy.zeros((naxis2, naxis1), numpy.float32)
    outdata.shape = naxis2, naxis1
    if varframe:
        vardata = outdata * 0
        bpmdata = outdata * 0 + 1

    # iterate over CCDs, stich them to produce a full image
    hdu = 0
    totxshift = 0
    for hdu in range(1, nccds + 1):

        # read DATASEC keywords
        ydsec, xdsec = tranhdu[hdu].shape

        # define size and shape of final image
        # tile CCDs to yield mosaiced image
        x1 = int((hdu - 1) * (xdsec + gap / xbin)) + int(totxshift)
        x2 = xdsec + x1
        y1 = int(0)
        y2 = int(ydsec)
        outdata[y1:y2, x1:x2] = tranhdu[hdu]
        totxshift += int(abs(xsh[hdu] / xbin)) + 1
        if varframe:
            vardata[y1:y2, x1:x2] = tranhdu[hdu + nccds]
            bpmdata[y1:y2, x1:x2] = tranhdu[hdu + 2 * nccds]

    # make sure to cover up all the gaps include bad areas
    if varframe:
        baddata = (outdata == 0)
        baddata = nd.maximum_filter(baddata, size=3)
        bpmdata[baddata] = 1
        

    # fill in the gaps if requested
    if fill:
        if varframe:
            outdata = fill_gaps(outdata, 0)
        else:
            outdata = fill_gaps(outdata, 0)

    # add to the file
    outlist[1] = fits.ImageHDU(outdata)
    if varframe:
        outlist[2] = fits.ImageHDU(vardata,name='VAR')
        outlist[3] = fits.ImageHDU(bpmdata,name='BPM')

    # create the image structure
    outstruct = fits.HDUList(outlist)

    # update the head informaation
    # housekeeping keywords
    saltkey.put('NEXTEND', 2, outstruct[0])
    saltkey.new('EXTNAME', 'SCI', 'Extension name', outstruct[1])
    saltkey.new('EXTVER', 1, 'Extension number', outstruct[1])
    if varframe:
        saltkey.new('VAREXT', 2, 'Variance frame extension', outstruct[1])
        saltkey.new('BPMEXT', 3, 'BPM Extension', outstruct[1])

    try:
        saltkey.copy(struct[1], outstruct[1], 'CCDSUM')
    except:
        pass

    # Add keywords associated with geometry
    saltkey.new('SGEOMGAP', gap, 'SALT Chip Gap', outstruct[0])
    c1str = '{:3.2f} {:3.2f} {:3.4f}'.format(xshift[0],
                                     yshift[0],
                                     rotation[0])
    saltkey.new('SGEOM1', c1str, 'SALT Chip 1 Transform', outstruct[0])
    c2str = '{:3.2f} {:3.2f} {:3.4f}'.format(xshift[1],
                                     yshift[1],
                                     rotation[1])
    saltkey.new('SGEOM2', c2str, 'SALT Chip 2 Transform', outstruct[0])

    # WCS keywords
    saltkey.new('CRPIX1', 0, 'WCS: X reference pixel', outstruct[1])
    saltkey.new('CRPIX2', 0, 'WCS: Y reference pixel', outstruct[1])
    saltkey.new(
        'CRVAL1',
        float(xbin),
        'WCS: X reference coordinate value',
        outstruct[1])
    saltkey.new(
        'CRVAL2',
        float(ybin),
        'WCS: Y reference coordinate value',
        outstruct[1])
    saltkey.new('CDELT1', float(xbin), 'WCS: X pixel size', outstruct[1])
    saltkey.new('CDELT2', float(ybin), 'WCS: Y pixel size', outstruct[1])
    saltkey.new('CTYPE1', 'pixel', 'X type', outstruct[1])
    saltkey.new('CTYPE2', 'pixel', 'Y type', outstruct[1])

    # cleanup temporary files
    if cleanup:
        for tfile in tranfile:
            if os.path.isfile(tfile):
                saltio.delete(tfile)
        if os.path.isfile(tilefile):
            status = saltio.delete(tilefile)

    # return the file
    return outstruct


def fill_gaps(data, mask):
    """Interpolate in the gaps in the data

       Parameters
       ----------
       data: np.ndarray
          data to have values filled in for

       mask: float or nd.ndarray
          If an nd.ndarray, it will be assumed to be a mask
          with values equal to 1 where they should be interpolated
          over.  If a float, pixels with that value will be replaced

    """
    ys, xs = data.shape
    if isinstance(mask, numpy.ndarray):
        mask = (mask == 0)
        for i in range(ys):
            x = numpy.arange(xs)
            rdata = data[i, :]
            rmask = mask[i, :]
            rmask = nd.minimum_filter(rmask, size=3)
            if rmask.any() == True:
                rdata = numpy.interp(x, x[rmask], rdata[rmask])
                data[i, rmask == 0] = rdata[rmask == 0]
    else:
        mask = (data != mask)
        for i in range(ys):
            x = numpy.arange(xs)
            rdata = data[i, :]
            rmask = mask[i, :]
            rmask = nd.minimum_filter(rmask, size=3)
            if rmask.any() == True:
                rdata = numpy.interp(x, x[rmask], rdata[rmask])
                data[i, rmask == 0] = rdata[rmask == 0]

    return data


def tran_func(a, xshift, yshift, xmag, ymag, xrot, yrot):
    xtran = ymag * a[0] * cos(yrot * pi / 180.0) \
        - xmag * a[1] * sin(xrot * pi / 180) \
        - yshift
    ytran = ymag * a[0] * sin(yrot * pi / 180.0) \
        + xmag * a[1] * cos(xrot * pi / 180) \
        - xshift
    return xtran, ytran


# -----------------------------------------------------------
# main code
if not iraf.deftask('saltmosaic'):
    parfile = iraf.osfn("saltred$saltmosaic.par")
    t = iraf.IrafTaskFactory(
        taskname="saltmosaic",
        value=parfile,
        function=saltmosaic,
        pkgname='saltred')
