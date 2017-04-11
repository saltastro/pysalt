#!/usr/bin/env python
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.                                                     #
"""
SPECSLIT is a program that reads in MOS data and uses an edge detection
algorithm to determine the position of slits from a flat field. The slit
positions are then used to extract the spectra from the object frames into
individual spectra which are stored in a MEF.

Author                 Version      Date
-----------------------------------------------
J. P. Kotze (SAAO)     0.2       21 Sept 2010
S. M. Crawford (SAAO)  0.3       17 Apr  2012

TODO
----
- add a checker to verify the fits header information. this will probably mean
adding a function to the saltio lib.
- add the slit position writing to the output fits file in fits table format
and also the file output output.



LIMITATIONS
-----------
1. Does not take any input, everything is calculated from the input image
2. Sampling of the spectrum is rough, therefore the spline fitting is not
very accurate
"""

from __future__ import with_statement

import os
import sys
import numpy as np

from astropy.io import fits

from pyraf import iraf

import saltsafeio as saltio
import saltsafekey as saltkey
from saltsafelog import logging

import mostools as mt
from spectools import SALTSpecError


debug = True


def specslit(image, outimage, outpref, exttype='auto', slitfile='', outputslitfile='',
             regprefix='', sections=3, width=25, sigma=2.2, thres=6, order=3, padding=5, yoffset=0,
             inter=False, clobber=True, logfile='salt.log', verbose=True):

    with logging(logfile, debug) as log:

        # check all the input and make sure that all the input needed is provided
        # by the user

        # read the image or image list and check if each in the list exist
        infiles = saltio.argunpack('Input', image)

        # unpack the outfiles
        outfiles = saltio.listparse(
            'Outimages',
            outimage,
            outpref,
            infiles,
            '')

        # from the extraction type, check whether the input file is specified.
        # if the slitfile parameter is specified then use the slit files for
        # the extraction. if the extraction type is auto then use image for the
        # detection and the slit extraction

        if exttype == 'rsmt' or exttype == 'fits' or exttype == 'ascii' or exttype == 'ds9':
            slitfiles = saltio.argunpack('Slitfile', slitfile)
            if len(slitfiles) == 1:
                slitfiles = slitfiles * len(infiles)
            saltio.comparelists(infiles, slitfiles, 'image', 'slitfile')
        elif exttype == 'auto':
            slitfiles = infiles
            log.message(
                'Extraction type is AUTO. Slit detection will be done from image')

        # read in if an optional ascii file is requested
        if len(outputslitfile) > 0:
            outslitfiles = saltio.argunpack('Outslitfiles', outputslitfile)
            saltio.comparelists(
                infiles,
                outslitfiles,
                'image',
                'outputslitfile')
        else:
            outslitfiles = [''] * len(infiles)

        # check if the width and sigma parameters were specified.
        # default is 25 and 2.2
        if width < 10.:
            msg = 'The width parameter needs be a value larger than 10'
            raise SALTSpecError(msg)

        if sigma < 0.0:
            msg = 'Sigma must be greater than zero'
            raise SaltSpecError(msg)

        # check the treshold parameter. this needs to be specified by the user
        if thres <= 0.0:
            msg = 'Threshold must be greater than zero'
            raise SaltSpecError(msg)

        # check to make sure that the sections are greater than the order
        if sections <= order:
            msg = 'Number of sections must be greater than the order for the spline fit'
            raise SaltSpecError(msg)

        # run through each of the images and extract the slits
        for img, oimg, sfile, oslit in zip(
                infiles, outfiles, slitfiles, outslitfiles):
            log.message('Proccessing image %s' % img)

            # open the image
            struct = saltio.openfits(img)
            ylen, xlen = struct[1].data.shape
            xbin, ybin = saltkey.ccdbin(struct[0], img)
            # setup the VARIANCE and BPM frames
            if saltkey.found('VAREXT', struct[1]):
                varext = saltkey.get('VAREXT', struct[1])
                varlist = []
            else:
                varext = None

            # setup the BPM frames
            if saltkey.found('BPMEXT', struct[1]):
                bpmext = saltkey.get('BPMEXT', struct[1])
                bpmlist = []
            else:
                bpmext = None

            # open the slit definition file or identify the slits in the image
            slitmask = None
            ycheck = False
            if exttype == 'rsmt':
                log.message('Using slits from %s' % sfile)
                if yoffset is None:
                   yoffset = 0 
                   ycheck = True
                slitmask = mt.read_slitmask_from_xml(sfile)
                xpos = -0.3066
                ypos = 0.0117
                cx = int(xlen / 2.0)
                cy = int(ylen / 2.0) + ypos / 0.015 / ybin + yoffset
                order, slit_positions = mt.convert_slits_from_mask(
                    slitmask, order=1, xbin=xbin, ybin=ybin, pix_scale=0.1267, cx=cx, cy=cy)
                sections = 1
            elif exttype == 'fits':
                log.message('Using slits from %s' % sfile)
                order, slit_positions = read_slits_from_fits(sfile)
            elif exttype == 'ascii':
                log.message('Using slits from %s' % sfile)
                order, slit_positions = mt.read_slits_from_ascii(sfile)
            elif exttype == 'ds9':
                log.message('Using slits from %s' % sfile)
                order, slit_positions, slitmask = mt.read_slits_from_ds9(
                    sfile, order=order)
                slitmask = None
                sections = 1
            elif exttype == 'auto':
                log.message('Identifying slits in %s' % img)
                # identify the slits in the image
                order, slit_positions = identify_slits(
                    struct[1].data, order, sections, width, sigma, thres)

                # write out the slit identifications if ofile has been supplied
                if oslit:
                    log.message('Writing slit positions to %s' % oslit)
                    mt.write_outputslitfile(slit_positions, oslit, order)

            if ycheck:
               slit_positions, dy = check_ypos(slit_positions, struct[1].data) 
               log.message('Using an offset of {}'.format(dy))

            # extract the slits
            spline_x = mt.divide_image(struct[1].data, sections)
            spline_x = 0.5 * (np.array(spline_x[:-1]) + np.array(spline_x[1:]))
            extracted_spectra, spline_positions = mt.extract_slits(slit_positions,
                                                                   spline_x, struct[1].data, order=order, padding=padding)
            if varext:
                extracted_var, var_positions = mt.extract_slits(slit_positions,
                                                                spline_x, struct[varext].data, order=order, padding=padding)
            if bpmext:
                extracted_bpm, bpm_positions = mt.extract_slits(slit_positions,
                                                                spline_x, struct[bpmext].data, order=order, padding=padding)

            # write out the data to the new array
            # create the new file
            hdulist = fits.HDUList([struct[0]])

            # log the extracted spectra if needed
            log.message('', with_stdout=verbose)

            # setup output ds9 file
            if regprefix:
                regout = open(
                    regprefix +
                    os.path.basename(img).strip('.fits') +
                    '.reg',
                    'w')
                regout.write('# Region file format: DS9 version 4.1\n')
                regout.write('# Filename: %s\n' % img)
                regout.write(
                    'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nphysical\n')

            # add each
            imglist = []
            nslits = len(spline_positions)
            for i in range(nslits):
                y1 = spline_positions[i][0].min()
                y2 = spline_positions[i][1].max()
                msg = 'Extracted Spectra %i between %i to %i' % (i + 1, y1, y2)
                # log.message(msg, with_header=False, with_stdout=verbose)
                sdu = fits.ImageHDU(
                    extracted_spectra[i],
                    header=struct[1].header)
                if varext:
                    vdu = fits.ImageHDU(
                        extracted_var[i],
                        header=struct[varext].header)
                    sdu.header['VAREXT'] = i + nslits + 1
                    varlist.append(vdu)
                if bpmext:
                    bdu = fits.ImageHDU(
                        extracted_bpm[i],
                        header=struct[bpmext].header)
                    sdu.header['BPMEXT']= i + 2 * nslits + 1
                    bpmlist.append(bdu)
                imglist.append(sdu)

                # add in some additional keywords
                imglist[i].header['MINY'] = (y1,
                    'Lower Y value in original image')
                imglist[i].header['MAXY'] = (y2,
                    'Upper Y value in original image')
                if regprefix:
                    xsize = struct[1].data.shape[1]
                    xsize = int(0.5 * xsize)
                    rtext = ''
                    if slitmask:
                        # rtext='%s, %8.7f, %8.7f, %3.2f' % (slitmask.slitlets.data[i]['name'], slitmask.slitlets.data[i]['targ_ra'], slitmask.slitlets.data[i]['targ_dec'], slitmask.slitlets.data[i]['slit_width'])
                        pass
                    regout.write('box(%i,%i, %i, %i) #text={%s}\n' % (
                        xsize, 0.5 * (y1 + y2), 2 * xsize, y2 - y1, rtext))

                # add slit information
                if slitmask:
                    imglist[i].header['SLITNAME'] = (slitmask.slitlets.data[i]['name'],
                        'Slit Name')
                    imglist[i].header['SLIT_RA'] = (slitmask.slitlets.data[i]['targ_ra'],
                        'Slit RA')
                    imglist[i].header['SLIT_DEC'] = (slitmask.slitlets.data[i]['targ_dec'],
                        'Slit DEC')
                    imglist[i].header['SLIT'] = (slitmask.slitlets.data[i]['slit_width'],
                        'Slit Width')

            # add to the hdulist
            hdulist += imglist
            if varext:
                hdulist += varlist
            if bpmext:
                hdulist += bpmlist

            # write the slit positions to the header
            # create the binary table HDU that contains the split positions
            tbhdu = mt.slits_HDUtable(slit_positions, order)
            bintable_hdr = tbhdu.header

            # add the extname parameter to the extension
            tbhdu.header['EXTNAME'] = 'BINTABLE'

            # add the extname parameter to the extension
            hdulist[0].header['SLITEXT'] = len(hdulist)
            hdulist.append(tbhdu)

            # add addition header information about the mask
            if slitmask:
                hdulist[0].header['MASKNAME'] = (slitmask.mask_name, 'SlitMask Name')
                hdulist[0].header['MASK_RA'] = (slitmask.center_ra,
                    'SlitMask RA')
                hdulist[0].header['MASK_DEC'] = ( slitmask.center_dec,
                    'SlitMask DEC')
                hdulist[0].header['MASK_PA'] = ( slitmask.position_angle,
                    'SlitMask Position Angle')

            # write out the image
            saltio.writefits(hdulist, oimg, clobber)


def identify_slits(slit_img, order=3, section=3, width=25, sigma=2.5, thres=6):
    """determine the x and y positions of the slits from identification in the image"""

    # split the image into sections determined from sect
    inc = mt.divide_image(slit_img, section)

    # determine where the slits are located, if option = True then only return the
    # edges, if option = False then return all the additional info.
    slits = mt.locate_slits(slit_img, inc, width, sigma, thres, False)

    # Check the detected slits for any variations in the amount of slits detected
    # for each image section and that the same amount of edges for left and right
    # have been detected for each image section.
    slit_positions = mt.get_slits(slits)

    return order, slit_positions


def read_slits_from_fits(simg):
    """Read the slit definitions from a FITS fiel where the slit
       definitions have been stored in a table in the FITS file under an
       extension with the name
    """

    # first check if the file exists
    saltio.fileexists(simg)

    # open the slit image
    struct = saltio.openfits(simg)

    # get the extension of the slit table
    slitext = saltkey.get('SLITEXT', struct[0])

    # extract the order of the fit and the positions of each of the slits
    order, slit_positions = mt.read_slits_HDUtable(struct[slitext])

    return order, slit_positions


def check_ypos(slit_positions, data):
    """Check the y-position"""
    from scipy.signal import find_peaks_cwt

    # determine the points that 
    y = data.sum(axis=1)
    yp = find_peaks_cwt(np.gradient(y), np.array([5]))
    x = np.arange(len(y))
    s_min= x.max()
    s_max= x.min()
    for n, ym, yx  in slit_positions: 
        if ym < s_min: s_min = ym
        if yx > s_max: s_max = yx
    dy = yp[yp>10].min()-s_min
    for i in range(len(slit_positions)):
        n, ym, yx = slit_positions[i]
        slit_positions[i] = [n, ym+dy, yx+dy]
    
    return slit_positions, dy

# main code
parfile = iraf.osfn("saltspec$specslit.par")
t = iraf.IrafTaskFactory(taskname="specslit", value=parfile, function=specslit,
                         pkgname='saltspec')
