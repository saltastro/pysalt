#!/usr/bin/env python
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.                                                     #
"""
specprepare is a program to read in SALT RSS spectroscopic data
and prepare the data for any processing.  The main steps right
now are to check that the data have been reduced and to add
a variance from to the data

# Author                 Version      Date
# -----------------------------------------------
# S. M. Crawford (SAAO)    0.1       20 Jun 2009

"""

from __future__ import with_statement


from pyraf import iraf
from pyraf.iraf import pysalt
import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging
import os
import string
import sys
import glob
import pyfits
import time
import numpy

from spectools import SALTSpecError

debug = True

# -----------------------------------------------------------
# core routine


def specprepare(images, outimages, outpref, badpixelimage='',
                clobber=True, logfile='salt.log', verbose=True):

    with logging(logfile, debug) as log:

        # Check the input images
        infiles = saltio.argunpack('Input', images)

        # create list of output files
        outfiles = saltio.listparse('Outfile', outimages, outpref, infiles, '')

        # verify that the input and output lists are the same length
        saltio.comparelists(infiles, outfiles, 'Input', 'output')

        # open the badpixel image
        if saltio.checkfornone(badpixelimage) is None:
            badpixelstruct = None
        else:
            try:
                badpixelstruct = saltio.openfits(badpixelimage)
            except saltio.SaltIOError as e:
                msg = 'badpixel image must be specificied\n %s' % e
                raise SALTSpecError(msg)

        # open each raw image file

        for img, oimg, in zip(infiles, outfiles):

            # open the fits file
            struct = saltio.openfits(img)

            # prepare file
            struct = prepare(struct, badpixelstruct)

            # write FITS file
            saltio.writefits(struct, oimg, clobber=clobber)
            saltio.closefits(struct)

            message = 'SPECPREPARE -- %s => %s' % (img, oimg)
            log.message(message)


# -----------------------------------------------------------
# prepare FITS file for processing

def prepare(struct, badpixelstruct):
    """Prepare a structure for spectroscopic processing

    """

    # identify instrument
    infile = struct._HDUList__file.name

    nextend = len(struct)

    if badpixelstruct is not None:
        bpfile = badpixelstruct._HDUList__file.name
        # Check that the image is the same number of extensions as the
        # badpixelmap
        if len(struct) != len(badpixelstruct):
            message = '%s and %s are not the same length' % (infile, bpfile)
            raise SALTSpecError(message)

    # Add variance frames and bad pixel maps
    # If no variance frame exists, add the file at the end of the
    # extensions.  If a variance file does exist, do nothing
    # Add a keyword to the Science extension indication what
    # extension the variance from is on

    j = 0
    for i in range(nextend):
        if struct[i].size() > 0 and struct[i].name == 'SCI':
            # Check to see if the variance frame exists
            try:
                key = struct[i].header['VAREXT']
            except:
                try:
                    hdu = createvariance(struct[i], i, nextend + j)
                    struct[i].header.update(
                        'VAREXT',
                        nextend + j,
                        comment='Extension for Variance Frame')
                    struct.append(hdu)
                    j = j + 1
                except Exception as e:
                    message = 'Could not create variance extension in ext %i of %s because %s' \
                        % (i, infile, e)
                    raise SALTSpecError(message)

            # Check to see if the BPM  frame exists
            try:
                key = struct[1].header['BPMEXT']
            except:
                try:
                    hdu = createbadpixel(
                        struct,
                        badpixelstruct,
                        i,
                        nextend +
                        j)
                except Exception as e:
                    message = 'Could not create bad pixel extension in ext %i of %s because %s' \
                        % (i, infile, e)
                    raise SALTSpecError(message)
                if (1):
                    struct[i].header.update(
                        'BPMEXT',
                        nextend + j,
                        comment='Extension for Bad Pixel Mask')
                    struct.append(hdu)
                    j = j + 1

    return struct


def createbadpixel(inhdu, bphdu, sci_ext, bp_ext):
    """Create the bad pixel hdu from a bad pixel hdu"""
    if bphdu is None:
        data = inhdu[sci_ext].data * 0.0
    else:
        infile = inhdu._HDUList__file.name
        bpfile = bphdu._HDUList__file.name
        for k in ['INSTRUME', 'CCDSUM', 'NAXIS1', 'NAXIS2']:
            if not saltkey.compare(
                    inhdu[sci_ext], bphdu[sci_ext], k, infile, bpfile):
                message = '%s and %s are not the same %s' % (infile, bpfile, k)
                raise SALTSpecError(message)
        data = bphdu[sci_ext].data.copy()

    header = inhdu[sci_ext].header.copy()
    header['EXTVER'] = bp_ext
    header.update('SCIEXT', sci_ext, comment='Extension of science frame')
    return pyfits.ImageHDU(data=data, header=header, name='BPM')


def createvariance(inhdu, sci_ext, var_ext):
    """Create a variance hdu from an input hdu"""

    # create the variance array
    data = inhdu.data.copy()
    if (data <= 0).any():
        j = numpy.where(data > 0)
        min_pos = data[j].min()
        j = numpy.where(data <= 0)
        data[j] = min_pos
    data = data ** 0.5

    header = inhdu.header.copy()
    header['EXTVER'] = var_ext
    header.update('SCIEXT', sci_ext, comment='Extension of science frame')
    return pyfits.ImageHDU(data=data, header=header, name='VAR')

# -----------------------------------------------------------
# main code

parfile = iraf.osfn("saltspec$specprepare.par")
t = iraf.IrafTaskFactory(
    taskname="specprepare",
    value=parfile,
    function=specprepare,
    pkgname='saltspec')
