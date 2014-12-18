#!/usr/bin/env python
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.                                                     #

"""
specreduce is a program to reduce a set of RSS spectroscopic data.
The task will calculate the wavelength calibration for any
arc images, rectify all science images, sky subtract the data,
and produce 1-D spectra assuming the data is longslit.

There are three main methods for reduction which are None, rss,
and line.  If None, it will do no wavelength calibration.  If the
method is rss, then it will only use the RSS model to calculate
and determine the wavelength calibrate.   If it is line, then it
will go through and create the wavelength solution first,
and then correct for that.

Author                 Version      Date
-----------------------------------------------
 S. M. Crawford (SAAO)    0.1       1 April 2011

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
from specprepare import prepare
from specrectify import rectify
from specsky import skysubtract
from specextract import extract, write_extract
from spectools import makesection

debug = True

# -----------------------------------------------------------
# core routine


def specreduce(images, badpixelimage=None, caltype='rss',
               function='polynomial', order=3,
               skysub=True, skysection=None, findobj=False,
               objsection=None, thresh=3.0,
               clobber=True, logfile='salt.log', verbose=True):

    with logging(logfile, debug) as log:

        # Check the input images
        infiles = saltio.argunpack('Input', images)

        # open the badpixelstruct
        if saltio.checkfornone(badpixelimage):
            badpixelstruct = saltio.openfits(badpixelimage)
        else:
            badpixelstruct = None

        # set up the section for sky estimate
        if skysection is not None:
            skysection = makesection(skysection)
        else:
            skysub = False

        # set up the section for sky estimate
        section = saltio.checkfornone(objsection)
        if section is not None:
            sections = saltio.getSection(section, iraf_format=False)
            objsection = []
            for i in range(0, len(sections), 2):
                objsection.append((sections[i], sections[i + 1]))

        # determine the wavelength solutions
        if caltype == 'line':
            calc_wavesol(infiles)

        # correct the images
        for img in infiles:
            # open the fits file
            struct = saltio.openfits(img)

            # prepare filep
            log.message('Preparing %s' % img)
            struct = prepare(struct, badpixelstruct)

            # rectify the spectrum
            log.message('Rectifying %s using %s' % (img, caltype))
            struct = rectify(
                struct,
                None,
                caltype=caltype,
                function=function,
                order=order)

            # sky subtract the spectrum
            # assumes the data is long slit and in the middle of the field
            if skysub:
                log.message('Subtracting the sky from %s' % img)

                struct = skysubtract(
                    struct,
                    method='normal',
                    section=skysection)

            # extract the spectrum
            log.message('Extracting the spectrum from %s' % img)
            print objsection
            aplist = extract(
                struct,
                method='normal',
                section=objsection,
                thresh=thresh)
            oimg = os.path.dirname(
                os.path.abspath(img)) + '/s' + os.path.basename(img.strip())
            ofile = oimg[:-5] + '.txt'
            write_extract(ofile, aplist, clobber=clobber)

            # write FITS file
            log.message('Writing 2-D corrected image as %s' % oimg)
            saltio.writefits(struct, oimg, clobber=clobber)
            saltio.closefits(struct)


def calc_wavesol(infiles):
    """For each file in infiles, calculate whether the file
       is an arc frame and if so, calculate the wavlength
       solution
    """
    return

# -----------------------------------------------------------
# main code

parfile = iraf.osfn("saltspec$specreduce.par")
t = iraf.IrafTaskFactory(
    taskname="specreduce",
    value=parfile,
    function=specreduce,
    pkgname='saltspec')
