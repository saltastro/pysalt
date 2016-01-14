#!/usr/bin/env python
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved. See LICENSE for more details                        #
"""
SPECWAVMAP is a program to read in SALT RSS spectroscopic data and
a wavelength solution for that data.  It will then create a wavemap, where
each pixel will have a corresponding wavelength

If a wavemap already exists and the solution for that data set only has a 
single line in it, then it updates the existing wave map based on that solution
and the values in the wavemap.

Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)    1.0       06 Jan 2015

TODO
----


LIMITATIONS
-----------

"""
# Ensure python 2.5 compatibility
from __future__ import with_statement

import os
import string
import sys
import glob
import math
import datetime
import numpy as np
from astropy.io import fits

from pyraf import iraf
import saltsafekey as saltkey
import saltsafeio
from saltsafelog import logging, SaltLog

from PySpectrograph.Models import RSSModel

import WavelengthSolution
import spectools as st
from spectools import SALTSpecError
import specrectify as sr


debug = True
import pylab as pl


# -------------------------------------------------------------------------
# core routine

def specwavemap(images, outimages, outpref, solfile=None, caltype='line',
                function='polynomial', order=3, blank=0, nearest=False,
                clobber=True, logfile='salt.log', verbose=True):

    with logging(logfile, debug) as log:

        # set up the variables
        infiles = []
        outfiles = []

        # Check the input images
        infiles = saltsafeio.argunpack('Input', images)

        # create list of output files
        outfiles = saltsafeio.listparse(
            'Outimages',
            outimages,
            outpref,
            infiles,
            '')

        # read in the wavelength solutions and enter them into
        # there own format
        if caltype == 'line':
            soldict = sr.entersolution(solfile)
        else:
            soldict = None

        # read in rectify each image
        for img, oimg, in zip(infiles, outfiles):
            if caltype == 'line':
                msg = 'Creating wave map image %s from image %s using files %s' % (
                    oimg, img, solfile)
            else:
                msg = 'Creating wave map image %s from image %s using RSS Model' % (
                    oimg, img)
            log.message(msg)
            hdu = saltsafeio.openfits(img)
            hdu = wavemap(hdu, soldict, caltype=caltype, function=function,
                          order=order,blank=blank, nearest=nearest, 
                          clobber=clobber, log=log, verbose=verbose)
            saltsafeio.writefits(hdu, oimg, clobber=clobber)


# -----------------------------------------------------------
# rectify data

def wavemap(hdu, soldict, caltype='line', function='poly', order=3, 
            blank=0, nearest=False, array_only=False, clobber=True, log=None,
            verbose=True):
    """Read in an image and a set of wavlength solutions.  Calculate the best
       wavelength solution for a given dataset and then apply that data set to the
       image

     return
    """

    # set up the time of the observation
    dateobs = saltkey.get('DATE-OBS', hdu[0])
    utctime = saltkey.get('TIME-OBS', hdu[0])
    exptime = saltkey.get('EXPTIME', hdu[0])
    instrume = saltkey.get('INSTRUME', hdu[0]).strip()
    grating = saltkey.get('GRATING', hdu[0]).strip()
    if caltype == 'line':
        grang = saltkey.get('GRTILT', hdu[0])
        arang = saltkey.get('CAMANG', hdu[0])
    else:
        grang = saltkey.get('GR-ANGLE', hdu[0])
        arang = saltkey.get('AR-ANGLE', hdu[0])
    filtername = saltkey.get('FILTER', hdu[0]).strip()
    slitname = saltkey.get('MASKID', hdu[0])
    slit = st.getslitsize(slitname)
    xbin, ybin = saltkey.ccdbin(hdu[0])

    timeobs = sr.enterdatetime('%s %s' % (dateobs, utctime))

    # check to see if there is more than one solution
    if caltype == 'line':
        if len(soldict) == 1:
            sol = soldict.keys()[0]
            slitid = None
            if not sr.matchobservations(
                    soldict[sol], instrume, grating, grang, arang, filtername, slitid):
                msg = 'Observations do not match setup for transformation but using the solution anyway'
                if log:
                    log.warning(msg)

    for i in range(1, len(hdu)):
        if hdu[i].name == 'SCI':
            if log:
                log.message('Correcting extension %i' % i)
            istart = int(0.5 * len(hdu[i].data))

            # open up the data
            # set up the xarr and initial wavlength solution
            xarr = np.arange(len(hdu[i].data[istart]), dtype='int64')

            # get the slitid
            try:
                slitid = saltkey.get('SLITNAME', hdu[i])
            except:
                slitid = None
           
            #check to see if wavext is already there and if so, then check update
            #that for the transformation from xshift to wavelength
            if saltkey.found('WAVEXT',hdu[i]):
                w_ext = saltkey.get('WAVEXT', hdu[i])-1
                wavemap = hdu[w_ext].data
                function, order, coef = sr.findlinesol(soldict, istart, nearest, 
                          timeobs, exptime, instrume, grating, grang, arang, 
                          filtername, slitid, xarr)
                ws = WavelengthSolution.WavelengthSolution(
                     xarr,
                     xarr,
                     function=function,
                     order=order)
                ws.set_coef(coef)
                for j in range(len(hdu[i].data)):
                    wavemap[j,:] = ws.value(wavemap[j,:])
                if array_only: return wavemap
                hdu[w_ext].data = wavemap
                continue
 
            # set up a wavelength solution -- still in here for testing MOS data
            try:
                w_arr = sr.findsol(xarr, soldict, istart, caltype, nearest, timeobs, exptime, instrume, grating, grang, arang, filtername,
                                slit, xbin, ybin, slitid, function, order)
            except SALTSpecError as e:
                if slitid:
                    msg = 'SLITID %s: %s' % (slitid, e)
                    if log:
                        log.warning(msg)
                    continue
                else:
                    raise SALTSpecError(e)

            if w_arr is None:
                w_arr = sr.findsol(xarr, soldict, istart, 'rss', nearest, timeobs, exptime, instrume, grating, grang, arang, filtername,
                                slit, xbin, ybin, slitid, function, order)

            # for each line in the data, determine the wavelength solution
            # for a given line in the image
            wavemap = np.zeros_like(hdu[i].data)
            for j in range(len(hdu[i].data)):
                # find the wavelength solution for the data
                w_arr = sr.findsol(xarr, soldict, j, caltype, nearest, timeobs, exptime, instrume, grating, grang, arang, filtername,
                                slit, xbin, ybin, slitid, function, order)
                if w_arr is not None: wavemap[j,:] = w_arr
            if array_only: return wavemap

            # write out the oimg
            hduwav = fits.ImageHDU(data=wavemap, header=hdu[i].header, name='WAV')
            hdu.append(hduwav)
            saltkey.new('WAVEXT', len(hdu)-1, 'Extension for Wavelength Map', hdu[i])

    return hdu

if not iraf.deftask('specwavemap'):
    parfile = iraf.osfn("saltspec$specwavemap.par")
    t = iraf.IrafTaskFactory(taskname="specwavemap", value=parfile,
                             function=specwavemap, pkgname='saltspec')

