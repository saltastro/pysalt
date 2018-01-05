#!/usr/bin/env python
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.                                                     #
"""
SPECARCSTRAIGHT is a program to read in an arc lamp and cross-correlate
it with itself to straighten the lines.   This will not wavelength
calibrate the data but it will determine the correction for spectral
curvature along the spatial dimension.

Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)    1.0       12 Feb 2012

TODO
----


LIMITATIONS
-----------

"""
# Ensure python 2.5 compatibility
from __future__ import with_statement

import os
import sys
import math
import time
import numpy as np

from pyraf import iraf
import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging
from salterror import SaltError, SaltIOError
from scipy import ndimage as nd

from PySpectrograph.WavelengthSolution import WavelengthSolution
from PySpectrograph.Spectra import apext

import spectools as st
import mostools as mt

from spectools import SALTSpecError
from AutoIdentify import getwsfromIS


debug = True


# -----------------------------------------------------------
# core routine

def specarcstraighten(images, outfile, function='poly', order=3, rstep=1,
                      rstart='middlerow', nrows=1, 
                      y1=None, y2=None, sigma=5, sections=3, niter=5,
                      startext=0, clobber=False, logfile='salt.log', verbose=True):

    with logging(logfile, debug) as log:

        # Check the input images
        infiles = saltio.argunpack('Input', images)

        # create list of output files
        outfiles = saltio.argunpack('Output', outfile)

        # Identify the lines in each file
        for img, ofile in zip(infiles, outfiles):

            # open the image
            hdu = saltio.openfits(img)

            # get the basic information about the spectrograph
            dateobs = saltkey.get('DATE-OBS', hdu[0])
            try:
                utctime = saltkey.get('UTC-OBS', hdu[0])
            except SaltError:
                utctime = saltkey.get('TIME-OBS', hdu[0])

            instrume = saltkey.get('INSTRUME', hdu[0]).strip()
            grating = saltkey.get('GRATING', hdu[0]).strip()
            grang = saltkey.get('GR-ANGLE', hdu[0])
            grasteps = saltkey.get('GRTILT', hdu[0])
            arang = saltkey.get('AR-ANGLE', hdu[0])
            arsteps = saltkey.get('CAMANG', hdu[0])
            rssfilter = saltkey.get('FILTER', hdu[0])
            specmode = saltkey.get('OBSMODE', hdu[0])
            masktype = saltkey.get('MASKTYP', hdu[0]).strip().upper()
            slitname = saltkey.get('MASKID', hdu[0])
            xbin, ybin = saltkey.ccdbin(hdu[0], img)

            for i in range(startext, len(hdu)):
                if hdu[i].name == 'SCI':
                    log.message('Proccessing extension %i in  %s' % (i, img))
                    # things that will change for each slit

                    if masktype == 'LONGSLIT':
                        slit = st.getslitsize(slitname)
                        objid = None
                    #elif masktype == 'MOS':
                        #slit = 1.5
                        # slit=saltkey.get('SLIT', hdu[i])

                        # set up the x and y positions
                        #miny = hdu[i].header['MINY']
                        #maxy = hdu[i].header['MAXY']
                        #ras = hdu[i].header['SLIT_RA']
                        #des = hdu[i].header['SLIT_DEC']
                        #objid = hdu[i].header['SLITNAME']

                        # Check the perfomance of masks at different PA
                        #rac = hdu[0].header['MASK_RA']
                        #dec = hdu[0].header['MASK_DEC']
                        #pac = hdu[0].header['PA']

                    else:
                        msg = '%s is not a currently supported masktype' % masktype
                        raise SALTSpecError(msg)

                    if instrume not in ['PFIS', 'RSS']:
                        msg = '%s is not a currently supported instrument' % instrume
                        raise SALTSpecError(msg)

                    # set up the data for the source
                    try:
                        data = hdu[i].data
                    except Exception as e:
                        message = 'Unable to read in data array in %s because %s' % (
                            img, e)
                        raise SALTSpecError(message)

                    # set up the center row
                    if rstart == 'middlerow':
                        ystart = int(0.5 * len(data))
                    else:
                        ystart = rstart


                    # set up the xarr array based on the image
                    xarr = np.arange(len(data[ystart]), dtype='int64')

                    # calculate the transformation
                    ImageSolution = arcstraight(data, xarr, ystart, function=function, order=order, 
                                                rstep=rstep, y1=y1, y2=y2, sigma=sigma, sections=sections,
                                                niter=niter, log=log, verbose=verbose)

                    if outfile and len(ImageSolution):
                        writeIS(ImageSolution, outfile, dateobs=dateobs, utctime=utctime, instrume=instrume,
                                grating=grating, grang=grang, grasteps=grasteps, arsteps=arsteps,
                                arang=arang, rfilter=rssfilter, slit=slit, xbin=xbin,
                                ybin=ybin, objid=objid,
                                filename=img, log=log, verbose=verbose)


def arcstraight(data, xarr, istart, function='poly', order=3,
                rstep=1, y1=None, y2=None, sigma=5, sections=3, niter=5, 
                log=None, verbose=True):
    """For a given image, assume that the line given by istart is the fiducial and then calculate
       the transformation between each line and that line in order to straighten the arc

       returns Wavlenght solution
    """

    #set up the edges 
    if y1 is None: y1=0
    if y2 is None: y2=data.shape[0]

    # extract the central row
    ys, xs = data.shape
    farr = data[istart,:]
    xp, xf = st.findpoints(xarr, farr, sigma, niter, sections=sections)

    #short function to return closest peak
    def get_closest_x(y_arr, i):
        y = np.arange(len(y_arr))
        j = (abs(y[y_arr >0] - i)).argmin()
        return y_arr[y[y_arr >0][j]]

    #loop through all lines and rows to trace the lines
    line_dict = {}
    yc = int(ys/2.0)
    xdiff = 5
    for x in xp:
        y_arr = np.zeros(ys)
        y_arr[yc] = x
        for i in range(1, yc):
            if yc+i < y2:
                y_arr[yc+i] = st.mcentroid(xarr, data[yc+i,:], kern=st.default_kernal, 
                                       xdiff=xdiff, xc=get_closest_x(y_arr, yc+i))
            if yc+i > y1:
                y_arr[yc-i] = st.mcentroid(xarr, data[yc-i,:], kern=st.default_kernal, 
                                       xdiff=xdiff, xc=get_closest_x(y_arr, yc-i))
        line_dict[x] = y_arr

    #find the solution for each row
    iws={}
    xp = np.array(line_dict.keys())
    for i in range(y1, y2):  
        wp = [line_dict[x][i] for x in xp]
        ws = WavelengthSolution.WavelengthSolution(np.array(wp), xp, order=3, function='polynomial')
        ws.fit()
        exit()
        iws[i] = ws

    return iws


def writeIS(ImageSolution, outfile, dateobs=None, utctime=None, instrume=None,
            grating=None, grang=0.0, grasteps=None, objid=None,
            arang=0.0, arsteps=None, rfilter=None, slit=None, xbin=2, ybin=2,
            filename=None, log=None, verbose=False):

    # set up the list of solutions to into an array
    key_arr = np.array(ImageSolution.keys())
    arg_arr = key_arr.argsort()

    # set up the wavelength solution
    ws = ImageSolution[key_arr[0]]
    ws_arr = np.zeros((len(arg_arr), len(ws.coef) + 1), dtype=float)

    # write the solution to an array
    for j, i in enumerate(arg_arr):
        if isinstance(ImageSolution[key_arr[i]],
                      WavelengthSolution.WavelengthSolution):
            function = ImageSolution[key_arr[i]].function
            order = ImageSolution[key_arr[i]].order
            ws_arr[j, 0] = key_arr[i]
            ws_arr[j, 1:] = ImageSolution[key_arr[i]].coef

    # write header to the file that should include the order and function
    if os.path.isfile(outfile):
        dout = open(outfile, 'a')
    else:
        dout = open(outfile, 'w')

    msg = '#WS: Wavelength solution for image %s\n' % filename
    msg += '#The following parameters were used in determining the solution:\n'
    msg += '#name=%s\n' % filename
    msg += '#time-obs=%s %s\n' % (dateobs, utctime)
    msg += '#instrument=%s\n' % instrume
    msg += '#grating=%s\n' % grating.strip()
    msg += '#graang=%s\n' % grang
    msg += '#gratilt=%s\n' % grasteps
    msg += '#arang=%s\n' % arang
    msg += '#artilt=%s\n' % arsteps
    msg += '#filter=%s\n' % rfilter.strip()
    if objid:
        msg += '#slitid=%s\n' % objid
    msg += '#Function=%s\n' % function
    msg += '#Order=%s\n' % order
    msg += '#Starting Data\n'
    dout.write(msg)

    for i in range(len(ws_arr)):
        if ws_arr[i, 0]:
            msg = '%5.2f ' % ws_arr[i, 0]
            msg += ' '.join(['%e' % k for k in ws_arr[i, 1:]])
            dout.write(msg + '\n')
    dout.write('\n')
    dout.close()

    return


# main code

parfile = iraf.osfn("saltspec$specarcstraighten.par")
t = iraf.IrafTaskFactory(
    taskname="specarcstraighten",
    value=parfile,
    function=specarcstraighten,
    pkgname='saltspec')
