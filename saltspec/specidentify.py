#!/usr/bin/env python
################################# LICENSE ##################################
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved. See LICENSE file for more details                   #
############################################################################
"""
SPECIDENTIFY  is a program to read in SALT RSS spectroscopic arc lamps and
determine the wavelength solution for that data.  The input data should be
a SALT arc lamp and a line list or another arc lamp image with high quality
wavelength solution. The lamp list can be either wavelengths, wavelengths
and fluxes, or an arc image with a high quality solution.  The line lamp
can also be left unspecified and the user will manually enter the data.

From there, the user has several different possible choices.  They can provide
a first guess of the coefficients for the wavelength solution or transformation,
indicate the model for the spectragraph for the first guess, or provide an
image with a solution already as the first guess.

Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)    1.0       10 Oct 2009

TODO
----


LIMITATIONS
-----------
1. Currently assumes that the linelist is of the form of an ascii file with
   either appropriate information in either one or two columns
2. Defaults to the central line

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
import WavelengthSolution


from PySpectrograph.Models import RSSModel

import spectools as st
import mostools as mt
from specrectify import readsolascii, findlinesol, enterdatetime

from spectools import SALTSpecError
from InterIdentify import InterIdentify
from AutoIdentify import AutoIdentify

debug = True


# -----------------------------------------------------------
# core routine

def specidentify(images, linelist, outfile, guesstype='rss', guessfile='',
                 automethod='Matchlines', function='poly', order=3, rstep=100,
                 rstart='middlerow', mdiff=5, thresh=3, niter=5, smooth=0,
                 subback=0, inter=True, startext=0,  clobber=False, 
                 textcolor='black', preprocess=False,  logfile='salt.log', verbose=True):

    with logging(logfile, debug) as log:

        # set up the variables
        infiles = []
        outfiles = []

        # Check the input images
        infiles = saltio.argunpack('Input', images)

        # create list of output files
        outfiles = saltio.argunpack('Output', outfile)

        # open the line lists
        slines, sfluxes = st.readlinelist(linelist)

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
                        xpos = -0.2666
                        ypos = 0.0117
                        objid = None
                    elif masktype == 'MOS':
                        slit = 1.5
                        #slit=saltkey.get('SLIT', hdu[i])

                        # set up the x and y positions
                        miny = hdu[i].header['MINY']
                        maxy = hdu[i].header['MAXY']
                        ras = hdu[i].header['SLIT_RA']
                        des = hdu[i].header['SLIT_DEC']
                        objid = hdu[i].header['SLITNAME']

                        # TODO: Check the perfomance of masks at different PA
                        rac = hdu[0].header['MASK_RA']
                        dec = hdu[0].header['MASK_DEC']
                        pac = hdu[0].header['PA']

                        # these are hard wired at the moment
                        xpixscale = 0.1267 * xbin
                        ypixscale = 0.1267 * ybin
                        cx = int(3162 / xbin)
                        cy = int(2050 / ybin)

                        x, y = mt.convert_fromsky(ras, des, rac, dec, xpixscale=xpixscale, ypixscale=ypixscale,
                                                  position_angle=-pac, ccd_cx=cx, ccd_cy=cy)
                        xpos = 0.015 * 2 * (cx - x[0])
                        ypos = 0.0117
                    else:
                        msg = '%s is not a currently supported masktype' % masktype
                        raise SALTSpecError(msg)

                    if instrume not in ['PFIS', 'RSS']:
                        msg = '%s is not a currently supported instrument' % instrume
                        raise SALTSpecError(msg)

                    # create RSS Model
                    rss = RSSModel.RSSModel(grating_name=grating.strip(), gratang=grang,
                                            camang=arang, slit=slit, xbin=xbin, ybin=ybin,
                                            xpos=xpos, ypos=ypos)
                    res = 1e7 * rss.calc_resolelement(rss.alpha(), -rss.beta())
                    dres = res / 10.0
                    wcen = 1e7 * rss.calc_centralwavelength()
                    R = rss.calc_resolution(
                        wcen / 1e7, rss.alpha(), -rss.beta())
                    logmsg = '\nGrating\tGR-ANGLE\tAR-ANGLE\tSlit\tWCEN\tR\n'
                    logmsg += '%s\t%8.3f\t%8.3f\t%4.2f\t%6.2f\t%4f\n' % (
                        grating, grang, arang, slit, wcen, R)
                    if log:
                        log.message(logmsg, with_header=False)

                    # set up the data for the source
                    try:
                        data = hdu[i].data
                    except Exception, e:
                        message = 'Unable to read in data array in %s because %s' % (
                            img, e)
                        raise SALTSpecError(message)

                    # set up the center row
                    if rstart == 'middlerow':
                        ystart = int(0.5 * len(data))
                    else:
                        ystart = int(rstart)

                    rss.gamma = 0.0
                    if masktype == 'MOS':
                        rss.gamma = 180.0 / math.pi * math.atan((y * rss.detector.pix_size * rss.detector.ybin
                                                                 - 0.5 * rss.detector.find_height()) / rss.camera.focallength)

                    # set up the xarr array based on the image
                    xarr = np.arange(len(data[ystart]), dtype='int64')

                    # get the guess for the wavelength solution
                    if guesstype == 'rss':
                        # set up the rss model
                        ws = st.useRSSModel(
                            xarr, rss, function=function, order=order, gamma=rss.gamma)
                        if function in ['legendre', 'chebyshev']:
                            ws.func.func.domain=[xarr.min(), xarr.max()]
                    elif guesstype == 'file':
                        soldict = {}
                        soldict = readsolascii(guessfile, soldict)
                        timeobs = enterdatetime('%s %s' % (dateobs, utctime))
                        exptime = saltkey.get('EXPTIME', hdu[0])
                        filtername = saltkey.get('FILTER', hdu[0]).strip()
                        try:
                            slitid = saltkey.get('SLITNAME', hdu[i])
                        except:
                            slitid = None

                        function, order, coef, domain = findlinesol(
                            soldict, ystart, True, timeobs, exptime, instrume, grating, grang, arang, filtername, slitid, xarr=xarr)
                        ws = WavelengthSolution.WavelengthSolution(
                            xarr, xarr, function=function, order=order)
                        ws.func.func.domain = domain
                        ws.set_coef(coef)
                    else:
                        raise SALTSpecError(
                            'This guesstype is not currently supported')


                    # identify the spectral lines
                    ImageSolution = identify(data, slines, sfluxes, xarr, ystart, ws=ws, function=function,
                                             order=order, rstep=rstep,  mdiff=mdiff, thresh=thresh, niter=niter,
                                             method=automethod, res=res, dres=dres, smooth=smooth, inter=inter, filename=img,
                                             subback=0, textcolor=textcolor, preprocess=preprocess, log=log, verbose=verbose)

                    if outfile and len(ImageSolution):
                        writeIS(ImageSolution, outfile, dateobs=dateobs, utctime=utctime, instrume=instrume,
                                grating=grating, grang=grang, grasteps=grasteps, arsteps=arsteps,
                                arang=arang, rfilter=rssfilter, slit=slit, xbin=xbin,
                                ybin=ybin, objid=objid, filename=img, log=log, verbose=verbose)


def identify(data, slines, sfluxes, xarr, istart, ws=None, function='poly', order=3,
             rstep=1, nrows=1, mdiff=5, thresh=3, niter=5, dc=3, ndstep=50, dsigma=5,
             method='Zeropoint', res=2, dres=0.2, filename=None, smooth=0, inter=True,
             subback=0, textcolor='green', preprocess=False, log=None, verbose=True):
    """For a given image, find the solution for each row in the file.  Use the appropriate first guess and
       guess type along with the appropriate function and order for the fit.

       Write out the new image with the solution in the headers and/or as a table in the multi-extension
       fits file

       data--2-D array of the arc image
       slines--wavelengths for known lines
       sfluxes--fluxes for known lines
       xarr--1-D array of x positions for data
       ws--WavelengthSolution expersion of the initial guess
       function--function form fit
       order--order of the fit
       rstep--rows to step when fitting the solution
       nrows--number of rows to average when extracting a line
       mdiff--limit in pixels for matching
       thresh--threshhold for detecting lines and for rejecting sources
       niter--maximum number of iterations to make
       dc--step size for zero point cross correlation
       dstep--number of
       method--method for automatic wavelength determination
       dsigma--detection threshhold for the lines
       rstart--first line to extract
       subback--order for function to subtract off background
       inter--run in interactive mode
       textcolor--color for wavelengths
       verbose--print out the results

       returns
    """
    ImageSolution = {}

    # update some of the parameters of the wavelength solution
    if ws is not None:
        ws.thresh = thresh
        ws.niter = niter

    # run in either interactive or non-interactive mode
    if inter:
        ImageSolution = InterIdentify(xarr, data, slines, sfluxes, ws, mdiff=mdiff, rstep=rstep,
                                      function=function, order=order, sigma=thresh, niter=niter,
                                      res=res, dres=dres, dc=dc, ndstep=ndstep, istart=istart,
                                      method=method, smooth=smooth, filename=filename,
                                      subback=subback, textcolor=textcolor, preprocess=preprocess, log=log, verbose=True)
    else:
        ImageSolution = AutoIdentify(xarr, data, slines, sfluxes, ws,
                                     rstep=rstep, method=method, istart=istart, nrows=nrows, mdiff=mdiff,
                                     dsigma=dsigma, res=res, dres=2 * dres, dc=dc, ndstep=ndstep, sigma=thresh,
                                     subback=subback, smooth=smooth, niter=niter, log=log, verbose=verbose)

    return ImageSolution


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
        if isinstance(ImageSolution[key_arr[i]], WavelengthSolution.WavelengthSolution):
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
    if ws.func.func.domain is not None: 
        msg += '#domain={},{}\n'.format(ws.func.func.domain[0], ws.func.func.domain[1])
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

if not iraf.deftask('specidentify'):
    parfile = iraf.osfn("saltspec$specidentify.par")
    t = iraf.IrafTaskFactory(
        taskname="specidentify", value=parfile, function=specidentify, pkgname='saltspec')
