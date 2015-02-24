#!/usr/bin/env python
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.   See LICENSE for more details			   #
"""
AutoIDENTIFY  is a program to automatically identify spectral lines in
an arc image.

Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)    1.0       21 Aug 2010

TODO
----


LIMITATIONS
-----------

"""
# Ensure python 2.5 compatibility

import os
import sys
import time
import numpy as np

from pyraf import iraf
import saltprint
import saltio
import saltkey
import saltsafekey
import saltsafeio
from saltsafelog import logging
from salterror import SaltError, SaltIOError


from PySpectrograph.Spectra import apext


import spectools as st
from spectools import SALTSpecError

debug = True

# autoidentify_options=['Zeropoint', 'Matchlines', 'MatchZero', 'FullXCor']
autoidentify_options = ['Zeropoint', 'Matchlines', 'MatchZero']


def AutoIdentify(xarr, specarr, slines, sfluxes, ws, method='Zeropoint',
                 rstep=1, istart=None, nrows=1, res=2, dres=0.1, dsigma=5,
                 sigma=5, smooth=0, niter=5, mdiff=20, dc=20, ndstep=20, farr=None,
                 subback=0, oneline=False, log=None, verbose=True):
    """Automatically find the wavlength solution for the entire image.  The following
       methods are used:

       Zeropoint--Assume that the form for the initial guess of the wavelength solution
                  is correct

       Matchlines--For each line, use the initial guess to match the lines and then find
                   the best fit

       MatchZero--First calculate the zeropoint, then match the lines.

       FullXCor--Provides a full cross correlation for all the coefficients



    """
    ImageSolution = {}

    # run it if only the zeropoint needs to be calculated
    if method == 'Zeropoint':
        func = st.findzeropoint
        ImageSolution = runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=False, oneline=oneline,
                                    rstep=rstep, istart=istart, nrows=nrows, res=res, smooth=smooth, dres=dres, farr=farr,
                                    dsigma=dsigma, dniter=niter, log=log, verbose=verbose, dc=dc, ndstep=ndstep)

    # use a line matching algorithm to match the lines
    # in the image with those in the line list
    if method == 'Matchlines':
        func = st.findwavelengthsolution
        # set wdiff
        try:
            wdiff = mdiff * ws.coef[1]
        except:
            wdiff = mdiff
        wdiff = 20
        ImageSolution = runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=True, oneline=oneline,
                                    rstep=rstep, istart=istart, nrows=nrows, res=res, subback=subback, smooth=smooth, dres=dres, farr=farr,
                                    dsigma=dsigma, dniter=niter, log=log, verbose=verbose, mdiff=mdiff, wdiff=wdiff, sigma=sigma, niter=niter)

    # first fit a zeropoint, then match the lines, and then
    # find the rest of the points by using only the zeropoint
    if method == 'MatchZero':
        func = st.findzeropoint
        ws = runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=False, oneline=True,
                         rstep=rstep, istart=istart, nrows=nrows, res=res, subback=subback, smooth=smooth, dres=dres, farr=farr,
                         dsigma=dsigma, dniter=niter, log=log, verbose=verbose, dc=10, ndstep=20)

        func = st.findwavelengthsolution
        ws = runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=True, oneline=True,
                         rstep=rstep, istart=istart, nrows=nrows, res=res, subback=subback, smooth=smooth, dres=dres, farr=farr,
                         dsigma=dsigma, dniter=niter, log=log, verbose=verbose, sigma=sigma, niter=niter)
        func = st.findzeropoint
        ImageSolution = runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=False, oneline=oneline,
                                    rstep=rstep, istart=istart, nrows=nrows, res=res, subback=subback, smooth=smooth, dres=dres, farr=farr,
                                    dsigma=dsigma, dniter=niter, log=log, verbose=verbose, dc=dc, ndstep=ndstep)

    if method == 'FullXcor':
        func = st.findxcor
        if ws is None:
            return None
        dcoef = ws.coef * 0.1
        dcoef[0] = dc
        ImageSolution = runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=True, oneline=oneline,
                                    rstep=rstep, istart=istart, nrows=nrows, res=res, subback=subback, smooth=smooth, dres=dres, farr=farr,
                                    dsigma=dsigma, dniter=niter, log=log, verbose=verbose, dcoef=dcoef, ndstep=ndstep)

    return ImageSolution


def runsolution(xarr, specarr, slines, sfluxes, ws, func, ivar=None,
                fline=True, oneline=False, farr=None, rstep=20,
                istart=None, nrows=1, dsigma=5, dniter=5, subback=0, smooth=0, res=2.0,
                dres=0.1, log=None, verbose=True, **kwargs):
    """Starting in the middle of the image, it will determine the solution
       by working its way out to either edge and compiling all the results into
       ImageSolution.  The image solution is only saved if the sigma is less than dres.

       xarr--Full range in x of pixels to solve for

       specarr--Input 2D flux

       func--function to use for the solution

       fline--True if spectral lines are in array format.  If False, spectral
              lines are assumed to be in line format

       oneline--whether to measure one line or all lines


    """
    # set up the variables
    ImageSolution = {}

    # Setup the central line if it isn't specified
    if istart is None:
        istart = int(0.5 * len(specarr))

    # set up the flux from the central line (or the line specified by the user
    # in istart)
    if farr is None:
        specext = apext.apext(xarr, specarr, ivar=ivar)
        farr = apext.makeflat(specarr, istart, istart + nrows)
        farr = st.flatspectrum(xarr, farr, mode='poly', order=subback)

    # smooth the data
    if smooth > 0:
        farr = st.smooth_spectra(xarr, farr, sigma=smooth)

    # detect the lines
    cxp = st.detect_lines(xarr, farr, dsigma, dniter)
    nlines = len(cxp)

    # first set up the artificial spectrum
    if fline:
        swarr = slines
        sfarr = sfluxes
    else:
        swarr, sfarr = st.makeartificial(
            slines, sfluxes, farr.max(), res, dres)

    # find the solution for the central wavelegnth
    k = istart
    min_lines = 0.1 * len(cxp)
    if oneline:
        mws = solution(xarr, farr, swarr, sfarr, ws, func,
                       min_lines=min_lines, dsigma=dsigma, dniter=dniter, **kwargs)
        for i in range(dniter - 1):
            mws = solution(xarr, farr, swarr, sfarr, mws, func,
                           min_lines=min_lines, dsigma=dsigma, dniter=dniter, **kwargs)

        if verbose and mws is not None:
            msg = "%5i %3i %3.2f" % (k,
                                     mws.func.mask.sum(),
                                     mws.sigma(
                                         mws.func.x,
                                         mws.func.y))
            if log is not None:
                log.message(msg)
        return mws

    # now loop through each step, and calculate the wavelengths for the given
    if log is not None:
        log.message('%5s %3s %4s' % ('Line', 'N', 'RMS'), with_header=False)

    for i in range(0, int(0.5 * len(specarr)), rstep):
        for k in [istart - i, istart + i]:

            if k in ImageSolution.keys():
                continue

            lws = getwsfromIS(k, ImageSolution, default_ws=ws)

            # set up the flux from the set of lines
            farr = apext.makeflat(specarr, k, k + nrows)

            if smooth > 0:
                farr = st.smooth_spectra(xarr, farr, sigma=smooth)

            # continuum correct the spectrum if possible
            try:
                farr = st.flatspectrum(xarr, farr, mode='poly', order=subback)
            except:
                continue

            # find the solution to the lines
            fws = solution(xarr, farr, swarr, sfarr, lws, func,
                           min_lines=min_lines, dsigma=dsigma, dniter=dniter, **kwargs)
            if fws is not None:
                if fws.sigma(fws.func.x, fws.func.y) < dres:
                    ImageSolution[k] = fws

                if verbose:
                    p_new = i * 100.0 / (0.5 * len(specarr))
                    # ctext='Percentage Complete: %d %d %f\r' % (i,p_new, time.clock()) #p_new
                    # sys.stdout.write(ctext)#
                    # sys.stdout.flush()
                    msg = "%5i %3i %3.2f" % (k,
                                             fws.func.mask.sum(),
                                             fws.sigma(
                                                 fws.func.x,
                                                 fws.func.y))
                    if log is not None:
                        log.message(msg, with_header=False)

    return ImageSolution


def solution(xarr, farr, sl, sf, ws, func, min_lines=2,
             dsigma=5, dniter=3, pad=50, **kwargs):
    """Extract a single line and calculate the wavelneght solution"""

    # check to see if there are any points
    xp = st.detect_lines(xarr, farr, dsigma, dniter)

    if len(xp) > min_lines and ws:
        # make the artificial list
        wmin = ws.value(xarr.min())
        wmax = ws.value(xarr.max())
        smask = (sl > wmin - pad) * (sl < wmax + pad)

        # fit the function
        try:
            fws = func(xarr, farr, sl[smask], sf[smask], ws, **kwargs)
            pass
        except SALTSpecError as e:
            return None
        except TypeError as e:
            return None
        except IndexError as e:
            return None
        except Exception as e:
            print e
            return None

        return fws

    return None


def getwsfromIS(k, ImageSolution, default_ws=None):
    """From the imageSolution dictionary, find the ws which is nearest to
       the value k

    """
    if len(ImageSolution) == 0:
        return default_ws
    ISkeys = np.array(ImageSolution.keys())
    ws = ImageSolution[ISkeys[abs(ISkeys - k).argmin()]]
    if ws is None:
        dist = abs(ISkeys[0] - k)
        ws = ImageSolution[ISkeys[0]]
        for i in ISkeys:
            if ImageSolution[i] and abs(i - k) < dist:
                dist = abs(i - k)
                ws = ImageSolution[i]
    return ws
