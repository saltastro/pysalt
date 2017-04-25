#!/usr/bin/env python
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved. See LICENSE for more details                        #
"""
SPECRECTIFY is a program to read in SALT RSS spectroscopic data and
a wavelength solution for that data.  It will then apply that wavelength
solution to the data set.

For the solution, the user should be able to specify a model, an ascii file, a fits
file, or a set of files.  The ascii file also might have multiple solutions
in it or the program might be given a set of fits files.  It should look
through the solutions and datetime of the solutions and create a solution
which is the interpolation over that time.

Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)    1.0       01 May 2010

TODO
----
1. Implement the fits files for wavelength solution
2. Implement multiple solutions/data


LIMITATIONS
-----------

1.  A mix of ascii and fits files are not supported

2.  The size of the output image is set by the input image

3.  Pixel scale is currently hardwired for RSS

4.  Linear interpolations between steps in the wavelength solution

"""
# Ensure python 2.5 compatibility
from __future__ import with_statement

import os
import string
import sys
import glob
import pyfits
import math
import datetime
import numpy as np

from pyraf import iraf
import saltsafekey as saltkey
import saltsafeio
from saltsafelog import logging, SaltLog

from PySpectrograph.Models import RSSModel

import WavelengthSolution
import spectools as st
from spectools import SALTSpecError


debug = True
import pylab as pl


# -------------------------------------------------------------------------
# core routine

def specrectify(images, outimages, outpref, solfile=None, caltype='line',
                function='polynomial', order=3, inttype='linear', w1=None,
                w2=None, dw=None, nw=None, blank=0, conserve=False, nearest=False,
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
            soldict = entersolution(solfile)
        else:
            soldict = None

        # read in rectify each image
        for img, oimg, in zip(infiles, outfiles):
            if caltype == 'line':
                msg = 'Creating rectified image %s from image %s using files %s' % (
                    oimg, img, solfile)
            else:
                msg = 'Creating rectified image %s from image %s using RSS Model' % (
                    oimg, img)
            log.message(msg)
            hdu = saltsafeio.openfits(img)
            hdu = rectify(hdu, soldict, caltype=caltype, function=function,
                          order=order, inttype=inttype, w1=w1, w2=w2, dw=dw, nw=nw,
                          pixscale=0.0, blank=blank, conserve=conserve, nearest=nearest,
                          clobber=clobber, log=log, verbose=verbose)
            # write out the oimg
            saltsafeio.writefits(hdu, oimg, clobber=clobber)


# -----------------------------------------------------------
# rectify data

def rectify(hdu, soldict, caltype='line', function='poly', order=3, inttype='interp',
            w1=None, w2=None, dw=None, nw=None, blank=0, pixscale=0.0, time_interp=False,
            conserve=False, nearest=False, clobber=True, log=None, verbose=True):
    """Read in an image and a set of wavlength solutions.  Calculate the best
       wavelength solution for a given dataset and then apply that data set to the
       image

     return
    """

    # set the basic values
    set_w1 = (w1 is None)
    set_w2 = (w2 is None)
    set_dw = (dw is None)
    set_nw = (nw is None)

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

    timeobs = enterdatetime('%s %s' % (dateobs, utctime))

    # check to see if there is more than one solution
    if caltype == 'line':
        if len(soldict) == 1:
            sol = soldict.keys()[0]
            slitid = None
            if not matchobservations(
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
            # set up a wavelength solution
            try:
                w_arr = findsol(xarr, soldict, istart, caltype, nearest, timeobs, exptime, instrume, grating, grang, arang, filtername,
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
                w_arr = findsol(xarr, soldict, istart, 'rss', nearest, timeobs, exptime, instrume, grating, grang, arang, filtername,
                                slit, xbin, ybin, slitid, function, order)

            # set up the output x-axis

            if set_w1:
                w1 = w_arr.min()
            if set_w2:
                w2 = w_arr.max()
            if set_nw:
                nw = len(xarr)
            if set_dw:
                dw = float(w2 - w1) / nw
            nw_arr = createoutputxaxis(w1, w2, nw)

            # setup the VARIANCE and BPM frames
            if saltkey.found('VAREXT', hdu[i]):
                varext = saltkey.get('VAREXT', hdu[i])
            else:
                varext = None

            # setup the BPM frames
            if saltkey.found('BPMEXT', hdu[i]):
                bpmext = saltkey.get('BPMEXT', hdu[i])
            else:
                bpmext = None

            # for each line in the data, determine the wavelength solution
            # for a given line in the image
            for j in range(len(hdu[i].data)):
                # find the wavelength solution for the data
                w_arr = findsol(xarr, soldict, j, caltype, nearest, timeobs, exptime, instrume, grating, grang, arang, filtername,
                                slit, xbin, ybin, slitid, function, order)

                # apply that wavelength solution to the data
                if w_arr is not None:
                    try:
                        hdu[i].data[
                            j,
                            :] = st.interpolate(
                            nw_arr,
                            w_arr,
                            hdu[i].data[
                                j,
                                :],
                            inttype,
                            left=blank,
                            right=blank)
                    except Exception as e:
                        hdu[i].data[j, :] = hdu[i].data[j, :] * 0.0 + blank
                        msg = 'In row %i, solution cannot be found due to %s' % (
                            i, e)

                    # correct the variance frame
                    if varext:
                        try:
                            hdu[varext].data[
                                j,
                                :] = st.interpolate(
                                nw_arr,
                                w_arr,
                                hdu[varext].data[
                                    j,
                                    :],
                                inttype,
                                left=blank,
                                right=blank)
                        except Exception as e:
                            msg = 'In row %i, solution cannot be found due to %s' % (
                                i, e)

                    # correct the BPM frame
                    if bpmext:
                        try:
                            hdu[bpmext].data[
                                j,
                                :] = st.interpolate(
                                nw_arr,
                                w_arr,
                                hdu[bpmext].data[
                                    j,
                                    :],
                                inttype,
                                left=blank,
                                right=blank)
                        except Exception as e:
                            msg = 'In row %i, solution cannot be found due to %s' % (
                                i, e)
                else:
                    hdu[i].data[j, :] = hdu[i].data[j, :] * 0.0 + blank

            if conserve:
                hdu[i].data = hdu[i].data / dw
                if varext:
                    hdu[varext].data = hdu[varext].data / dw

            # Add WCS information
            saltkey.new('CTYPE1', 'LAMBDA', 'Coordinate Type', hdu[i])
            saltkey.new('CTYPE2', 'PIXEL', 'Coordinate Type', hdu[i])
            saltkey.new(
                'CD1_1',
                dw,
                'WCS: Wavelength Dispersion in angstrom/pixel',
                hdu[i])
            saltkey.new('CD2_1', 0.0, 'WCS: ', hdu[i])
            saltkey.new('CD1_2', 0.0, 'WCS: ', hdu[i])
            saltkey.new('CD2_2', ybin * pixscale, 'WCS: ', hdu[i])
            saltkey.new('CRPIX1', 0.0, 'WCS: X Reference pixel', hdu[i])
            saltkey.new('CRPIX2', 0.0, 'WCS: Y Reference pixel', hdu[i])
            saltkey.new('CRVAL1', w1, 'WCS: X Reference pixel', hdu[i])
            saltkey.new('CRVAL2', 0.0, 'WCS: Y Reference pixel', hdu[i])
            saltkey.new('CDELT1', 1.0, 'WCS: X pixel size', hdu[i])
            saltkey.new('CDELT2', 1.0, 'WCS: Y pixel size', hdu[i])
            saltkey.new('DC-FLAG', 0, 'Dispesion Corrected image', hdu[i])

    return hdu


def makeinterpsolution(xarr, soldict, timeobs, exptime,
                       instrume, grating, grang, arang, filtername, slitid):
    """Given a set off solutions taken at different times, interpolate between the two best solutions
       and determine the wavelength solution at each row
    """

    # determine the closest solution and then for each row in there, determine
    # the coefficient
    btime = 1e10
    bsol = None
    for sol in soldict:
        if matchobservations(
                soldict[sol], instrume, grating, grang, arang, filtername, slitid):
            dtime = (soldict[sol][0] - timeobs).seconds
            if dtime == 0.0:
                return None
            else:
                if dtime < btime:
                    btime = dtime
                    bsol = sol
    if bsol is None:
        return None

    # now loop through that one and create the coefficentis for each y_i and
    # row of that value
    wsrow = soldict[bsol][9]
    wscoef = []
    for y_i in wsrow:
        function, order, coef, domain = findlinesol(
            soldict, y_i, nearest, timeobs, exptime, instrume, grating, grang, arang, filtername, slitid, xarr)
        wscoef.append(coef)

    return [timeobs, instrume, grating, grang, arang,
            filtername, None, function, order, wsrow, wscoef, domain]


def findsol(xarr, soldict, y_i, caltype, nearest, timeobs, exptime, instrume, grating, grang, arang, filtername,
            slit, xbin, ybin, slitid, function, order):
    """Find the wavelength solution.  Either from a database containing all the
       solutions adn calculate the best one or calculate based on the model
       for RSS
    """

    if caltype == 'line':
        function, order, coef, domain = findlinesol(
            soldict, y_i, nearest, timeobs, exptime, instrume, grating, grang, arang, filtername, slitid, xarr)
        if function is None:
            msg = 'No solution matches for %s, %s, %s, %s, %s, %s' % (
                instrume, grating, grang, arang, filtername, slitid)
            raise SALTSpecError(msg)
        if coef is None:
            return None
        order = int(order)
        ws = WavelengthSolution.WavelengthSolution(
            xarr,
            xarr,
            function=function,
            order=order) 
        ws.func.func.domain = domain
        ws.set_coef(coef)
        w_arr = ws.value(xarr)

    elif caltype == 'rss':
        w_arr = calcsol(
            xarr,
            y_i,
            instrume,
            grating,
            grang,
            arang,
            filtername,
            slit,
            xbin,
            ybin,
            function,
            order)
    else:
        message = 'SALTRECTIFY--Invalid caltype'
        raise SALTSpecError(message)

    return w_arr

# -----------------------------------------------------------------------------
# set the output axis wavelenght values


def createoutputxaxis(wstart, wend, nw):
    """Create the output array for the axis for the data to be transformed onto.
    """

    try:
        warr = np.linspace(wstart, wend, num=nw)
    except Exception as e:
        msg = 'No output array was create because %s' % e
        raise SALTSpecError(msg)
    return warr


def calcsol(xarr, y_i, instrume, grating, grang, arang,
            filtername, slit, xbin, ybin, function, order):
    rss = RSSModel.RSSModel(grating_name=grating.strip(), gratang=grang,
                            camang=arang, slit=slit, xbin=xbin, ybin=ybin)
    gamma = 180.0 / math.pi * math.atan((y_i * rss.detector.pix_size * rss.detector.ybin
                                         - 0.5 * rss.detector.find_height()) / rss.camera.focallength)
    # y_i/rssmodel.detector.height*rssmodel
    d = rss.detector.xbin * rss.detector.pix_size * \
        (xarr - rss.detector.get_xpixcenter())
    alpha = rss.alpha()
    beta = -rss.beta()
    dbeta = -np.degrees(np.arctan(d / rss.camera.focallength))
    w_arr = 1e7 * rss.calc_wavelength(alpha, beta + dbeta, gamma=gamma)

    return w_arr


# -----------------------------------------------------------------------------
# Find the best wavelength solution

def findlinesol(soldict, yc, nearest, timeobs, exptime,
                instrume, grating, grang, arang, filtername, slitid, xarr=None):
    """Find the best wavelength solution given a datetime of the observation
       and a wavelenght line.  We will assume that the best solution is found
       between the two solutions which are on either side of the timeobs or
       otherwise the closest time.

    """

    # if there is only one return that one
    if len(soldict) == 1:
        sol = soldict.keys()[0]
        function = soldict[sol][7]
        order = soldict[sol][8]
        coef = findcoef(yc, soldict[sol][9], soldict[sol][10])
        domain = soldict[sol][11]
        return function, order, coef, domain

    # first find the closest wavelength solution in time and use that for
    # a base sample
    time_list = []
    function_list = []
    order_list = []
    domain_list = []
    coef_list = []
    for sol in soldict:
        if matchobservations(
                soldict[sol], instrume, grating, grang, arang, filtername, slitid):
            # right at the same time
            if (soldict[sol][0] - timeobs).seconds == 0.0:
                function = soldict[sol][7]
                order = soldict[sol][8]
                coef = findcoef(yc, soldict[sol][9], soldict[sol][10])
                domain = soldict[sol][11]
                return function, order, coef, domain
            else:
                t = timeobs + datetime.timedelta(seconds=exptime)
                time_list.append(subtracttime(soldict[sol][0], t))
                function_list.append(soldict[sol][7])
                order_list.append(soldict[sol][8])
                domain_list.append(soldict[sol][11])
                coef_list.append(
                    findcoef(
                        yc,
                        soldict[sol][9],
                        soldict[sol][10]))

    if len(time_list) == 1:
        function = function_list[0]
        order = order_list[0]
        domain = domain_list[0]
        coef = coef_list[0]
    elif len(time_list) < 1:
        return None, None, None, None
    else:
        # determine the points closest in time to the observation
        t_arr = np.array(time_list)
        id_min = t_arr.argmin()
        function = function_list[id_min]
        order = order_list[id_min]
        domain = domain_list[id_min]
        coef = coef_list[id_min]

        # If xarr is None, return the solution closest in time to the
        # obsevation
        if coef is None:
            return function, order, coef, domain

        # if nearest, return that
        if nearest:
            return function, order, coef, domain

        # If xarr, then calculation w_arr for each of the solutions available,
        # After calculating the solution, calculated the time weighted average
        # wavelength calibration for that pixel array.  Then recalculate the solution
        # using the same funcitonal form as the best function
        w_ave = xarr * 0.0
        wei = 0.0
        for i, t in enumerate(time_list):
            if t == 0:
                return function, order, coef
            ws = WavelengthSolution.WavelengthSolution(
                xarr,
                xarr,
                function=function_list[i],
                order=order_list[i])
            if function in ['legendre', 'chebyshev']:
                ws.func.func.domain=domain_list[i] 
            ws.set_coef(coef_list[i])
            try:
                w_ave += ws.value(xarr) / t ** 2
                wei += 1.0 / t ** 2
            except:
                pass
        w_ave = w_ave / wei
        # for the purposes of speed, we can sparsely sample the arrays

        j = int(len(xarr) / order / 8)
        ws = WavelengthSolution.WavelengthSolution(
            xarr[
                ::j], w_ave[
                ::j], function=function, order=order)
        if function in ['legendre', 'chebyshev']:
           ws.func.func.domain=[xarr.min(), xarr.max()]
        ws.fit()
        coef = ws.coef

    return function, order, coef, domain

# -----------------------------------------------------------------------------
# Determine if two observations are the same


def matchobservations(
        sol, instrume, grating, grang, arang, filtername, slitid):
    if sol[1] != instrume:
        return False
    if sol[2] != grating:
        return False
    if abs(sol[3] - grang) > 0.01:
        return False
    if abs(sol[4] - arang) > 0.01:
        return False
    if sol[5] != filtername:
        return False
    if slitid is not None:
        if sol[6] != slitid:
            return False
    return True

# -----------------------------------------------------------------------------
# Find the best coefficient solution


def findcoef(yc, wsrow, wscoef):
    """Find the coefficient given a row.  If yc can be found in wsrow, then that
       coefficient is returned.  Otherwise, it will linearly interpoloate between
       the rows and find the best coefficient
    """
    if len(wscoef) == 1:
        return wscoef[0]

    # If it is off the edges, then return None
    if yc < wsrow.min():
        return None
    if yc > wsrow.max():
        return None

    # if it is in the list, then return that coef
    dist = abs(wsrow - yc)
    if dist.min() == 0:
        pass
        # return wscoef[dist.argmin()]

    # if it isn't in the list, then compute by taking
    # the interpolation between the two end points
    coefarr = np.array(wscoef)
    order = len(coefarr[0])
    coef = np.zeros(order, dtype=float)
    for k in range(order):
        tk = np.polyfit(wsrow, coefarr[:, k], 2)
        coef[k] = np.polyval(tk, yc)

    return coef


def entersolution(solfiles):
    """For each solution in the list, enter the information about that solution
       into the file.  The following information will be entered in:
       name, time, grating, grating angle, camera angle, function, order, solution

       return dict
    """
    soldict = {}

    if isinstance(solfiles, str):
        solfile = solfiles
        if solfile[-4:] == 'fits':
            print 'Not supported yet'
        else:
            soldict = readsolascii(solfile, soldict)
    else:
        for solfile in solfiles:
            if solfile[-4:] == 'fits':
                print 'Not supported yet'
            else:
                soldict = readsolascii(solfile, soldict)

    return soldict


def readsolascii(solfile, soldict):
    """Read in an ascii file that has the wavelength solution in formation in it.
       The ascii file is expected to be in the same format as the output from
       saltidentify

       return soldict
    """
    # open up the ascii file
    fin = saltsafeio.openascii(solfile, 'r')
    data = fin.read()

    # Each new solution in the file should have WS at the start of it.  We will split
    # the file on ws
    sol_list = data.split('#WS')
    for sol in sol_list[1:]:
        name = findwskeyword('name', sol)
        timeobs = enterdatetime(findwskeyword('time-obs', sol))
        instrume = findwskeyword('instrument', sol)
        grating = findwskeyword('grating', sol)
        graang = float(findwskeyword('gratilt', sol))
        arang = float(findwskeyword('artilt', sol))
        filtername = findwskeyword('filter', sol)
        function = findwskeyword('Function', sol)
        order = int(findwskeyword('Order', sol))
        try:
            slitid = findwskeyword('slitid', sol)
            name = name + '_' + slitid
        except:
            slitid = None
        try:
            domain = findwskeyword('domain', sol)
            domain = [float(x) for x in domain.split(',')]
        except:
            domain = None
        wsrow = []
        wscoef = []
        j = finddata(sol)
        for ws in sol[j:].splitlines():
            if ws.strip():
                ws = ws.split()
                wsrow.append(float(ws[0]))
                wscoef.append(np.array(ws[1:], dtype=float))
        soldict[name] = [
            timeobs,
            instrume,
            grating,
            graang,
            arang,
            filtername,
            slitid,
            function,
            order,
            np.array(wsrow),
            wscoef,
            domain]

    return soldict


def finddata(sol):
    """Find where the data starts in the solution"""
    i = sol.index('#Starting Data')
    j = sol[i:].index('\n')
    return i + j + 1


def findwskeyword(keyword, sol):
    """Find and return a value for a keyword in the list of the wavelength solution"""
    i = sol.index(keyword)
    j = sol[i:].index('\n')
    return sol[i:i + j].split('=')[1].strip()


def enterdatetime(dstr):
    """Break up the datetime string to create a datetime object
       return datetime
    """
    dlist = dstr.split()
    year, month, day = dlist[0].split('-')
    hour, minute, second = dlist[1].split(':')

    return datetime.datetime(
        int(year), int(month), int(day), int(hour), int(minute), int(float(second)))


def subtracttime(d1, d2):
    """Return the difference in two dates in seconds"""
    dt = max(d1, d2) - min(d1, d2)
    return 86400 * dt.days + dt.seconds


# main code
if not iraf.deftask('specrectify'):
    parfile = iraf.osfn("saltspec$specrectify.par")
    t = iraf.IrafTaskFactory(taskname="specrectify", value=parfile, 
                             function=specrectify, pkgname='saltspec')
