# Copyright (c) 2009, South African Astronomical Observatory (SAAO)       #
# All rights reserved.   See LICENSE file for more details                #
"""
SPECSELFID will preform simple zeropoint rectification on the data itself.
It will extract one line to count as the fiducial and then use SPECIDENTIFY
to calculate zeropoint offsets from that line and SPECRECTIFY to produce
a straighten image.

Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)    1.0       01 Oct 2013
"""
import sys
import os
import pyfits
import numpy as np

from pyraf import iraf
from iraf import pysalt

import saltsafeio as saltio
import saltsafekey as saltkey
from saltsafelog import logging, SaltLog


from specidentify import identify
import spectools as st
from PySpectrograph.WavelengthSolution import WavelengthSolution


import spectools as st

debug = False


def specselfid(images, outimages, outpref, refimage=None, ystart='middlerow',
               rstep=3, clobber=False, logfile='salt.log', verbose=True):

    with logging(logfile, debug) as log:

        # set up the variables
        infiles = []
        outfiles = []

        # Check the input images
        infiles = saltio.argunpack('Input', images)

        # create list of output files
        outfiles = saltio.listparse(
            'Outimages',
            outimages,
            outpref,
            infiles,
            '')

        # set up defaults
        if saltio.checkfornone(refimage) is not None:
            rhdu = saltio.openfits(refimage)
        else:
            refimage = None

        # read in rectify each image
        for img, oimg, in zip(infiles, outfiles):
            hdu = saltio.openfits(img)
            log.message(
                'Performing self-identification and rectification on %s' %
                img)
            for i in range(1, len(hdu)):
                if hdu[i].name == 'SCI':
                    if refimage is None:
                        sdata = hdu[i].data
                    else:
                        sdata = rhdu[i].data
                    hdu[i].data = selfid(
                        hdu[i].data,
                        sdata,
                        ystart=ystart,
                        rstep=rstep)
                    if saltkey.found('VAREXT', hdu[i]):
                        varext = saltkey.get('VAREXT', hdu[i])
                        hdu[varext].data = selfid(
                            hdu[varext].data,
                            sdata,
                            ystart=ystart,
                            rstep=rstep)
                    if saltkey.found('BPMEXT', hdu[i]):
                        bpmext = saltkey.get('BPMEXT', hdu[i])
                        hdu[bpmext].data = selfid(
                            hdu[bpmext].data,
                            sdata,
                            ystart=ystart,
                            rstep=rstep)

            # write out the oimg
            saltio.writefits(hdu, oimg, clobber=clobber)


def selfid(data, sdata, ystart='middlerow', rstep=3):
    """Create a wavelength identification and rectification based on the data
       itself.   The first step is
       to extract the line itself and then run specidentify on the data
    """
    # setup the fiducial line
    if ystart == 'middlerow':
        ystart = int(0.5 * len(data))
    else:
        ystart = int(ystart)

    sfluxes = sdata[ystart, :]
    xarr = np.arange(len(data[0]))
    function = 'polynomial'
    order = 2
    blank = 0
    inttype = 'interp'

    ws = WavelengthSolution.WavelengthSolution(
        xarr,
        xarr,
        function=function,
        order=order)
    ws.fit()
    ndstep = 50
    nw_arr = xarr
    for j in range(0, ystart):
        dcoef = [2.0, 0.0, 0.0]
        nws = st.findxcor(
            xarr,
            sdata[
                j,
                :],
            xarr,
            sfluxes,
            ws,
            dcoef=dcoef,
            ndstep=ndstep,
            best=False,
            inttype=inttype,
            debug=False)
        data[
            j,
            :] = zeroshift(
            data[
                j,
                :],
            xarr,
            nws,
            blank=blank,
            inttype=inttype)
    for j in range(ystart, len(data) - 1):
        dcoef = [2.0, 0.0, 0.0]
        nws = st.findxcor(
            xarr,
            sdata[
                j,
                :],
            xarr,
            sfluxes,
            ws,
            dcoef=dcoef,
            ndstep=ndstep,
            best=False,
            inttype=inttype,
            debug=False)
        data[
            j,
            :] = zeroshift(
            data[
                j,
                :],
            xarr,
            nws,
            blank=blank,
            inttype=inttype)
    return data


def zeroshift(data, xarr, nws, blank=0, inttype='interp'):
    """Calculate zeropoint shift and apply to the data
    """
    if nws is not None:
        w_arr = nws.value(xarr)
        data = st.interpolate(
            xarr,
            w_arr,
            data,
            inttype,
            left=blank,
            right=blank)
    return data

if not iraf.deftask('specselfid'):
    parfile = iraf.osfn("saltspec$specselfid.par")
    t = iraf.IrafTaskFactory(
        taskname="specselfid",
        value=parfile,
        function=specselfid,
        pkgname='saltspec')
