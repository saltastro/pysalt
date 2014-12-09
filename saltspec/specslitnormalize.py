#!/usr/bin/env python
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved. See LICENSE file for more details                   #
"""
SPECSLITNORMALIZE is a program to apply a correction to the data according
to the response of slit lines in the verticle direction.   The response
is either calculated via the slit lines in the image or they can also
be supplied via a file

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

import spectools as st
import mostools as mt

from spectools import SALTSpecError


debug = True


# -----------------------------------------------------------
# core routine

def specslitnormalize(images, outimages, outpref, response=None,
                      response_output=None, order=2, conv=1e-2, niter=20,
                      startext=0, clobber=False, logfile='salt.log',
                      verbose=True):

    with logging(logfile, debug) as log:

        # Check the input images
        infiles = saltio.argunpack('Input', images)

        # create list of output files
        outfiles = saltio.listparse('Outfile', outimages, outpref, infiles, '')

        # read in the response function
        response = saltio.checkfornone(response)
        if response:
            log.message('Loading response from %s' % response)
            response = readresponse(response)

        # Identify the lines in each file
        for img, ofile in zip(infiles, outfiles):

            # open the image
            hdu = saltio.openfits(img)

            for i in range(startext, len(hdu)):
                if hdu[i].name == 'SCI':
                    log.message('Normalizing extension %i in  %s' % (i, img))
                    # things that will change for each slit

                    # set up the data for the source
                    try:
                        data = hdu[i].data
                    except Exception as e:
                        message = \
                            'Unable to read in data array in %s because %s' % \
                            (img, e)
                        raise SALTSpecError(message)

                    if response is None:
                        response = create_response(
                            data,
                            spatial_axis=1,
                            order=order,
                            conv=conv,
                            niter=niter)
                        if response_output:
                            write_response(response, clobber=clobber)
                    else:
                        # add a check that the response is the same shape as
                        # the data
                        if len(response) != data.shape[0]:
                            raise SALTSpecError(
                                'Length of response function does not equal size of image array')

                    # correct the data
                    data = data / response

                    # correct the variance frame
                    if saltkey.found('VAREXT', hdu[i]):
                        vhdu = saltkey.get('VAREXT', hdu[i])
                        hdu[vhdu].data = hdu[vhdu].data / response

                saltio.writefits(hdu, ofile, clobber=clobber)


def readresponse(response_file):
    """Read in the response file"""
    try:
        response = np.loadtxt(response_file, usecols=[0], unpack=True)
    except Exception as e:
        raise SaltIOError('Could not read %s due to %s' % (response_file, e))
    return response.reshape(len(y), 1)


def write_response(response, response_output, clobber):
    if os.path.isfile(response_output) and clobber:
        os.remove(response_output)
    savetxt(response_output, response)


def create_response(data, spatial_axis=1, order=2, conv=1e-2, niter=20):
    """Create the slit response function.   This fits a low order polynomial
       across the spatial dimension of the image and returns the result

       Parameters
       ----------
       data: numpy array
          Two dimensional image array

       spatial_axis: int
          Axis of the spatail data

       order: int
          Order of polynomial fit

       conv: float
          Two dimensional image array

       order: int
          Order of polynomial fit

       Returns
       --------
       response: numpy array
           Response function for the slit illumination

    """
    y = data.sum(axis=spatial_axis)
    x = np.arange(len(y))
    coef = fit_response(x, y, order, conv, niter)
    response = np.polyval(coef, x)
    response = response / response.mean()

    return response.reshape(len(y), 1)


def fit_response(x, y, order=2, conv=1e-2, niter=20):
    """ Fit a polynomial to the response curve
    """
    coef = np.polyfit(x, y, order)
    s = (y - np.polyval(coef, x)).std()
    for i in range(niter):
        mask = (abs(y - np.polyval(coef, x)) < 3 * s)
        coef = np.polyfit(x[mask], y[mask], order)
        s1 = (y[mask] - np.polyval(coef, x[mask])).std()
        if abs(s1 - s) < conv:
            break
        s = s1
    return coef


# main code

parfile = iraf.osfn("saltspec$specslitnormalize.par")
t = iraf.IrafTaskFactory(
    taskname="specslitnormalize",
    value=parfile,
    function=specslitnormalize,
    pkgname='saltspec')
