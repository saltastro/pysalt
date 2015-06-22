#!/usr/bin/env python
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.                                                     #
"""
SPECSENS calulates the calibration curve given an observation, a standard star,
and the extinction curve for the site.  The task assumes a 1-D spectrum that
has already been sensed from the original observations.

Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)    1.0       21 Mar 2011

TODO
----

LIMITATIONS
-----------

"""
# Ensure python 2.5 compatibility
from __future__ import with_statement

import os
import sys
import time
import numpy as np
import pyfits

from matplotlib.pyplot import *

from pyraf import iraf
import saltstat
import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging
import spectools as st
from spectools import SALTSpecError

from PySpectrograph.Spectra import Spectrum
from saltfit import interfit

from pylab import *

debug = True


# -----------------------------------------------------------
# core routine

def specsens(specfile, outfile, stdfile, extfile, airmass=None, exptime=None,
             stdzp=3.68e-20, function='polynomial', order=3, thresh=3, niter=5,
             fitter='gaussian', clobber=True, logfile='salt.log', verbose=True):

    with logging(logfile, debug) as log:

        # read in the specfile and create a spectrum object
        obs_spectra = st.readspectrum(specfile.strip(), error=True, ftype='ascii')

        # smooth the observed spectrum
        # read in the std file and convert from magnitudes to fnu
        # then convert it to fwave (ergs/s/cm2/A)
        std_spectra = st.readspectrum(stdfile.strip(), error=False, ftype='ascii')
        std_spectra.flux = Spectrum.magtoflux(std_spectra.flux, stdzp)
        std_spectra.flux = Spectrum.fnutofwave(
            std_spectra.wavelength, std_spectra.flux)

        # Get the typical bandpass of the standard star,
        std_bandpass = np.diff(std_spectra.wavelength).mean()
        # Smooth the observed spectrum to that bandpass
        obs_spectra.flux = st.boxcar_smooth(obs_spectra, std_bandpass)
        # read in the extinction file (leave in magnitudes)
        ext_spectra = st.readspectrum(extfile.strip(), error=False, ftype='ascii')

        # determine the airmass if not specified
        if saltio.checkfornone(airmass) is None:
            message = 'Airmass was not supplied'
            raise SALTSpecError(message)

        # determine the exptime if not specified
        if saltio.checkfornone(exptime) is None:
            message = 'Exposure Time was not supplied'
            raise SALTSpecError(message)

        # calculate the calibrated spectra
        log.message('Calculating the calibration curve for %s' % specfile)
        cal_spectra = sensfunc(
            obs_spectra, std_spectra, ext_spectra, airmass, exptime)

        # plot(cal_spectra.wavelength, cal_spectra.flux * std_spectra.flux)
        # fit the spectra--first take a first cut of the spectra
        # using the median absolute deviation to throw away bad points
        cmed = np.median(cal_spectra.flux)
        cmad = saltstat.mad(cal_spectra.flux)
        mask = (abs(cal_spectra.flux - cmed) < thresh * cmad)
        mask = np.logical_and(mask, (cal_spectra.flux > 0))

        # now fit the data
        # Fit using a gaussian process.
        if fitter=='gaussian':
            from sklearn.gaussian_process import GaussianProcess
            #Instanciate a Gaussian Process model

            dy = obs_spectra.var[mask] ** 0.5
            dy /= obs_spectra.flux[mask] / cal_spectra.flux[mask]
            y = cal_spectra.flux[mask]
            gp = GaussianProcess(corr='squared_exponential', theta0=1e-2,
                                 thetaL=1e-4, thetaU=0.1, nugget=(dy / y) ** 2.0)
            X = np.atleast_2d(cal_spectra.wavelength[mask]).T
            # Fit to data using Maximum Likelihood Estimation of the parameters
            gp.fit(X, y)
    
            x = np.atleast_2d(cal_spectra.wavelength).T
            # Make the prediction on the meshed x-axis (ask for MSE as well)
            y_pred = gp.predict(x)

            cal_spectra.flux = y_pred

        else:
            fit=interfit(cal_spectra.wavelength[mask], cal_spectra.flux[mask], function=function, order=order, thresh=thresh, niter=niter)
            fit.interfit()
            cal_spectra.flux=fit(cal_spectra.wavelength)

        # write the spectra out
        st.writespectrum(cal_spectra, outfile, ftype='ascii')


def sensfunc(obs_spectra, std_spectra, ext_spectra, airmass, exptime):
    """Given an observe spectra, calculate the calibration curve for the
       spectra.  All data is interpolated to the binning of the obs_spectra.
       The calibrated spectra is then calculated from
       C =  F_obs/ F_std / 10**(-0.4*A*E)/T/dW
       where F_obs is the observed flux from the source,  F_std  is the
       standard spectra, A is the airmass, E is the
       extinction in mags, T is the exposure time and dW is the bandpass

    Parameters
    -----------
    obs_spectra--spectrum of the observed star (counts/A)
    std_spectra--know spectrum of the standard star (ergs/s/cm2/A)
    ext_spectra--spectrum of the extinction curve (in mags)
    airmass--airmass of the observations
    exptime--exposure time of the observations
    function
    """

    # re-interpt the std_spectra over the same wavelength
    std_spectra.interp(obs_spectra.wavelength)

    # re-interp the ext_spetra over the same wavelength
    ext_spectra.interp(obs_spectra.wavelength)

    # create the calibration spectra
    cal_spectra = Spectrum.Spectrum(
        obs_spectra.wavelength, obs_spectra.flux.copy(), stype='continuum')

    # set up the bandpass
    bandpass = np.diff(obs_spectra.wavelength).mean()

    # correct for extinction
    cal_spectra.flux = cal_spectra.flux / \
        10 ** (-0.4 * airmass * ext_spectra.flux)

    # correct for the exposure time and calculation the sensitivity curve
    cal_spectra.flux = cal_spectra.flux / exptime / bandpass / std_spectra.flux

    return cal_spectra


# main code
parfile = iraf.osfn("saltspec$specsens.par")
t = iraf.IrafTaskFactory(
    taskname="specsens", value=parfile, function=specsens, pkgname='saltspec')
