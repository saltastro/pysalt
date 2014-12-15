#!/usr/bin/env python
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.                                                     #
"""
SPECCAL corrections a given observation by a calibration curve
and the extinction curve for the site.  The task assumes a 1-D spectrum that
has already been caled from the original observations.


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

from pyraf import iraf
import saltstat
import saltsafekey as saltkey
import saltsafeio as saltio
from saltsafelog import logging
import spectools as st
from spectools import SALTSpecError

from PySpectrograph.Spectra import Spectrum
from PySpectrograph.Utilities.fit import interfit


debug = True


# -----------------------------------------------------------
# core routine

def speccal(specfile, outfile, calfile, extfile, airmass=None, exptime=None,
            clobber=True, logfile='salt.log', verbose=True):

    with logging(logfile, debug) as log:

        # read in the specfile and create a spectrum object
        obs_spectra = st.readspectrum(specfile, error=True, ftype='ascii')

        # read in the std file and convert from magnitudes to fnu
        # then convert it to fwave (ergs/s/cm2/A)
        cal_spectra = st.readspectrum(calfile, error=False, ftype='ascii')

        # read in the extinction file (leave in magnitudes)
        ext_spectra = st.readspectrum(extfile, error=False, ftype='ascii')

        # determine the airmass if not specified
        if saltio.checkfornone(airmass) is None:
            message = 'Airmass was not supplied'
            raise SALTSpecError(message)

        # determine the exptime if not specified
        if saltio.checkfornone(airmass) is None:
            message = 'Exposure Time was not supplied'
            raise SALTSpecError(message)

        # calculate the calibrated spectra
        log.message('Calculating the calibration curve for %s' % specfile)
        error = False
        try:
            if obs_spectra.var is not None:
                error = True
        except:
            error = False
        flux_spectra = calfunc(
            obs_spectra,
            cal_spectra,
            ext_spectra,
            airmass,
            exptime,
            error)

        # write the spectra out
        st.writespectrum(flux_spectra, outfile, ftype='ascii', error=error)


def calfunc(obs_spectra, std_spectra, ext_spectra,
            airmass, exptime, error=False):
    """Given an observe spectra, calculate the calibration curve for the
       spectra.  All data is interpolated to the
       binning of the obs_spectra.  The calibrated spectra is then calculated
       from:
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

    # re-interp the ext_spetra over the sam ewavelength
    ext_spectra.interp(obs_spectra.wavelength)

    # create the calibration spectra
    cal_spectra = Spectrum.Spectrum(
        obs_spectra.wavelength,
        obs_spectra.flux.copy(),
        stype='continuum')

    # set up the bandpass
    bandpass = np.diff(obs_spectra.wavelength).mean()

    # correct for extinction
    cal_spectra.flux = obs_spectra.flux / \
        10 ** (-0.4 * airmass * ext_spectra.flux)

    # correct for the exposure time and calculation the calitivity curve
    cal_spectra.flux = cal_spectra.flux / exptime / bandpass / std_spectra.flux

    # correct the error calc
    if error:
        cal_spectra.var = obs_spectra.var * cal_spectra.flux / obs_spectra.flux

    return cal_spectra


# main code
if not iraf.deftask('speccal'):
    parfile = iraf.osfn("saltspec$speccal.par")
    t = iraf.IrafTaskFactory(
        taskname="speccal",
        value=parfile,
        function=speccal,
        pkgname='saltspec')
