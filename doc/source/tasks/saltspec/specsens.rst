.. _specsens:

********
specsens
********


Name
====

specsens-- Calculate sensitivity curve

Usage
=====

specsens specfile outfile stdfile extfile (airmass) (exptime)
(stdzp) (function) (order) (thresh) (clobber) (logfile) (verbose)

Parameters
==========


*specfile*
    String. ASCII file contain the spectra of the calibration source.  This
    should be the output from specextract.

*outfile*
    String. Name of an output file to write the calibrated sensitivity curve.

*stdfile*
    String. ASCII file that contains the calibrated magnitudes for the
    source.

*extfile*
    String. ASCII file that contains the extinction curve for the observing
    site.

*(airmass)*
    Real.  Airmass of the observations

*(Exptime)*
    Real.  Exposure time for the observation

*(stdzp)*
    Real.  Zeropoint for the magnitudes listed in the stdfile for converting into fluxes.

*(function)*
    String.  Functional form to fit to the observation.

*(order)*
    Int.  Order of the function to be fit to the observations.

*(thresh)*
    Real.  Threshold for rejecting values in the obserations of the source.

*(clobber)*
    Hidden boolean. If set to 'yes' files contained within the outpath
    directory will be overwritten by newly created files of the same
    name.

*(logfile)*
    String. Name of an ascii file for storing log and error messages
    from the tool. The file may be new, or messages can also be appended to a
    pre-existing file.

*(verbose)*
    Hidden Boolean. If verbose=n, log messages will be suppressed.

Description
===========


SPECSENS calulates the calibration curve given an observation, a standard star,
and the extinction curve for the site.  The task assumes a 1-D spectrum that
has already been sensed from the original observations.

EXAMPLES


Time and disk requirements
==========================



Bugs and limitations
====================


Send feedback and bug reports to salthelp@saao.ac.za

See also
========

 :ref:`speccal`