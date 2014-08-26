.. _specal:

******
specal
******


Name
====

speccal -- Apply sensitivity curve to observations

Usage
=====

specal   specfile outfile calfile extfile (airmass) (exptime)
(clobber) (logfile) (verbose)

Parameters
==========


*specfile*
    String. ASCII file contain the spectra of the calibration source.  This
    should be the output from specextract.

*outfile*
    String. Name of an output file to write the calibrated sensitivity curve.

*calfile*
    String. ASCII file that contains the calibrated magnitudes for the
    source.

*extfile*
    String. ASCII file that contains the extinction curve for the observing
    site.

*(airmass)*
    Real.  Airmass of the observations

*(exptime)*
    Real.  Exposure time for the observation

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


SPECCAL corrections a given observation by a calibration curve
and the extinction curve for the site.  The task assumes a 1-D spectrum that
has already been caled from the original observations.


EXAMPLES


Time and disk requirements
==========================



Bugs and limitations
====================


Send feedback and bug reports to salthelp@saao.ac.za

See also
========

 :ref:`speccal`