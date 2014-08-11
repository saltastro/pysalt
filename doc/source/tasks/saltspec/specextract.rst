.. _specextract:

***********
specextract
***********


Name
====

specextract -- Extract 1-D spectra from 2-D image

Usage
=====

specextract images outfiles (method) (section) (thresh) (minsize)
(outformat) (clobber) (logfile) (verbose)

Parameters
==========


*images*
    String. List of FITS images to prepare. Data can be provided as a
    comma-separated list, or a string with a wildcard
    (e.g. 'images=S20061210*.fits'), or a foreign file containing an ascii
    list of image filenames. For the ascii list option, the filename
    containing the list must be provided preceded by a '@' character,
    e.g. 'images=@listoffiles.lis'.

*outfile*
    String. An output file for the extracted spectra.  If outformat is
    set to 'FITS', the output format will be a fits file.  If outformat
    is set to 'ascii', the output format will be an ascii file.  Either
    a single outfile can be supplied or a list of outfiles can be
    supplied.   Data can be provided as a
    comma-separated list, or a string with a wildcard
    (e.g. 'images=S20061210*.fits'), or a foreign file containing an ascii
    list of image filenames. For the ascii list option, the filename
    containing the list must be provided preceded by a '@' character,
    e.g. 'images=@listoffiles.lis'.

*(method)*
    String.  The extraction method to use for creating the 1D spectra.  If
    this is set to 'normal', a straight sum of the pixels will be used to
    extract the counts.  If set to 'weighted', a inverse variance weighted
    average will be used.

*(section)*
    String.  Section to extract.  The should be in the format of the upper
    and lower edges separated by a colon and enclosed in brackets, e.g. '[900:1100]'.
    If None, the task will find all objects in the image and extract them.  Multiple
    extraction windows can be given by separating each one by a comma, e.g.
    '[900:1000, 1050:1150, 1200:1300]'.

*(thresh)*
    Real.  Threshold for detecting sources in the 2D spectra.

*(minsize)*
    Real.  Minimum size in pixels of a detect source.   Sources smaller than this
    will be rejected.

*(outformat)*
    String.  Format for outfile.   It can either be set to 'FITS' or 'ascii'.

*(ext)*
    .  Extension to extract spectrum.

*(convert)*
    Hidden boolean. If set to 'yes' wavelength will be converted from
    air to vacuum wavelengths.

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


SPECEXTRACT will create 1D spectra from a 2D image.  The task will
extract the spectrum from coordinates which are provided by the user
or will identify objects and extract them from the spectra.

If section is set to None, SPECTEXTRACT will identify objects in the image.
It will calculate the median array along the spectra axis, and then detect
sources that are three sigma above the background along that array.   The
program will also deblend sources that are detected.

For each region specified by section or each object detected, SPECEXTRACT
will extract that region from the image.   In 'normal' mode, the task
will extract the summed counts from tehe region.   At this time,
this is the only mode supported.

The user has two options for the output format.  Outformat can
either be set to an ascii file or a FITS table.  If it is an ascii file,
a single file will be created with the first column being wavelength, the
second column being counts, and the third column being the error on those
counts.   The error is either determined either from the variance frame
or from the square root of the counts if no variance frame exists.   If more
than one source exists is identified in the image, than

EXAMPLES
1. To extra a spectrum from a 2D image with specextact::

    --> specextract images='pmbxpP201104250029.fits' outfile='pmbxpP201104250029.txt'
    method='normal' section='[900:1100]'

Time and disk requirements
==========================

Individual unbinned raw full-frame RSS files can be 112MB in size. It is
recommended to use workstations with a minimum of 512MB RAM. On a
linux machine with 2.8 Ghz processor and 2 Gb of RAM, one 2051x2051 image
in 0.15 sec.

Bugs and limitations
====================

The weighted option is not currently implimented.

Send feedback and bug reports to salthelp@saao.ac.za

See also
========

 :ref:`specidentify`