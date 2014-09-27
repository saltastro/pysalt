.. _slotutcfix:

**********
slotutcfix
**********


Name
====

slotreadtimefix -- Correct UTC times in slotmode data headers

Usage
=====

slotreadtimefix images outimages outpref
(clobber) (logfile) (verbose)

Parameters
==========


*images*
    String. List of images to reduce. Data can be provided as a comma-delineated
    list, or a string with a wildcard (e.g. 'images=S20061210*.fits'), or
    a foreign file containing an ascii list of image filenames. For ascii
    list option, the filename containing the list must be provided
    preceded by a '@' character, e.g. 'images=@listoffiles.lis'. Note
    that SLOT mode fits files often contain more than one exposed frame.

*outimage*
    String. A list of images. Data can be provided as a comma-separated
    list, or a string with a wildcard (e.g. 'outimages=rS20061210*.fits'), or
    a foreign file containing an ascii list of image filenames. For ascii
    list option, the filename containing the list must be provided
    preceded by a '@' character, e.g. 'outimages=@listoffiles.lis'. This list
    must be of the same size as the images argument list.

*outpref*
    String. If the outpref string is non-zero in length and contains
    characters other than a blank space, it will override any value of the
    outimages argument. Output file names will use the name list provided
    in the images argument, but adding a prefix to the basename of
    each  output file defined by outpref. An absolute or relative directory
    path can be included in the prefix, e.g. 'outpref=/Volumes/data/p'.

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

.ih
===

DESCRIPTION

SLOTREADTIME is a tool to fix the UTC time.  The time in the binary
files is the readtime and not the start time of the exposure.  The actually time of the exposure is 7xEXPTIME less than the time of the read out as the image is shifted down 7 places on the CCD before being read out.

This program will read in a file, and for each extension correct the UTC-TIME keyword by the time it takes to shift the exposure so the UTC-TIME corresponds to the start of the exposure and adds a keyword with READTIME.

It will not run on files with the READTIME header keyword already in place.


Examples
========

1. To automatically correct the UTC values in the image``::

    --> slotreadtimefix images="*.fits" outimages='' outpref='r'
    clobber=True logfile=salt.log verbose=y

Time requirements
=================

A linux machine with 2 GB of RAM and a 2.8 Ghz processer was able to
process 12 SALTICAM slotmode exposures with 200 extensions in 16 seconds.

Bugs and limitations
====================


Send feedback and bug reports to salthelp@saao.ac.za

See also
========

 :ref:`saltslot` :ref:`slotview`