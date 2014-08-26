.. _masktool:

********
masktool
********


Name
====

masktool -- Create a MOS Mask file

Usage
=====

masktool catalog images  (logfile) (verbose)

Parameters
==========


*catalog*
    String. Catalog of objects taht should be used to create mask.

*images*
    String. FITS file that should be used to create a mask

*(logfile)*
    String. Name of an ascii file for storing log and error messages
    from the tool. The file may be new, or messages can also be appended to a
    pre-existing file.

*(verbose)*
    Hidden Boolean. If verbose=n, log messages will be suppressed.

Description
===========


PyRAF wrapper for the PySLITMASK tool.   No input needs to be given.


Examples
========

1. To prepare raw FITS files residing in the /Volumes/data1/ and create
new files with paths and names stored in the ascii list outimages.lis::

    --> masktool images='P*.fits'
    outimages='@outimages.lis' outpref='' logfile='salt.log'

2. To prepare raw FITS files residing in the /Volumes/data1/ and create
new files in the directory /Volumes/data2 with prefix 'p', over-write
existing files, create variance frames and badpixelmasks,  and suppress
log messages::

    --> masktool images='/Volumes/data2/P*.fits'
    outimages='' outpref='/Volumes/data2/p'
    createvar='yes', badpixelmask='bpm.fits'
    clobber='yes' logfile='salt.log' verbose='no'

Time and disk requirements
==========================

Individual unbinned raw full-frame RSS files can be 112MB in size. It is
recommended to use workstations with a minimum of 512MB RAM. On a
linux machine with 2.8 Ghz processor and 2 Gb of RAM, one 2051x2051 image
in 0.15 sec.

Bugs and limitations
====================

Currently there is no statistical error propagation in the SALT
pipeline tasks. However, masktool can
create multiple new image extensions containing a) bad
pixel maps and b) varaiance maps for all amplifiers.

Send feedback and bug reports to salthelp@saao.ac.za

See also
========

 :ref:`saltpipe` :ref:`saltclean`