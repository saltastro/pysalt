.. _bvittofits:

**********
bvittofits
**********


Name
====

bvittofits -- Converts BVIT sparse data files to standard FITS images

Usage
=====

bvittofits images (prefix) (maxdelay) (tbin) (xbin) (ybin)
(out32bit) (clobber) (logfile) (verbose) (debug)

Parameters
==========


*images*
    String. List of images to reduce. Data can be provided as a comma-delineated
    list, or a string with a wildcard (e.g. 'images=S20061210*.fits'), or
    a foreign file containing an ascii list of image filenames. For ascii
    list option, the filename containing the list must be provided
    preceded by a '@' character, e.g. 'images=@listoffiles.lst'.

*prefix*
    String. Prefix to apply to output numbered FITS images. E.g. prefix='bvit-' will produce ``bvit-1.fits``, ``bvit-2.fits``, etc.

*maxdelay*
    Int. The maximum delay allowed between files in miliseconds. If the delay is larger the timestamp from the header will be taken as a new reference point.

*tbin*
    Real. Binsize in seconds.

*xbin*
    Int. Binsize to use for spatial x-direction. If a value greater than 1 is supplied the image is binned by summing xbin pixels in the x-direction.

*ybin*
    Int. Binsize to use for spatial y-direction. If a value greater than 1 is supplied the image is binned by summing ybin pixels in the y-direction.

*out32bit*
    Boolean. Write output image using 32 bit unsigned integers instead of default 16 bit unsigned integers for higher dynamic range at the cost of doubling the file size.
    By default when summing the counts in a time bin each pixel has a maximum value of ``2**16-1``, when this option is set that value becomes ``2**32-1``.

*clobber*
        Hidden Boolean. If clobber=y the tool is permitted to overwrite an exisiting
        file with name outfile.

*logfile*
        String. Name of an ascii file for storing log and error messages
        from the tool. The file may be new, or messages can also be appended to a
        pre-existing file.

*verbose*
        Boolean. If verbose=n, log messages will be suppressed.

*debug*
        Boolean. If debug=y, will give more debug information if an error occurs (use this option to gather information when reporting a bug).

Description
===========

BVIT uses a sparse matrix format to store data.
`bvittofits` converts files from this format to standard FITS images by summing the counts in each time bin optionally rebinning the resulting image in the spatial direction.
This allows subsequent annalysis of BVIT data using standard photometry tools.
It also calculates the julian date at the center of each time bin.


Examples
========

1. To convert all files listed in images.lst using a 1 second binning, setting a new reference point if the delay is larger than 1 second, use::

    --> bvittofits images="@images.lst" prefix='bvit-' maxdelay=1000
    tbin=1.0 out32bit=no clobber=yes logfile=salt.log verbose=y

Time requirements
=================

A macbook with 2 GB of RAM and a dual core 2.16 Ghz processer was able to
process a dataset at a rate of one 1 second bin per second.

Bugs and limitations
====================

The current version of BVITTOFITS has not been fully tested yet.

Send feedback and bug reports to salthelp@saao.ac.za

See also
========

 :ref:`bvitphot`