.. _salt2iraf:

*********
salt2iraf
*********


Name
====

salt2iraf -- Convert from MEF to single extension FITS files

Usage
=====

salt2iraf images outimages outpref (clobber) (logfile) (verbose)

Parameters
==========


*images*
    String. List of FITS images to prepare. Data can be provided as a
    comma-separated list, or a string with a wildcard
    (e.g. 'images=S20061210*.fits'), or a foreign file containing an ascii
    list of image filenames. For the ascii list option, the filename
    containing the list must be provided preceded by a '@' character,
    e.g. 'images=@listoffiles.lis'. The list can contain data files from
    multiple SALT instruments.

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
    in the images argument, but adding a prefix to each output file
    defined by outpref. An absolute or relative directory path can be
    included in the prefix, e.g. 'outpref=/Volumes/data/p'.

*(ext)*
    Hidden Int.   The extension number to copy to the Primary Extension.

*(clobber)*
        Hidden boolean. If set to 'yes' files contained within the outpath
        directory will be overwritten by newly created files of the same
        name.

*(logfile)*
        Hidden String. Name of an ascii file for storing log and error messages
        from the tool. The file may be new, or messages can also be appended to a
        pre-existing file.

*(verbose)*
        Hidden Boolean. If verbose=n, log messages will be suppressed.

Description
===========

SALT2IRAF converts PySALT Multi-Extension FITS (MEF) files into
single extension FITS files.   Given a list of images to convert,
it will copy the Primary HDU header information and the data file
appearing in the first extension to the Primary HDU of the output
file.  It will also add any additional information in the first
extension that does not appear in the Primary HDU.  These files
will then be compatible with basic IRAF tasks and other tools.

EXAMPLES
1. To convere raw FITS files residing in the /Volumes/data1/ and create
new files with paths and names stored in the ascii list outimages.lis::

    --> salt2iraf images='mbxpP*.fits' outimages='' outpref='i'
    clobber=y logfile='salt.log' verobse=y

TIME AND DISK REQUIREMENTS
Individual unbinned raw full-frame RSS files can be 112MB in size. It is
recommended to use workstations with a minimum of 512MB RAM. On a
linux machine with 2.8 Ghz processor and 2 Gb of RAM, one 2051x2051 image
in 0.08 sec.

Bugs and limitations
====================

The task will raise an error if more than one 'SCI' extension exists
in the images.   If the data is already in single extension, then the
task will copy the image to the output image.

Send feedback and bug reports to salthelp@saao.ac.za

See also
========

 :ref:`.endhelp`