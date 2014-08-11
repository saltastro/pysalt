.. _saltfpcalring:

*************
saltfpcalring
*************


Name
====

saltfpcalring -- Create the calibration data file

Usage
=====

saltfpcalring images outfile waves (method) (thresh) (niter) (minsize) (axc) (ayc)
(clobber) (logfile) (verbose)

Parameters
==========


*images*
    String. Name of the input ring image, or the name of a file containing a list of ring input images. If a list file is used, it should contain one ring input file per line. In this case the name of the list file should be preceded with an ``@`` character to indicate that a file list is being input. No wildcards are accepted. The ring image/s must have been processed by saltprepare or saltfpprep.

*outfile*
    String. Output file which stores the fitting results. If this is a new file a comment line (up to 80 characters) may be added. If the output file already exists, the results of SALTFPCALRING's fits will be appended to the end of the file.

*waves*
    String. List of wavelengths for the calibrate arcs.   For each image, there should be one value for the wavelength of
    the arc line in the image

*method*
    String. Method for determining ring center.  If None, then axc and ayc must be specified.

*thresh*
    Float. Threshhold for detecting ring features above the background in the image

*minsize*
    Int.  Minimize size for detecting a ring in the image.

*niter*
    Integer. Maximum number of iterations to be conducted when determining the ring centre. If the routine does not converge upon a ring centre during itmax iterations the routine will stop recalculating the centre position at itmax iterations.

*conv*
    Float. This value determines the convergence criteria for the ring fitting. Convergence is said to occur when the correction required for the ring centre position between successive iterations is less than conv times the maximum of the estimated error in the ring centre.

*axc*
    Float. X-value for ring center.   If specified, it will be fixed with this value.  Otherwise, leave blank to have the task determine the ring x-center.  Only one x-center can be specified for all images.

*ayc*
    Float. Y-value for ring center.   If specified, it will be fixed with this value.  Otherwise, leave blank to have the task determine the ring y-center.  Only one y-center can be specified for all images.

*clobber*
    boolean.  If true, overwrite outfile.

*logfile*
    String. Name of an ascii file for storing log and error messages
    written by the task. The file may be new, or messages can also be
    appended to a pre-existing file.

*verbose*
    Boolean. Verbose (yes) or quiet (no) execution of the routine.

Description
===========

The optical geometry of the SALT FPI produces a wavelength gradient over the field; this gradient is circularly symmetric about the etalon's optical axis and varies quadratically with the radial distance from the axis. Thus a uniform monochromatic illumination of the entrance aperture will produce an image with a ring-shaped geometry. These calibration rings simplify the wavelength calibration procedure for the instrument. The program SALTFPCALRING analyses one or more calibration ring images and determines the radius and centre coordinates for each ring using saltfpringfind. These data are recorded in a file for subsequent analysis by the SALTFPZERPOINT program.


Examples
========

1. To determine the ring centres and radii of files listed in the filelist rplist. The config file 'pfp.par' containing initial guesses for axc, ayc, etc is used in this example:

--> saltfpcalring images=@rplist outfile=calring.out waves='6655.5,6655.5,6655.5' logfile=salt.log


Time and disk requirements
==========================

Individual unbinned full frame RSS image files can be 112MB in
size. It is recommended to use workstations with a minimum of 512MB
RAM. On a linux machine with 2.8 Ghz processor and 2 Gb of RAM, one
2051x2051 image in 0.85 sec.


Bugs and limitations
====================

Send feedback and bug reports to salthelp@saao.ac.za

See also
========

 :ref:`saltfpringfind`