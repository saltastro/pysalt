.. _saltcrclean:

***********
saltcrclean
***********


Name
====

saltcrclean -- Clean  multiple CCD extensions for cosmic rays

Usage
=====

saltcrclean images outimages outpref (crtype) (thresh) (mbox)
(flux_ratio) (bbox) (bthresh) (gain) (rdnoise) (fthresh)
(bfactor) (gbox) (maxiter) (multithread) (clobber) (logfile) (verbose)

Parameters
==========


*images*
    String. List of input images including, if necessary, absolute or
    relative paths to the data. Data can be provided as a comma-separated
    list, or a string with a wildcard (e.g. 'images=S20061210*.fits'), or
    a foreign file containing an ascii list of image filenames. For ascii
    list option, the filename containing the list must be provided
    preceded by a '@' character, e.g. 'images=@listoffiles.lis'. Do not
    provide any files that have already been mosaiced.

*outimages*
    String. A list of images. Data can be provided as a comma-separated
    list, or a string with a wildcard (e.g. 'outimages=rS20061210*.fits'),
    or a foreign file containing an ascii list of image filenames. For the
    ascii list option, the filename containing the list must be provided
    preceded by a '@' character, e.g. 'outimages=@listoffiles.lis'. This
    list must be of the same size as the images argument list. If the
    output is intended for a different directory the absolute or relative
    path must be supplied with the file name.

*outpref*
    String. If the outpref string is non-zero in length and contains characters other than a blank space, it will override any value of the outimages argument. Output file names will use the name list provided in the images argument, but adding a prefix to the basename of each  output file defined by outpref. An absolute or relative directory path can be included in the prefix, e.g. 'outpref=/Volumes/data/p'.

*crtype*
    String. The type of cosmic ray cleaning to perform.  The options are
    fast, median, and edge.  Fast cleaning identifies cosmic rays from
    statistical deviations above the global background.  Median cleaning
    smooths the image first and uses local background determinations Edge
    cleaning looks for cosmic rays using a Laplacian filter.

*(thresh)*
    Real. The detection threshhold for the identification of cosmic rays.
    This is a factor above the background deviation used to detect cosmic
    rays.

*(mbox)*
    Int.  Window size for searching for cosmic rays.  For fast detection,
    a detected pixel will be compared to neighbors within this box.
    For median detection, it is the size of median filter. It is not used
    for edge detection.

*(flux_ratio)*
    Real.   The flux_ratio is used to confirm cosmic rays in the fast
    method.  If the ratio of the flux in the peak pixel to the median
    flux from the neighboring pixels within a window of size mbox around
    that pixel is greater than the flux_ratio, the pixel is identified
    as a cosmic ray.

*(bbox)*
    Int.  For median cleaning, the bbox is the background box size used to
    determine the local background statistics.  The bbox should be equal
    or greater than 10 in order to have sufficient pixel statistics.

*(bthresh)*
    
    Real.  Sigma-clipping threshhold for the determination of the
    background statistics.  Used in the fast and median types of cleaning.

*(gain)*
    Real.  Gain of the CCD image.  This is only used for calculation of
    noise statistics in the edge identification scheme.

*(rdnoise)*
    Real.  Read noise of the CCD image.This is only used for calculation of
    noise statistics in the edge identification scheme.

*(fthresh)*
    Real.  Threshhold for excluding compact sources in the edge detection
    scheme.  Any sources below this threshhold will not be identified as
    cosmic rays.

*(bfactor)*
    Int.  Factor to sub-pixel sample the images.  This is only
    used in edge detection scheme.

*(gbox)*
    Int.  Window size for growing sources.  Additional sources will not be
    searched for if gbox=0.
    #Growth of sources are dependent on the
    #detection method, but sources which are above gthresh will also be
    #flagged as cosmic rays.

*(maxiter)*
    Int.  Maximum number of times to repeat the cosmic ray detection

*(multithread)*
    Bool. If set to 'yes', the program will use a different thread
    for processing each of the SCI extensions in the image

*(clobber)*
    Hidden boolean. If set to 'yes' files contained within the outpath
    directory will be overwritten by newly created files of the same
    name.

*(logfile)*
    String. Name of an ascii file for storing log and error messages
    written by the task. The file may be new, or messages can also be
    appended to a pre-existing file.

*(verbose)*
    Boolean. If verbose=n, log messages will be suppressed.



Description
===========


SALTCRCLEAN is a task capable of cleaning cosmic rays form single or
multi-extension fits data.  The task gives users a choice of three
different methods of cleaning cosmic rays that provide different
levels of performance.  The three choices the user has are between
fast, median and edge cosmic ray detection.

Fast cosmic ray detection performs very fast identification of cosmic
rays in a CCD image.  The algorithm will calculate the sigma-clipped
background statistics, where all pixels above bthresh are ignored from
the calculations.  Next, it will detect all pixels above thresh times
the background deviation.  These pixels will only be identified as
cosmic rays if: (1) They are the peak pixel inside the window size of
mbox, and (2) if the ratio of the median flux of all the pixels within
the window of size mbox to the flux of the identified pixel is less
than flux_ratio.  Pixels identified as cosmic rays will be replaced
with the median value of all of the neighbors within the window
specified by mbox.  The process will repeat for maxiter.
If the grow parameter is set, it will identify all the pixels around an
identified cosmic ray and will then replace all of
the cosmic rays identified with the median of their neighbors.

Median cosmic ray detection detects cosmic rays through median
smoothing the image.  The first thing the algorithm does is measure
the local background statistics for each pixel.  Within a window of
bbox size, the background statistics are calculated for each pixel.
Then the image is median smoothed with a bin size given by mbox.
The median smoothed image is subtracted from the original image, and
any pixel thresh times greater than its local background deviation is
flagged as a cosmic ray.  If gbox > 0, cosmic rays above
gthresh will be identified if they are neighboring an already
identified cosmic ray.  The process will repeat for the number of
iterations set by the maxiter parameter.

Edge cosmic ray detection uses a Laplacian filter to detect cosmic
rays in an image.  This algorithm uses the process outlines in van
Dokkem (2001).  The first step is that the image is sub-pixel sampled.
Then the image is convolved with a Laplacian filter and then returned
to its original pixel sampling.  A noise image is created by median
smoothing the image and applying a realistic estimate for the gain and
rdnoise.  Cosmic rays are then identified from a signal to noise image
after extended sources have been removed.  Finally, compact sources
are removed by setting the fthresh parameter.  Once complete,
additional cosmic rays can be found by setting gbox>0.  The whole
process will iterate until maxiter or until no more cosmic rays are
removed.


Examples
========

1. To clean a sequence of images::

    --> saltcrclean images='/Volumes/data1/bxpP*.fits' outimages=''
    outpref='/Volumes/data2/m'
    clobber='yes' logfile='salt.log' verbose='yes'

Time and disk requirements
==========================

For a fast linux machine with a 2.8 Ghz processor and 2 Gb of RAM,
one 1043x1024 FITS image can be processed with the following times:
fast (4 sec), median (120 sec), and edge (45 sec).  For this test,
the parameters for SALTCRCLEAN were thresh=5, gbox=3, and maxiter=3
for all of the different methods.


Bugs and limitations
====================

None

Send feedback and bug reports to salthelp@saao.ac.za

See also
========

 :ref:`iraf.noao.imred.crutil.cosmicrays` :ref:`iraf.noao.imred.crutil.crmedian`