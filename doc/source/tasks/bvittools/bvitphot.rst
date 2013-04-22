.. _bvitphot:

********
bvitphot
********


Name
====

bvitphot -- Photometry from SALT BVIT image data

Usage
=====

bvitphot images outfile (outtype) srcfile imstart imend
tstep maxdelay tgt_btype cmp_btype
(clobber) (logfile) (verbose) (debug)

Parameters
==========


*images*
    String. List of images to bin. Data can be provided as a comma-delineated
    list, or a string with a wildcard (e.g. 'images=S20061210*.fits'), or
    a foreign file containing an ascii list of image filenames. For ascii
    list option, the filename containing the list must be provided
    preceded by a '@' character, e.g. 'images=@listoffiles.lis'.

*outfile*
    String. Name of the output file. The output file is a 14-column table.
    There is one row in the table for each time bin, and the output follows that from slotphot.
    The target and comparison flux values are the total signal divided by aperture area minus the corresponding background values. The flux ratio is the target flux/comparison flux.
    Errors assume photon noise for raw counts with standard error propagation for normalization and background subtraction.
    Normalized target, comparison, and background region errors are given by Sqrt(number of counts)/(aperture area).
    Background subtracted target and comparison errors are given by Sqrt(norm target error^2+norm comp error^2).
    Flux ratio error is given by the following: flux ratio*Sqrt((bkgdsub target error/bkgdsub target)^2+(bkgdsub comp error/bkgdsub comp)^2)
    The format is::

        id_number reference_time flux_ratio flux_ratio_error targetX target Y target_flux target_flux_error compX compY comp_flux comp_flux_error target_bkgd comp_bkgd

*outtype = [ascii|fits]*
    String. Output can be either plain ASCII (default) or a FITS table. While ASCII is generally easier to work with using other programs, FITS tables offer increased speed at readin and a convenient bundling of all relevant information such as object positions, binning and reference times used into it's header.

*srcfile*
        String. Name of the file defining the raw CCD pixel position of two sources
        in the image data and defining the size of circular aperture for the
        extraction of source counts and and a circular annulus for the
        extraction of background counts. The file format is an ascii table
        with 12 columns and 2 rows. An example is provided below::

            # star  x  y    x_end y_end  r  r_bkg1 r_bkg2 x1 y1 x2 y2
            # -------------------------------------------------------
            1     3  193  11    120    13 14     10     0  0  10 10
            2     3  37   12    133    13 14     20     20 20 50 40

        Rows beginning with the character '#' are ignored. Column 1 contains
        a numeric flag for the star which is either 1 or 2. The tool will
        follow any source drift image over time by linear interpolation between
        start and end position.
        The output file will always contain the source count ratio object 1 / object2, never the inverse ratio.
        Column 3 is the position of a source in (binned) CCD pixel
        units in the x-direction, or if you prefer the column number. Column
        4 refers to the corresponding y- or row coordinate. Both x and y can
        be routintely read off of your ds9 viewer or set using the bvitpreview tool.
        Columns 5 and 6 contain the end positions of the objects used for drift
        correction by linear interpolation.
        Column 7 contains the radius of the circular extraction region for the source counts. Columns 8 and 9 contain the inner and outer radii of the annulus used to extract background counts for each source. Any of the radial parameters are permitted to fall off of the edges of the image without breaking
        the tool but the tool does not take the image edges into account during calculations so this will influence the counts and it is reccomended to make sure that none of the regions goes outside the image.
        The last four columns specify the rectangular region used for background subtraction used if this mode is selected.

*imstart*
        String. Image used for selecting object start positions. The header time of this image is used for setting the reference time and for linear drift correction.

*imend*
        String. Image used for selecting object end position. The header time of this image is used for linear drift correction.

*tstep*
        Real. Size of time bin to use in seconds.

*maxdelay*
        Maximum delay between files in ms. If this delay is exceeded the header time of the current file will be used as a new reference.

*tgt_btype = [annulus|region]*
        Background type to use for target.

*cmp_btype = [annulus|region]*
        Background type to use for comparison.

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

This task does photometry optimized for the BVIT sparse data format. It reads a list of data files and a configuration file describing the positions of the source comparison and background regions (output from BVITPREVIEW) and produces a lightcurve.

Each BVIT fits table contains a list of up to 5 milion datapoints. A datapoint consists of the (x,y) position hit by a photon the counts produced in the phototube and it's arrival time in terms of a clock pulse count. The clock pulses are given off every 25 nanoseconds and the count is reset by a GPS signal every second.

This task first obtains a reference time from the header of the file specified in imstart. It subsequently loads the first fits file containing BVIT stream data and reads the data. The jumps corresponding to a GPS trigger are identified and data before the first jump is ignored. Subsequently the pulses are converted to time in seconds since the header time in imstart. This ensures correct handling of linear drift since the start image (which is used for position specification) was written.

The counts are then binned in time and photometry is performed. The mid-time of each bin is converted to seconds from noon on the previous day, starting with the reference time from the first data file. The background subtracted flux ratio is returned. This is (normalized target - normalized target background) / (normalized comparison - normalized comparison background). Normalization here refers to total counts in an aperture over a time bin divided by the area of the aperture. Flux for the target and the comparison (the numerator or denominator of the above equation, respectively), and normalized target and comparison background are returned. The positions of both target and comparison are corrected for linear drift, and these values are also provided.

Finally the outputs are returned in the file specified in outfile. They can then be viewed using BVITVIEW.


Examples
========

1. To get a lightcurve from a series of BVIT images specified in images.lst::

    --> bvitphot images="@images.lst" outfile="result.dat" outtype="ascii"
    srcfile="bvit.conf" imstart="OYCar_start" imend="OYCar_end"
    tstep=1.0 maxdelay=1000.0
    tgt_btype="annulus" cmp_btype="annulus"
    clobber=y logfile=salt.log verbose=y debug=n

Time requirements
=================

A Macbook with a 2.16 GHz Intel Core 2 Duo processor took about 10 seconds to process a 500 second dataset.

Bugs and limitations
====================

The current version of BVITPHOT has only been tested on one dataset.

Send feedback and bug reports to salthelp@saao.ac.za

See also
========

 :ref:`bvitpreview` :ref:`bvittofits`