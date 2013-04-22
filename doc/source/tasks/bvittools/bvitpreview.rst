.. _bvitpreview:

***********
bvitpreview
***********


Name
====

bvitpreview -- Preview BVIT data and select object and background regions

Usage
=====

bvitpreview imstart imend outfile (tgt_col) (cmp_col) (tgt_lw) (cmp_lw)
(cmap) (scale) (contrast)
(clobber) (logfile) (verbose) (debug)

Parameters
==========


*imstart*
    String. Image used to select object start positions.

*imend*
    String. Image used to select object end positions.

*outfile*
    String. Name of the output file. The output file is a 12-column ascii table
    with the object positions, radii and selected background regions. This is used as input for bvitphot.

*tgt_col = [b|g|r|c|m|y|k|w]*
    String. Line color used for target markers. This can be any valid Matplotlib color.

*cmp_col = [b|g|r|c|m|y|k|w]*
    String. Line color used for comparison markers. This can be any valid Matplotlib color.

*tgt_lw*
    Integer. Line width used for target markers.

*cmp_lw*
    Integer. Line width used for comparison markers.

*cmap = [autumn|bone|cool|copper|flag|gray|hot|hsv|jet|pink|prism|spring|summer|winter|spectral]*
    String. Colormap used for displaying image. This can be any valid Matplotlib colormap name.

*scale = [minmax|zscale]*
    String. Algorithm used to set dynamic range of image display, in general zscale gives the best results.

*contrast*
    Real. Contrast to use with zscale.

*clobber*
    Hidden Boolean. If clobber=y the tool is permitted to overwrite an existing
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


This tool is used for selecting the target and comparison positions and their and radii as well as their corresponding background regions.

All position parameters may be entered manually or captured by clicking on the image.

The output is an ascii table with 12 columns and 2 rows. An example is provided below::

    # star  x  y    x_end y_end  r  r_bkg1 r_bkg2 x1 y1 x2 y2
    # -------------------------------------------------------
    1     3  193  11    120    13 14     10     0  0  10 10
    2     3  37   12    133    13 14     20     20 20 50 40

Column 1 contains a numeric flag for the star which is either 1 or 2.
Column 3 is the position of a source in (binned) CCD pixel units in the x-direction, or if you prefer the column number. Column 4 refers to the corresponding y- or row coordinate.
Columns 5 and 6 contain the end positions of the objects used for drift correction by linear interpolation.
Column 7 contains the radius of the circular extraction region for the source counts.
Columns 8 and 9 contain the inner and outer radii of the annulus used to extract background counts for each source.
The last four columns specify the rectangular region used for background subtraction used if this mode is selected.





    --> bvitpreview imstart = "OYCar_start" imend = "OYCar_end"
    outfile = "bvit.conf" tgt_col = "g" cmp_col = "r"
    tgt_lw = 2 cmp_lw = 2 cmap="jet" scale="zscale"
    contrast=0.1
    clobber=y logfile=salt.log verbose=y debug=no

Then click on the capture button next to target position and click on the target.
Click on the target radius position capture button and click next to the target to the point where you want the target region to extend.
Click on the target end position capture button and click on the target in the second image.
Now click the target background tab and select a background type.
If annulus is selected click on the annulus capture button and click twice next to the target in the image.
The first click selects the inner radius of the target background annulus and the second click the outer radius.
If the region type is selected click on the region capture button and also click on the image twice to select the opposite corners of the rectangle.

Repeat this procedure for the comparison object.

In order to fine tune the parameters just repeat the desired step (in any order) or type the values in manually.

Optionally you may enable the lock button (broken chain icon) to ensure the target and comparison (or their background regions) always have the same size.

Finally click save to save the configuration file and quit the task.

If anything goes wrong just click cancel to ignore changes and quit.

Time and disk requirements
==========================

Should load in maximum of 5 seconds on any recent machine.
Main time consumption is in loading the fits images.


Bugs and limitations
====================

None

Send feedback and bug reports to salthelp@saao.ac.za

See also
========

 :ref:`bvitphot`