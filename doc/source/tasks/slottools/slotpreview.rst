.. _slotpreview:

***********
slotpreview
***********


Name
====

slotpreview -- Preview SALT slotmode data and select object and background regions

Usage
=====

slotpreview images outfile ignorexp (tgt_col) (cmp_col) (tgt_lw)
(cmp_lw) (clobber) (logfile) (verbose) (debug)

Parameters
==========


*images*
    String. List of images to display.

*outfile*
    String. Name of the output file. The output file is a 11-column ascii table
    with the object positions, radii and selected background regions. This is used as input for slotphot.

*ampperccd*
    Integer. Number of amplifiers used for the ccd. This will be used in combination with ignorexp to select the first exposure to load and to calculate the amplifier containing the target and/or comparison. This information is passed to slotphot as the second column of outfile.

*ignorexp*
    Integer. Number of exposures to ignore. Usually the first couple of exposures are filled with nonsense data. The first exposure loaded is ignorexp*amperccd+1.

*recenter_radius*
    Integer. Radius to use for automatic recentering algorithm. The algorithm, which returns the centroid of the image subsection, generally works best if this radius is not to large.

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

The output is an ascii table with 11 columns and 2 rows. An example is provided below::

    # object amp x    y    r   r1 r2 x1 y1 x2 y2
    # ------------------------------------------
    1      3   193  11   12  13 14 0  0  0  0
    2      2   37   12   13  0  0  20 20 40 50

Column 1 contains a numeric flag for the star which is either 1 or 2.
Column 2 is the amplifier containing the object, for slotmode this is a number between 1 and 4.
Column 3 is the position of a source in (binned) CCD pixel units in the x-direction, or if you prefer the column number. Column 4 refers to the corresponding y- or row coordinate.
Column 4 contains the radius of the circular extraction region for the source counts.
Columns 5 and 6 contain the inner and outer radii of the annulus used to extract background counts for each source.
The last four columns specify the rectangular region used for background subtraction used if this mode is selected.





    --> slotpreview images = "@images.lst" outfile = "slot.conf"
    ignorexp = 6 tgt_col = "g" cmp_col = "r" tgt_lw = 2
    cmp_lw = 2 clobber=y logfile=salt.log verbose=y

First click next until you find the exposure that contains the target.
Then click on the capture button next to target position and click on the target.
The marker is automatically recentered on the centroid of the object, if this is not desired deselect the recenter checkbox.
Click on the target radius position capture button and click next to the target to the point where you want the target region to extend.
Now click the target background tab and select a background type.
If annulus is selected click on the annulus capture button and click twice next to the target in the image.
The first click selects the inner radius of the target background annulus and the second click the outer radius.
If the region type is selected click on the region capture button and also click on the image twice to select the opposite corners of the rectangle.

Repeat this procedure for the comparison object.

In order to fine tune the parameters just repeat the desired step (in any order) or type the values in manually.

Optionally you may enable the lock button (broken chain icon) to ensure the target and comparison (or their background regions) always have the same size.
Because of strong vignetting in the slotmode data it may be usefull to center the background region vertically, to do this automatically enable the centering button on the background selection tab.
This changes the capture behaviour slightly, instead of specifying the vertices directly the distance of the cursor to the center line of the image is used on both sides of the line to set the y coordinate of the vertices the behaviour in the x direction is left unchanged.

Finally click save to save the configuration file and quit the task.

If anything goes wrong just click cancel to ignore changes and quit.

Time and disk requirements
==========================

Should load in maximum of 5 seconds on any recent machine.
Main time consumption is in loading the fits images.


Bugs and limitations
====================

Virtually untested, please help here.

Send feedback and bug reports to salthelp@saao.ac.za

See also
========

 :ref:`slotphot`