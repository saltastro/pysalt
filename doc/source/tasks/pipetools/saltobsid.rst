.. _saltobsid:

*********
saltobsid
*********


Name
====

saltobsid -- Collate data files from specific proposals

Usage
=====

saltobsid pinames obslog rawpath prodpath outpath (clobber) logfile verbose (status)

Parameters
==========


*pinames*
    String. List of Principle Investigator surnames. Surnames can be
    provided as a comma-separated list, or a string with a wildcard, or a
    foreign file containing an ascii list of surnames. For the ascii list
    option, the filename containing the list must be provided preceded by
    a '@' character, e.g. 'images=@listofnames.lis'. All surnames in the
    list must match surnames stored in the PROPOSERS FITS keyword of SALT
    image data, although string-matching within the task is
    case-insensitive. Alternatively if pinames='all' the task will collate
    data individually for all PI names that are tabulated within the
    observation log.

*obslog*
        String. The name of a FITS table file containing output from the
        saltlog task. If the file resides in a separate directory then the
        absolute or relative path must be supplied with the file name.

*rawpath*
        String. Absolute or relative path to the directory which contains the
        raw FITS files.

*prodpath*
        String. Absolute or relative path to the directory which contains the
        reduced FITS files.

*outpath*
        String. Absolute or relative path to the directory where new
        directories are created. Raw and reduced PI-dependent FITS files are
        linked symbolically from within the new diretories.

*(clobber)*
        Hidden boolean. If clobber='yes', files contained within the working
        directory will be overwritten by newly created files of the same
        name.

*logfile*
        String. Name of an ascii file for storing log and error messages
        written by the task. The file may be new, or messages can also be
        appended to a pre-existing file.

*(verbose)*
        Hidden Boolean. If verbose='no', log messages will be suppressed from
        both the terminal and the log file.  Error messages are excluded from
        this rule.

*(status)*
        Integer. Provided status=0 is passed to saltobsid, a successful run of
        the task will yield status=0 at completion of the task.  Any other
        value for the status flag at completion indicates failed execution.

Description
===========

saltobsid is called within the SALT pipeline after raw image files
have been cleaned and reduced (saltclean), and before data are copied
to the ftp server for PI retrieval (saltftp). The task simply
organizes the data according to the science program it was obtained
for. In the future, organization will occur based on the observation
ID of science images but currently programs are performed without
OBSIDs because the SALT database is not operational. Currently this
task bases its organization on the PI name stored in the primary HDU
keyword 'PROPOSER' of each science image.

saltobsid can also be run manually to extract all relavant data files
for a specific program from the archive. This is a useful feature when
a handful of files are required from a large directory of data.

While the data archiving format is not particularly rigid, saltiobid
does assume that data are archived accordng to the SALT pipeline model
in the sense that raw and reduced data files are found in separate
directories and that data from separate nights are stored in separate
directory trees. Typically the task will be called twice each day.
Once to collate the raw and reduced SALTICAM data and once again to
collate the raw and reduced RSS data.

As a prerequisite to running saltobid, the task saltlog must be run in
order to create a FITS-format observation log of the raw data
contained in a single directory. In terms of pipeline data it is
critical that the log contains raw data because there is not a
one-to-one correspondence between files in the raw directory and files
in the reduced products directory. saltobsid uses the raw data to
predict the contents of the product directory and will throw a flag if
a file is missing because it will mean an error has occured earlier in
the pipeline.

Within the outpath directory, saltobsid will create a directory for
each PI in the piname list. The name of the directory will be
identical to the content of the PROPOSER FITS keyword. Below these
directories, the task will write directories named 'raw' and 'product'
and then iterate through the list of files in the observation log,
writing symbolic links for each raw and associate reduced data file in
any appropriate PI directories. For data containing science data the
selection is trivial and based upon the value of the PROPOSER keyword.
For calibration exposures such as biases, flats and arcs, the selection
is based on various combinations of keywords listed here::

    CCDTYPE  -- type of exposure [OBJECT|ARC|BIAS|FLAT] (string)
    CCDSUM   -- on-chip binning of the instrument CCDs (string)
    GAINSET  -- CCD gain setting (string)
    ROSPEED  -- CCD reeadout speed (string)
    FILTER   -- instrument filter name (string)
    GRATING  -- grating name (string, RSS only)
    GR-ANGLE -- rotation angle of grating in deg (float, RSS only)
    AR-ANGLE -- articulation angle of spectrograph in deg (float, RSS only)

To couple a bias frame to science images in a program the CCDSUM,
GAINSET and ROSPEED keywords must match exactly. To couple a FLAT to a
science image, the CCDSUM, GAINSET, ROSPEED, FILTER and GRATING
keywords must match exactly and the GR-ANGLE and AR-ANGLE keywords
must match within a mechanical tolerance of plus or minus 0.01
deg. Arc exposures must meet the same criteria as the flats above
except the filter matching is dropped.

Generally, object exposures will belong to a single program,
calibration exposures can belong to multiple programs and so multiple
symbolic links will potentially be created to a single calibration
data file. There are exceptions to this rule, e.g. a standard star
observation may belong to several programs. In these case SALT
Astronomers generally reference all PIs in the FITS keyword,
separated by commas.  saltobsid recognises this format and will create
symbolic links for all PIs listed.

saltobsid currently pays no heed to the PROPOSER keyword in
calibration data. Calibration file coupling is performed solely on the
basis of the other keywords in the list above.

If relative paths to the input and output directories are provided by
the user, saltobsid will convert them to absolute paths before
creating symbolic links to data files.

In the PI directories, SALTICAM and RSS data files are not written to
separate directories. The only discrimator at the sub-PI directory
level is whether data is raw or reduced.

The use of the clobber argument should be used with care. If clobber=y
old symbolic links will be overwritten if a new one is being
created. But this is not the same procecdure as deleting all symbolic
links from a directory before creating new ones.

Because the pipeline orders and stores data by instrument name, if a
program contains data from multiple instruments then saltobsid must be
run multiple times, once per instrument, to create a full complement
of symbolic links.

Standard star observations are given no special keyword flags by the
telescope software. They are treated as normal object images in
starobsid and it is the users perogative to make sure that the
PROPOSER keyword references the correct PI in standard star data.

Examples
========

1. To create three new directories, ../data4/Bill, ../data4/Ben and
../data4/Weed, asscoiated with three PIs that are referenced in the
observation log /Volumes/data1/obslog.fits. Raw image data stored
in /Volumes/data2 and reduced data stored in /Volumes/data3 will be
symbolically linked within the outpath directories::

    --> saltobsid pinames='Bill,Ben,Weed' obslog='/Volumes/data1/obslog.fits'
    rawpath='Volumes/data2' prodpath='/Volumes/data3' outpath='../data4'
    clobber='yes' logfile='salt.log' verbose='yes'
    

Time and disk requirements
==========================

saltobsid creates only directories and symbolic links. Therefore there
are no significant, CPU, memory or disk capacity constraints.

Bugs and limitations
====================

Until the SALT database is released, data files will not contain an
observation ID within keywords. Before then, saltobsid will filter and
organize data according to PI name.

Send feedback and bug reports to salthelp@saao.ac.za

See also
========

 :ref:`saltpipe` :ref:`saltclean` :ref:`saltftp`