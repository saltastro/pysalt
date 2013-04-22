.. _salthtml:

********
salthtml
********


Name
====

salthtml -- Collate data files from specific proposals

Usage
=====

salthtml scamobslog rssobslog htmlpath nightlog notes (clobber)
logfile verbose (status)

Parameters
==========


*scamobslog*
    String. The name of a FITS table file containing SALTICAM output from
    the saltlog task. If the file resides in a separate directory then the
    absolute or relative path must be supplied with the file name. If no
    SALTICAM data exists from the night of interest, use
    scamobslog='None'.

*rssobslog*
    String. The name of a FITS table file containing RSS output from the
    saltlog task. If the file resides in a separate directory then the
    absolute or relative path must be supplied with the file name. If no
    RSS data exists from the night of interest, use rssobslog='None'.

*htmlpath*
    String. Absolute or relative path to the directory which will contain
    the html output. Typically this path will be identical to the outpath
    argument given to the salthtml task. If this so, the output will be
    copied to the multiple PI directories and tailored individually.

*nightlog*
    String. An ascii file containing a log of the night's activities,
    handwritten by the duty astronomer. In the future, this log will be
    copied from the SAMMI machine in the SALT control room, but is
    currently emailed across the operations team by the astronomer each
    morning. It is the duty of Cape Town personnel to extract the night
    log email into an ascii file before the pipeline starts. If the night
    log does not exist, the task will continue after providing a warning.

*readme*
    String. An ascii file containing relavant information for the data
    PI. A canned version of file is provided in
    /iraf/extern/salt/pysalt/html/readme.template which can be copied and
    edited by users. If the readme file does not exist, the task will
    continue after providing a warning.

*(clobber)*
    Hidden boolean. If clobber='yes', files contained within the working
    directory will be overwritten by newly created files of the same
    name.

*logfile*
    String. Name of an ascii file for storing log and error messages
    written by the task. The file may be new, or messages can also be
    appended to a pre-existing file. Generally this file should be the
    same as the log file for the entire pipeline because it is appended
    to the html documentation as reference material.

*(verbose)*
    Hidden Boolean. If verbose='no', log messages will be suppressed from
    both the terminal and the log file.  Error messages are excluded from
    this rule.

*(status)*
    Integer. Provided status=0 is passed to salthtml, a successful run of
    the task will yield status=0 at completion of the task.  Any other
    value for the status flag at completion indicates failed execution.

Description
===========

salthtml creates html pages documenting the observing and pipeline
procedures undergone to obtain and reduce a series of data. The basis
of these pages are::

    1. The night log maintained by the duty astronomer.
    2. The pipeline log created during automated data cleaning.
    3. Instrument diagnostics contained in individual image file data.
    4. Telescope, dome and atmosphere environment data.
    5. An exposure log derived from image data keywords.

The current html documents are for illustrative purposes only. The
SALT project has yet to decide upon requirements from the SALT
pipeline. Therefore this help file will not list the output in detail.

If the htmlpath also contains subdirectories named after each PI (see
the saltobsid task) then salthtml will copy the html documentation to
each subdirectory and tailor them individually.

Examples
========

1. To create a directory named /Volumes/data3/doc, containing html
documentation for all observations from the night starting 16 Aug
2006::

    --> salthtml scamobslog='/Volumes/data1/S20060816OBSLOG.fits'
    rssobslog='/Volumes/data1/S20060816OBSLOG.fits'
    htmlpath='/Volumes/data3' nightlog='night_report_20060810'
    notes='/iraf/extern/salt/pysalt/html/readme.template'
    clobber='yes' logfile='saltlog' verbose='yes'

Time and disk requirements
==========================

salthtml creates only directories and small html documents. Therefore
there are no significant, CPU, memory or disk capacity constraints.

Bugs and limitations
====================

Parsing night logs and edited readme files is not currently an
automatic procedure.  The content of the html output is not currently
under document control. They are examples constructed by the
programmer. Formal decisions concerning the content of the pipeline
output are pending.

Send feedback and bug reports to salthelp@saao.ac.za

See also
========

 :ref:`saltpipe`