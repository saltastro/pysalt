.. _saltedtky:

*********
saltedtky
*********


Name
====

saltedtky -- Edit primary header keywords in multiple files

Usage
=====

saltedtky inpath outpath keyfile (record) (recfile) clobber logfile (verbose)

Parameters
==========


*inpath*
    String. Path to the directory containing the files which are to be edited.
    The path can be the full path or the relative path. The files within the
    directory must all have been taken during the same night of observing.
    The tool will exit with a warning if this is not the case. Also the
    file names must also conform to the SALT standard, e.g. P200701080001.fits,
    with any number of additional character prefixes at the *beginning* of
    the name.

*outpath*
    String. Path to the directory which will contain new files with edited
    or additional primary keywords. If the output directory does not exist
    the tool will attempt to create it, but this will only be successful
    under the usual constraints of unix directory creation, e.g. the
    parent directory must already exist. outpath and inpath can have the
    same value, but only if clobber=y. This will allow the tool to
    overwrite the input files. Users should decide for themselves whether
    this is a good idea.  Each new file will contain the date it was
    created by the task in the SAL_EDT keyword and the arguments used to
    create the file in the HISTORY keywords.

*keyfile*
    String. The name (with path if necessary) of an ascii file containing
    a list of the files that require editing, the keywords that require
    adding or adapting and the new value of the keywords. The structure of the
    file will look like the following::

        P1-2     NOTES="junk"
        S1       PROPOSAL="Kniazev" MASKID="1.0" NOTES="B image"
        S2       PROPOSAL="Kniazev" MASKID="1.0" NOTES="slit view"
        P3-4     PROPOSAL="Kniazev" MASKID="1.0" NOTES="PA=+100"
        P5       PROPOSAL="Kniazev" MASKID="1.0" NOTES="CuAr"
        P6-10    PROPOSAL="Kniazev" MASKID="1.0" NOTES="QTH1"
        P11      PROPOSAL="Kniazev" MASKID="1.0"
        P12      NOTES="light leak?"
        P17      PROPOSAL="Kniazev,Albrow" MASKID="1.0" NOTES="standard"
        P18      PROPOSAL="Kniazev" MASKID="1.0"
        P19	 NOTES="junk"
        S3       PROPOSAL="Albrow" NOTES="slit view"
        P20-22   PROPOSAL="Albrow" MASKID="1.0" NOTES="PA=+6"
        P23      PROPOSAL="Albrow" MASKID="1.0"
        P24-26   PROPOSAL="Albrow" MASKID="1.0" NOTES="PA=+6"
        P27      PROPOSAL="Albrow" MASKID="1.0"
        P28-32   PROPOSAL="Albrow" MASKID="1.0" NOTES="QTH1"
        P38-40   NOTES="junk"
        S4       PROPOSAL="Hughes" NOTES="slit view"
        P41-43   PROPOSAL="Hughes" MASKID="2.0" NOTES="PA=+93.5"
        P44      PROPOSAL="Hughes" MASKID="2.0"
        P45-49   PROPOSAL="Hughes" MASKID="2.0" NOTES="QTH1"
        S5       PROPOSAL="Nordsieck" MASKID="0.6" NOTES="B image"
        S6       PROPOSAL="Nordsieck" MASKID="0.6" NOTES="slit view"
        P50      PROPOSAL="Nordsieck" MASKID="0.6"
        P51-53   PROPOSAL="Nordsieck" MASKID="0.6" NOTES="HD37903"
        P54-55   PROPOSAL="Nordsieck" MASKID="0.6" NOTES="NGC 2023"
        P56      PROPOSAL="Nordsieck" MASKID="0.6" NOTES="NGC 2023, lost interf."
        P57      PROPOSAL="Nordsieck" MASKID="0.6"
        P58-62   PROPOSAL="Nordsieck" MASKID="0.6"
        S6       PROPOSAL="Hughes" NOTES="slit view"
        P68-69   PROPOSAL="Hughes" MASKID="2.0" NOTES="twilight starts 04:13"
        P70      PROPOSAL="Hughes" MASKID="2.0"

    Currently, the main role of this tool is within the data pipeline, in
    adding vital keywords to data that are missing, such as the PIs
    identity and size of the slit-mask, and to add useful historical notes
    to the files.  Typically the ascii file is created by the SALT
    Astronomer as he/she observes at the telescope before it is
    transferred and stored in Cape Town. The first column refers to the
    instrument (S = SALTICAM, P = RSS) and the data file numbers, i.e. the
    last four characters in the file name before the .fits extension. The
    format of file names is typically IYYYYMMDDnnnn.fits.  I is the
    instrument - either S or P, YYYYMMDD is the date of the observation
    and nnnn refers to the chronological data number. In the keyfile
    format P78-102 is an acceptable abbreviation of P0078-P0102.  The
    hypen refers to a continuous range of files. Occasionally this scheme
    has to be extended to five characters when > 10,000 files are created
    over a 24 hour period. The PROPOSAL keyword is now redundant. It will
    continue to be written by saltedtky for a short while, but the correct
    keywords of PROPID and PROPOSER are also added. Data currently do not
    have PROPID numeric codes propagated through the telescope system, so
    the name of the PROPOSER is temporarily used as a dummy value instead.

*(record)*
    Boolean. If record=y, a new FITS table will be constructed which
    records a log of the changes.

*(recfile)*
    String. The path and name of the FITS table to record keyword changes. The
    file will contain an empty primary header wit the table populating the
    first extension, called 'NEWKEYS'. The table will contain one row for
    each keyword change performed by the task. There will be four table
    columns::

        FILE - the file in which the change has occured
        KEYWD - the name of the keyword which has changed
        OLD_VAL - the old value of the keyword
        NEW_VAL - the new value of the keyword

    If a new keyword is added to the ouput file then the OLD_VAL column
    will contain the string 'None'.  The SAL_EDT keyword written to the
    table extension contains the time and date at which the file was
    created. The HISTORY keyword contains the arguments that were used to
    create the file.

*clobber*
    Boolean. If set to 'yes' files contained within the outpath
    directory will be overwritten by newly created files of the same
    name. This is a required option when the inpath and outpath arguments
    point to the same directory.

*logfile*
    String. Name of an ascii file for storing log and error messages
    from the tool. The file may be new, or messages can also be appended to a
    pre-existing file.

*(verbose)*
    Hidden Boolean. If verbose=n, log messages will be suppressed.

Description
===========

saltedtky was created to perform the function of adding important
primary extension keywords to science files passing through the data
pipeline. This task will be required until the telescope systems
mature to the state where all required keywords are propagated
correctly through the telescope. Currently major missing keywords are
any information coupling a data file to the the Proposer, and
information pertaining to the nature of the RSS slit mask. However
this tool may also find some role in keyword writing outside of the
pipeline because it provides a relatively simple method of adding or
changing keyword data over a large number of file simultaneously.


Examples
========

1. To add or change primary keywords in science files residing
in the /Volumes/data1/ directory and creating new files in the
/Volumes/data2 directory::

    --> saltedtky inpath='/Volumes/data1' outpath='/Volumes/data2'
    keyfile='night20061012.dat' clobber='n' logfile='salt.log'
    verbose='yes'

2. To add or change primary keywords in science files residing
in the /Volumes/data1/ and overwrite the old files::

    --> saltedtky inpath='/Volumes/data1' outpath='/Volumes/data1'
    keyfile='night20061012.dat' clobber='y' logfile='salt.log'
    verbose='yes'

Time and disk requirements
==========================

Unbinned raw full-frame RSS images are 111MB. It is recommended
to use workstations with a minimum of 512MB RAM. On a contemporary linux
machine, one file can be processed in a few seconds.

Bugs and limitations
====================

saltedtky does not have the flexibility to update files obtained on
multiple days. It is the users responsibility to ensure that the ascii
file containing keyword updates is formatted correctly.

The task will only perform keyword editing on files with names
==============================================================

consistent with the raw naming structure, e.g., IYYYYMMDDNNNN.fits.

The record file is created anew each time the task is executed.
===============================================================


The task has no functionality to append to old record files.
============================================================


Send feedback and bug reports to salthelp@saao.ac.za
====================================================



See also
========

 :ref:`saltpipe`