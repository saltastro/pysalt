.. _saltpipe:

********
saltpipe
********


Name
====

saltpipe -- SALT data reduction pipeline

Usage
=====

saltpipe obsdate archive email ftp emserver emuser empasswd qcpcuser
qcpcpasswd ftpserver ftpuser ftppasswd median function order rej_lo
rej_hi niter interp (clobber) logfile verbose (status)

Parameters
==========


*obsdate*
    String. Observing date of the data. The pipeline expects all data in a
    sequence to have been obtained during a single night of observing. The
    convention in the FITS keywords is that the date is the sunset date at
    the start of a night. The format is YYYYMMDD, i.e. data taken on the
    night starting August 20, 2007 can be reduced using obdate='20070820'.
    saltpipe expects data to be organized according to SALT standard
    directory paths, e.g. /salt/scam/data/20070823/raw/ or
    /salt/data/2007/0823/scam/raw/. If raw data directories do not contain
    the file 'disk.file' then the pipeline will stop without reducing
    data. obsdate also defines the name of the temporary working directory
    which is created under the current directory.

*archive*
    Boolean. If archive='yes' saltpipe will copy raw and reduced files to
    the SALT data archive. If archive='no', raw and reduced files will
    remain in the working directory.

*ftp*
    Boolean. If ftp='yes' data will be tarred and bzipped into archive
    files, one for each PI. The archive files are copied to the public
    SALT FTP server for retrieval by the PI and the copy of the archive
    file on the pipeline machine is deleted.

*email*
    Boolean. If email='yes' saltpipe will email various status messages to
    the Cape Town support team during processing. Typical messages
    e.g. report that all data from the telescope has been downloaded or
    that the pipeline has completed with or without errors. The default
    email for these reports is sa@salt.ac.za. If email='yes' the arguments
    below: emserver, emuser and empasswd must have valid
    values. Furthermore, if email='yes', PIs will be emailed by the
    pipeline and notified of data availability on the FTP server. The
    email argument cannot be true if the ftp argument is false.

*emserver*
    String. The address of the email server.

*emuser*
    String. A valid username on the email server.

*empasswd*
    String. The valid password assocciated with the username on the email
    server.

*qcpcuser*
    String. A valid username for the qcpc.salt machine at the
    telescope. This machine contains ascii files consisting of any keyword
    changes required by the observer. It is the duty of Cape Town support
    to copy the keyword file over to the directory
    saltastro.saao:/home/sa/newheadfiles and name it
    list_newhead_YYYYMMDD. If the file is not found, saltpipe will attempt
    to ftp it directly from the qcpc.salt machine. Until the SALT database
    is ready, the keyword file's most important job is to couple data
    files with PIs. If no keyword file is found either on saltastro.saao
    or qcpc.salt then the pipeline will proceed but files will not be
    organized according to PI/owner.

*qcpcpasswd*
    String. The valid password assocciated with the username on qcpc.salt.

*ftpserver*
    String. The address of the public FTP server used as a temporary
    repository for PIs to wnload their data.

*ftpuser*
    String. A valid username on the FTP server.

*ftppasswd*
    String. The valid password assocciated with the username on the FTP
    server.

*median*
        Boolean. If median='yes' the columns in the overscan region will be
        median averaged before fitting a single function to characterize the
        row-dependent structure of the bias. If median='no' the overscan
        columns will be mean-averaged before fitting the function.

*function*
        String. The functional form of the fit intended to characterize the
        the bias structure in the overscan region. The user has variety of
        function options to choose from::

            chebyshev  - Chebyshev polynomial
            polynomial - standard polynomial
            legendre   - Legendre polynomial
            spline1    - linear splines
            spline3    - cubic splines

        If the chebyshev, legendre of spline functions are called then
        saltbias will use the IRAF task colbias to subtract the overscan bias
        from science frames and trim the images.

*order*
        Integer. The order of the polynomial, or number of spline knots, in
        the overscan function defined above.

*rej_lo*
        Float. The overscan fit is an iterative sigma-clipping procedure
        employed to remove the biasing effects of data outliers. After the
        first fit iteration, any data below the threshold of rej_lo (in units
        of the sigma deviation between data and fit) will be rejected and the
        fit re-performed. The iterations will continue until no more data
        points are rejected or the number of iteration exceeds the limit
        defined by the niter argument.

*rej_hi*
        Float.  After the first fit iteration, any data above the threshold of
        rej_hi (in units of the sigma deviation between data and fit) will be
        rejected and the fit re-performed. The iterations will continue until
        no more data points are rejected or the number of iteration exceeds
        the limit defined by the niter argument.

*niter*
        String. The maximum number of iteration to perform during the overscan
        sigma clipping procedure.

*interp*
        String. The interpolation scheme used to rebin pixel values during the
        process of mosaicing CCD amplifier images into one single detector
        image. This requires both translations and rotations. The choices are::

            linear  -- linear function
            nearest -- nearest pixel center
            poly3   -- 3rd order polynomial function
            poly5   -- 5th order polynomial function
            spline3 -- cubic spline function
            sinc    -- sinc function

        'nearest' is the least expensive for CPU processing but the least
        accurate, the 'sinc' function is the most expensive.

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
        Hidden integer. Provided status=0 is passed to saltpipe, a successful
        run of the task will yield status=0 at completion of the task.  Any
        other value for the status flag at completion indicates failed
        execution.

Description
===========

By nature, the saltpipe task is a linear series of individual pipeline
tasks. To understand the principles of the underlying sub-tasks, the
user is referred to the individual help documents for saltfixsec,
saltedtky, saltlog, saltclean, salthtml, saltobsid, saltftp and
saltemail. The user can run all of these steps individually using the
sub-tasks if required or re-perform parts of the pipeline
manually. The processing sequence is as follows:

1. Create a working direcrtory below the current directory with a name
defined by obsdate, e.g. YYYYMMDD.

2. Test for the existence of both RSS and SALTICAM data directories on
ctfileserver, consistent with obsdate. The absense of one
instrument directory does not results in an error, but the absence
of both does and the pipeline will stop, sending the user a warning
email.

3. Check that RSS and SALTICAM data directories are complete. If they
are not, the pipeline will stop, sending a warning email to the
user.

4. Notify the user by email that the pipeline has started (if
email='yes').

5. Copy all raw data to the temporary workspace.

6. Convert any SLOT mode data to FITS format using the task
saltbin2fits. FITS files are written to the raw directories in the
temporary workspace.

7. If data were obtained before 12/08/06, fix the erroneous SEC
keywords in raw data using the task saltfixsec. The procedure
overwrites raw files in the temporary workspace.

8. Identify a file containing header keyword changes requested by the
duty Astronomer at the telescope. These files will be of the name
/home/sa/newheadfiles/list_newhead_YYYYMMDD. If an appropriate file
does not exist, saltpipe will attempt to identify and FTP a
suitable file from the machine at the telescope qcpc.salt. If that
operation is unsuccessful, the pipeline will continue without
keyword edits. At the current phase of the project, the most
critical information that is missing from keywords is PI
names. Without these edits, data cannot be coupled to PIs by the
pipeline.

9. Perform header keyword edits using the task saltedtky. This
procedure overwrites raw files in the temporary workspace. FITS
tables containining a log of the keyword changes is stored in the
product directory in the temporary workspace.

10. Create observation logs for the RSS and SALTICAM data. These logs
are FITS tables stored in the product directory in the temporary
workspace. The logs contain one row for each data file and one
column for each primary header keyword. This step is probably
temporary and the functions of the log file in data processing
will be replaced by the SALT database once it is online.

11. Process and clean the raw images. The steps are keyword
preparation, gain correction, crosstalk correction, bias
subtraction and amplifier mosaicing. Tasks for flat fielding,
fringe subtraction and cosmic ray rejection are pending. All steps
are performed by calling the task saltclean, which internally
calls, saltprepare, saltgain, saltxtalk, saltbias, saltmosaic and
saltslot. Cleanded data are stored in the product directory of the
temporary workspace.

12. From the observation logs, identify all PIs associated with the
data. Collate data files, both raw and cleaned, according to PI,
create new directories on the temporry workspace named after the
PIs and populate them with symbolic links to the appropriate raw
and cleaned data files. All steps are performed by the task
saltobsid.

13. Generate HTML documentation containing details, of the night log,
exposure sequence, pipeline log etc. Each PI receives a copy of
the documentation which is deposited in the doc/ directory. This
procedure is performed by the task salthtml.

14. For each PI, archive the data in a bzipped tar file and tranfer it
to the public FTP server for retrieval by the PI. The tar file is
deleted from the temporary workspace directly after transfer to
the FTP server. The procedure is performed by the task saltftp.

15. Using saltemail, if email='yes', email each PI notfication that
they have data to retrive.

16. If archive='yes', move the temporary working directory to
ctfileserver for permanent storage.

17. If email='yes', email the user notification of the
successful/unsuccessful completion of the pipeline with procesing
statistics.

Examples
========

1. To execute the full pipeline on all observations from the night
starting 16 Aug 2006::

    --> saltpipe obsdate=20060816 archive=yes email=yes ftp=yes
    emserver=smtp.saao.ac.za emuser=**** empasswd=****
    qcpcuser=**** qcpcpasswd=**** ftpserver=www.saao.ac.za
    ftpuser=**** ftppasswd=**** median=n function=polynomial
    order=3 rej_lo=3.0 rej_hi=3.0 niter=10 interp=linear
    clobber=y logfile=saltpipe.log verbose=y
    
    The user and password arguments have been swapped for the string
    '****' in this example to avoid security issues.
    

Time and disk requirements
==========================

Individual unbinned full frame RSS image files can be 112MB in size. It is
recommended to use workstations with a minimum of 512MB RAM. Depending on
file numbers and sizes, the pipeline may take many hours to complete. The
goal is to keep the running times shorter than real time observing so that
the pipeline never becomes backlogged.

Bugs and limitations
====================

Ideally data should be reduced using data extracted from the SALT
database. Until the SALT database is released, data files will not
contain an observation ID within keywords. Before then, saltpipe will
operate upon PI names supplied by the observer.

See individual tools below for individual bugs and limitations.

Send feedback and bug reports to salthelp@saao.ac.za

See also
========

 :ref:`saltfixsec` :ref:`saltedtky` :ref:`saltlog` :ref:`saltclean` :ref:`salthtml` :ref:`saltobsid` :ref:`saltftp` :ref:`saltemail`