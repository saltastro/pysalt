.. _saltfast:

********
saltfast
********


Name
====

saltfast -- Notify the PI that data were observed

Usage
=====

saltftp pinames scamobslog rssobslog datapath server username password cleanup
(clobber) logfile (verbose) (status)

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

*scamobslog*
    String. The name of a FITS table file containing SALTICAM output from
    the saltlog task. If the file resides in a separate directory then the
    absolute or relative path must be supplied with the file name. If no
    SALTICAM data exists from the night of interest, use scamobslog='None'.

*rssobslog*
    String. The name of a FITS table file containing RSS output from the
    saltlog task. If the file resides in a separate directory then the
    absolute or relative path must be supplied with the file name. If no
    RSS data exists from the night of interest, use rssobslog='None'.

*datapath*
    String. Absolute or relative path to the directory containing data
    organized by PI or OBSID.  Typically this path will be identical to
    the outpath argument given to the saltobsid task. Sub-directories must
    exist with names identical to the names listed in the pinames
    argument.

*server*
    String. The destination machine for the data to be transferred
    to. With the current version of the pipeline, transfer to the FTP
    machine is provided by server='www.saao.ac.za'.

*username*
    String. A valid username on the destination server.

*password*
    String. The valid password assocciated with username on the detination
    server.

*cleanup*
    Boolean. If cleanup='yes' the constructed tar files are deleted
    on the local machine after the ftp transfer. if cleanup='no', the
    local tar files are retained. The default value is cleanup='no'.

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
    Hidden Integer. Provided status=0 is passed to saltftp, a successful
    run of the task will yield status=0 at completion of the task.  Any
    other value for the status flag at completion indicates failed
    execution.

Description
===========

saltftp is called within the SALT pipeline after raw image files have
been cleaned and reduced (saltclean), separated into individual
programs (saltobsid) and documentation for each program has been
generated (salthtml). The task can also be executed standalone to
update existing data on the FTP server or process individual programs
rather than a full night of programs.

saltftp uses the RSS and SALTICAM observation logs created with
saltlog as the basis for it's tasks. In sequential order, saltftp will
read the pre-prepared observation logs and compile a list of PIs who
own data contained in the log. If pinames='all', data for every PI in
the observation log will have tar files transferred to the FTP
server. Otherwise pinames must contain a subset of the PIs in the
observation log. saltlog will ensure that the subset is
fully-contained in the observation log files or stop with an error.

The task will then ensure that each PI in the piname list has an
identically named sub-directory under the path defined by the datapath
argument. The matching is case-insensitve. If a mis-match is found
then the task will stop with an error. Each sub-directory will
normally contain data organized: datapath/PI/raw/ (raw data),
datapath/PI/product/ (reduced data) and datapath/PI/doc (observation
and reduction documentation).

For each matching sub-directory, all files contained within it will be
tarred and compressed using bzip2 compression into an archive called
YYYYMMDD.tar.bz2, where YYYY is the numerical year, MM the numerical
month and DD the numerical day at sunset of the night of
observations. The archive will be stored within the PI's
sub-directory, i.e. if the PI's surname is Dopey, the archive will be
stored in the file datapath/DOPEY/YYYYMMDD.tar.bz2.

Tar files are transferred to the public server using FTP. If there is
no directory on the FTP server with a name identical to the PI a new
one will be created. The tar file will be deposited in the PI-named
directory.

If cleanup='yes' the tar files created on the local machine will be
deleted before the task ends.

Examples
========


1. To create three bzipped2 tar files of PI-collated data called
BILL.tar.bz2, BEN.tar.bz2 and WEED.tar.bz2 in the existing
directories called /volumes/data2/BILL, /Volumes/data2/BEN and
/Volumes/data2/WEED and ftp the tar files to the server
ftp.xyz.ac.za, after which the tar files on the local machine will
be deleted::

    --> saltftp pinames='Bill,Ben,Weed' scamobslog='/Volumes/data1/Sobslog.fits'
    rssobslog='/Volumes/data1/Pobslog.fits' datapath='/Volumes/data2'
    server='ftp.xyz.ac.za' username='someone' password='!@#$%^&*'
    logfile='salt.log' verbose='yes'
    
    2. To create three bzipped2 tar files of PI-collated data called
    BILL.tar.bz2, BEN.tar.bz2 and WEED.tar.bz2 in the existing
    directories called /volumes/data2/BILL, /Volumes/data2/BEN and
    /Volumes/data2/WEED and ftp the tar files to the server
    ftp.xyz.ac.za, without deleting the tar files on the local machine::

        --> saltftp pinames='Bill,Ben,Weed' scamobslog='/Volumes/data1/Sobslog.fits'
        rssobslog='/Volumes/data1/Pobslog.fits' datapath='/Volumes/data2'
        server='ftp.xyz.ac.za' username='someone' password='!@#$%^&*'
        cleanup='no' logfile='salt.log' verbose='yes'


Time and disk requirements
==========================

saltftp creates large tar files, often containing many GB each. Space
for these files are required both on the working machine and the
detination server. Tarring GB of data may often require an hour or
more of processing.

Bugs and limitations
====================

Until the SALT database is released, data files will not contain an
observation ID within keywords. Before then, saltftp tar and transfer
data according to PI name.

Send feedback and bug reports to salthelp@saao.ac.za

See also
========

 :ref:`saltpipe`