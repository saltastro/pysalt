.. _saltemail:

*********
saltemail
*********


Name
====

saltemail -- email adapted standard letter to PIs

Usage
=====

saltemail pinames scamobslog rssobslog emailfile readme logfile (verbose) (status)

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

*emailfile*
    String. The name of an ascii file containing PI names and their email
    addresses. If the file is contained in a different directory then the
    relative or absolute path must also be included. The format of the
    file is a two-column table. Th first column is an exact match for PI
    names contained in data keywords (although case is not important). The
    second column is hte email address and the two columns are separated
    by a comma. Below is an example::

        bloke,abloke@saao.ac.za
        geezer,ugeezer@rutger.edu
        doofus,doof@soton.ac.uk

*readme*
    String. An ascii file containing relavant information for the data
    PI. This information is converted to a PI-tailored email message.  A
    canned version of the file is provided in
    /iraf/extern/salt/pysalt/html/readme.template which can be copied and
    edited by users. If the readme file does not exist, the task will
    continue after providing a warning.

*server*
    String. The address of the email server.

*username*
    String. A valid username on the email server.

*password*
    String. The valid password assocciated with username on the email
    server.

*logfile*
    String. Name of an ascii file for storing log and error messages
    written by the task. The file may be new, or messages can also be
    appended to a pre-existing file.

*(verbose)*
    Hidden Boolean. If verbose='no', log messages will be suppressed from
    both the terminal and the log file.  Error messages are excluded from
    this rule.

*(status)*
    Hidden Integer. Provided status=0 is passed to saltemail, a successful
    run of the task will yield status=0 at completion of the task.  Any
    other value for the status flag at completion indicates failed
    execution.

Description
===========

saltemail is the final task called within the SALT automated
pipeline. After data has been reduced, archived and transferred to the
public FTP server, saltemail notifies PIs of the availability of their
data by email.

Before execution, a pre-requisite is that the email addresses of the
PIs are recorded in the table referred to by the emailfile
argument. It is the repsonsibility of the SALT Astronomers to ensure
that the table remains up-to-date. Since the file is based upon PI
name there is potential for name duplication and email address
ambiguity.

The PI names referenced in the piname input argument are
cross-correlated with the PI names recoded in the RSS and SALTICAM
observation logs. If there are mis-matches then the task will stop
with an error message. The task will also stop if there are mutliple
matches of the same name in the email table.

The email message itself is copied from a text template, pointed to
with the readme argument. The template is copied exactly into the body
of the message with three exceptions. Any instance of the string 'SALT
user' will be replaced by the PI's name, any instance of the string
'yourname' will be replaced by the PI's name and any instance of the
string 'yyyymmdd' will be replaced by the date of the observation.

It is exected that saltemail will be run manually for a period of
time. Once confidence has been gained in the performance of the
pipeline, saltemail will be run antomatically.

Examples
========

1. To send email to three PIs based upon standard text kept in the
template file /iraf/extern/salt/pysalt/html/readme.template::

    --> saltemail pinames='Bloke,Geezer,Doofus'
    scamobslog='/Volumes/data1/Sobslog.fits'
    rssobslog='/Volumes/data1/Pobslog.fits'
    emailfile='/Volumes/data2/pi_emails.dat'
    readme='/iraf/extern/salt/pysalt/html/readme.template'
    logfile='salt.log' verbose='yes'

Time and disk requirements
==========================

saltemail requires no significant disk space. Execution time is
typically 1 sec.

Bugs and limitations
====================

Ideally email addresses should be extracted from the SALT
database. Until the SALT database is released, data files will not
contain an observation ID within keywords. Before then, saltemail will
operate upon PI names.

Send feedback and bug reports to salthelp@saao.ac.za

See also
========

 :ref:`saltpipe`