**********
Installing
**********
Instructions for the installation of the PySALT package are below.
In order to install PySALT, several other packages have to be installed as well including python, PyRAF, and IRAF.
Descriptions and links to each of those packages are included under dependencies.
All of the dependent programs are available on the SALT `ftp site <ftp://www.salt.ac.za/pub/>`_.
If you are downloading the package from outside South Africa, it is suggested that the dependencies be downloaded from their original source.

Prerequisites
=============
* This package requires `IRAF <http://iraf.noao.edu/>`_ version 2.12.2 or later.
* It also requires `PyRAF <http://www.stsci.edu/resources/software_hardware/pyraf>`_ version 1.7.1 or later.
    * tcl: Version 8.3.5 or higher
    * tk:  v8.3.5 or higher
    * qt: Version 4.0 or later
    * python: v2.5.1 or higher
    * readline: v4.3
    * Pmw: v1.2
    * libf2c

* Additional required python packages:
    * `matplotlib <http://matplotlib.sourceforge.net/>`_: v0.91 or later
    * `pyfits <http://www.stsci.edu/resources/software_hardware/pyfits>`_: v1.3
    * `numpy <http://numpy.scipy.org/>`_: v1.04 or later
    * `scipy <http://www.scipy.org/>`_: v0.5.2 or later
    * `pyqt4 <http://www.riverbankcomputing.co.uk/software/pyqt/download>`_: v4.6 or later
    * `sphinx <http://sphinx.pocoo.org/>`_ (optional): to compile the documentation

Downloading
===========
The PySALT user package v0.40 was released on 01 February 2010 and can be downloaded `here <http://www.salt.ac.za/~crawford/pysalt/pysalt.v0.40.tar.gz>`_.
`Release notes <http://www.salt.ac.za/science-support/salt-data-reduction/pysalt-users-package/release-notes/>`_ are available.

Installation of the PySALT package
==================================
To install the package as root, create a directory to contain the PySALT external package files.
This directory should be outside the IRAF directory tree and must be owned by the IRAF account.
In the following example, this root directory is named::

  /iraf/extern/pysalt

Make the appropriate file name substitutions for your site.

#. Log in as IRAF and edit the *extern.pkg* file in the *hlib*
   directory to define the package to the CL.
   From the IRAF account, outside the CL, you can move to this directory with the commands::

     % cd /iraf/iraf/unix/hlib/

   Define the environment variable salt to be the pathname to the salt
   root directory. UNIX pathnames must be terminated with a ``/``.
   Edit *extern.pkg* to include::

     reset pysalt     = /iraf/extern/pysalt/
     task  pysalt.pkg = pysalt$pysalt.cl

   Near the end of the *hlib$extern.pkg* file, update the definition of *helpdb*
   so it includes the salt help database, copying the syntax already used
   in the string. Add this line before the line containing a closing quote::
        
     ,pysalt$lib/helpdb.mip\

#. Change directories to the PySALT root directory created above
   and unpack the download file *pysalt.tar.gz*::

     % cd /iraf/extern/pysalt/
     % tar zxvf <path>/pysalt.tar.gz

   where ``<path>`` is the relative path from ``/iraf/extern/pysalt`` to the downloaded file.

#. To install the package under your home directory, unpack the
   downloaded file *pysalt.tar.gz*.
   Then add the following line to your user *login.cl*::
    
     reset pysalt   = "[your path]"
     task  pysalt.pkg   = "pysalt$pysalt.cl"
     reset helpdb    = (envget("helpdb") // ",pysalt$lib/helpdb.mip")

#. Build the accelerated and GUI modules::

     % cd /iraf/extern/pysalt/
     % make

#. Finally to build the documentation
   make sure you have *sphinx* installed and execute::

     % cd /iraf/extern/pysalt/doc
     % make html

   The start point for the documentation is in ``build/index.html``.
   For a manual in *pdf* format (you need LaTeX as well) use::

     % cd /iraf/extern/pysalt/doc
     % make latex
     % cd /iraf/extern/pysalt/doc/build/latex
     % pdflatex PySALT.tex
     % pdflatex PySALT.tex

The manual will than be in ``PySALT.pdf``

