###########
Development
###########
Information for PySALT developers.

******************
Coding conventions
******************
For Python, `PEP 8 <http://www.python.org/dev/peps/pep-0008/>`_ has emerged as the style guide that most projects adhere to; it promotes a very readable and eye-pleasing coding style. All PySALT modules should follow these conventions. The most important points are: 

* Use 4-space indentation, and no tabs.
  4 spaces are a good compromise between small indentation (allows greater nesting depth) and large indentation (easier to read).
  Tabs introduce confusion, and are best left out.
* Wrap lines so that they don’t exceed 79 characters.
  This helps users with small displays and makes it possible to have several code files side-by-side on larger displays.
* Use blank lines to separate functions and classes, and larger blocks of code inside functions.
* When possible, put comments on a line of their own.
* Use docstrings.
* Use spaces around operators and after commas, but not directly inside bracketing constructs: ``a = f(1, 2) + g(3, 4)``.
* Name your classes and functions consistently; the convention is to use ``CamelCase`` for classes and ``lower_case_with_underscores`` for functions and methods.
  Always use ``self`` as the name for the first method argument.
* Don’t use fancy encodings if your code is meant to be used in international environments.
  Plain ASCII works best in any case.

*********************
Documenting your code
*********************

PySALT uses `Sphinx <http://sphinx.pocoo.org/index.html>`_ for automatic documentation generation.
In order to do this please follow these conventions.

* Every module, class, flunction and method needs to have a *docstring* describing it. At the minimum this should be a simple plain text string. For example::

    """This is a sample module."""

    class Test:
        """This class generates test objects.

        Some more documentation.
        """
        def __init__(self):
            """some initialization."""
            pass

    def tester():
        """This function prints test."""
        print test

* For better structure use `reStructuredText <http://sphinx.pocoo.org/rest.html>`_ in your docstrings. This has a simple *wiki* like syntax for:
    * extensive `markup <http://sphinx.pocoo.org/markup/inline.html>`_ such as sectioning, bullet lists, numbered lists and descriptions;
    * source code examples, `notes <http://sphinx.pocoo.org/markup/para.html#dir-note>`_, `warnings <http://sphinx.pocoo.org/markup/para.html#dir-warning>`_, `todo <http://sphinx.pocoo.org/ext/todo.html>`_ items;
    * including `images <http://sphinx.pocoo.org/rest.html#images>`_ and LaTeX `formatted math <http://sphinx.pocoo.org/ext/math.html?highlight=math#module-sphinx.ext.mathbase>`_;
    * cross `referencing <http://sphinx.pocoo.org/markup/inline.html>`_, modules, classes, functions and methods.

****************************
Documenting IRAF/PyRAF tasks
****************************

IRAF uses it's own documentation format for documenting tasks.
Every PySALT, PyRAF task should have a help page *taskname.hlp* in the ``doc`` directory of it's containing package.
The :ref:`tasks` help page in this documentation is automatically generated from these ``.hlp`` files.
In order to ensure a proper conversion please use the following template::

    # South African Astronomical Observatory
    # PO Box 9
    # Observatory 7935
    # South Africa
    # ******************************************************************
    .help taskname date pysalt.container_package
    .ih
    NAME
    taskname -- A short description
    .ih
    USAGE
    "taskname(argument1, argument2, argument3)"
    .ih
    PARAMETERS
    .ls argument1
    type (e.g. String.). Description of parameter.
    .le
    .ls argument2
    String. Description of parameter.
    .le
    .ls argument3
    String. Description of parameter.
    .le
    .ih
    DESCRIPTION
    A long description of the program.
    .ih
    EXAMPLES
    1. Example
    .nf
        some source code
        some more source code
    .fi
    .ih
    TIME REQUIREMENTS
    Estimate of load on system.
    .ih
    BUGS AND LIMITATIONS
    Known limitations
    
    Send feedback and bug reports to salthelp@saao.ac.za
    .ih
    SEE ALSO
    taskname1, taskname2, taskname3, taskname4, taskname5
    .endhelp

Also please be aware of the following limitations:

#. The tasks referenced in ``SEE ALSO`` need to have documentation themselves.
#. All tasks in ``SEE ALSO`` must be on the same line in order for the parser to correctly generate links.
#. Literal blocks need to be enclosed by quotes::

    This is an example "/path/to/some/file.txt", and "parameter='something'" option.

**************
Error handling
**************

PySALT error handling is based on Pythons build in `exception handling <http://docs.python.org/tutorial/errors.html>`_.
This basically works in the following way.
If an error occurs an *Exception* is *raised* this will cause the program to fall through untill either the exception is caught by an *except* statement or the program terminates.
The *except* statement can catch an error and do something with it.
It may simply do some cleanup and quit or it can undertake a specific action (such as replacing division by zero error results by Inf values).

Implementing error handling in PySALT
=====================================

In order for errors to be correctly handled by PySALT every task needs to have the following structure::

    from saltsafelog import logging
    from salterror import SaltError, SaltIOError

    def fragile_function(x):
        if x<0:
           raise SaltIOError('cannot open file with number smaller than zero')
        return open('test-'+str(x)+'.fits')

    def sampletask(parameter1,parameter2,clobber,logfile,verbose,debug):
        # Enable logging
        with logging(logfile,debug) as log:
            # Some work that produces a_result

            if a_result==None:
                raise SaltError('cannot handle null results')

            try:
                b_result=fragile_function(a_result)
            except SaltIOError:
                b_result=fragile_function(abs(a_result))

            # Some more work that needs b_result
            log.message('everything is fine here')

    # -----------------------------------------------------------
    # main code 
    parfile = iraf.osfn("package$sampletask.par") 
    t = iraf.IrafTaskFactory(taskname="sampletask",value=parfile,function=sampletask,pkgname='package')

Let's go over this block by block.
The error handling mechaninsm for PySALT is implemented in two modules::

    from saltsafelog import logging
    from salterror import SaltError, SaltIOError

* *salterror* contains a tree structure of exceptions derived from *Exception*. The class *SaltError* is the top level of this tree and an error of this type should be raised if it does not fit into any of it's child classes.
* *saltsavelog* provides the mechanism for automatic logging of errors, warnings and messages in a simple way using the *with statement*

When your code finds an error it should *raise* an exception which is (if not caught in the current function) automatically passed to the containing function::

    if x<0:
        raise SaltIOError('cannot open file with number smaller than zero')


The following line sets up the context manager which simplifies logging::

    with logging(logfile,debug) as log:

It creates a *log* object which has two methods:

* *warning()* to log warnings
* *message()* to log messages

A header with the date and time the method was called is automatically added and the text is displayed in the terminal and written to the specified logfile.
Errors derived from *SaltError* which reach the end of the *with* block are also caught and logged after which the program automatically exits.
When the *debug* parameter is *True* a backtrace containing usefull information for developers is written to the logfile as well.
Furthermore the context manager automatically logs the call to the function in which it is contained allong with the parameters used to call it.


The following code shows you how to handle non-fatal errors (e.g. errors for which you do not want the routine to exit)::

    try:
        b_result=fragile_function(a_result)
    except SaltIOError:
        b_result=fragile_function(abs(a_result))

The program tries to execute the statements in the *try* block and if an error of the type *SaltIOError* occurs the code in the *except* block get's executed.
See the `python exception handling pages <http://docs.python.org/tutorial/errors.html>`_ for more information.

The last part is not related to the error handling but tells pyraf which function to execute when the task is run::

    # -----------------------------------------------------------
    # main code 
    parfile = iraf.osfn("package$sampletask.par") 
    t = iraf.IrafTaskFactory(taskname="sampletask",value=parfile,function=sampletask,pkgname='package')

Background information
======================

This section explains a bit more about the idea behind the error handling mechanism in PySALT.
The basic premisis is that there are three kinds of errors.

#. Errors that are expected to occur and should not cause the program to stop. For instance when one of the stars is not found in a particular frame you just want the frame to be skipped.
#. Errors that are expected to occur but the program cannot continue correctly. For instance when a user gives a non existing filename as input or the file cannot be opened because of wrong permissions. This should log an error the user can comprehend and quit.
#. Errors that are unexpected. By definition these errors constitute programming bugs and should cause the program to stop with some type of error message usefull to the developer.

In order to handle all three cases in a consistent manner the *salterror* introduces a tree of exception classes.
As in the whole of Python the top level error class for PySALT e.g. *SaltError* is derived from Pythons build in *Exception* class.
All errors of the first two categories should raise an exception derived from *SaltError* this means that giving *SaltError* as an exception type to accept will cause all expected errors to be caught leaving only the unexpected errors which will still cause the program to crash as required.
Catching also all unexpected errors is considered bad programming because the containing program can than no longer insure if everything went ok (this includes the shell which was used to launch pyraf in pipeline mode).

Try, except blocks can be used to catch errors that should not cause the program to crash and handle them accordingly.
All *SaltError* derived errors that remain are caught and handled by the *logging* context manager.

==About the context manager==
The line::

    with logging(logfile,debug) as log:

starts the logging context manager.
Uppon reaching this line the function *saltsafelog.logging()* is executed.

The code in this function first creates a *SaltLog* object.
This object has methods which write messages, errors and warnings to a logfile adding a header with the current date and time automatically.

Then *logging()* uses object introspection from the Python standard library *inspect* module to get the names and values of all the parameters of the function in which *logging()* is called (so the function containing the with statement).
This information is written to the logfile using the *message()* method from *SaltLog*.

Then a try block is started and the control is passed back to the calling function by yielding the *SaltLog* object as *log* (a function that yields a value instead of returns it is called a generator function and this function returns to the point right after the yield statement when called again).
The whole program just runs and may use the methods from the yielded *log* object to also write messages or warnings to the logfile (avoid writing errors directly using the *error()* method, rather raise errors).

When an exception derived from *SaltError* is raised which is not caught by some except block (so this should be case 2) the program will fall through and return to the try block in the *logging* context manager.
The error will be caught by it's except statement and the error message plus an additional backtrace (if debug=True) will be written to the logfile.
After that the program terminates normally.

***************
Working with Qt
***************

As of PySALT version 0.40 the standard toolkit used to build Graphical User Interface or GUI programs is `PyQt4 <http://www.riverbankcomputing.co.uk/software/pyqt/intro>`_.
This is a set of Python wrappers arround the commonly used `Qt <http://qt.nokia.com>`_ UI framework.
Qt, developed by Nokia, is the free and open souce cross platform GUI toolkit used to build KDE as well as many commercial applications.

The following tasks are currently based on PyQt4:

* slotpreview
* bvitpreview
* bvitview

In order to work with these packages you need a working installation of Qt version 4.0 or higher as well as PyQt4 itself.
All major Linux distributions include binaries for both systems and there is a simple DMG installer for Mac OSX.

Basic concepts for developers
=============================

PyQt is very easy to work with especially if one is used to other GUI programming toolkits.
A graphical user interface consists of a set of objects called widgets which are combined to form new widgets (for instance a text label widget ("do you want to quit?") combined with two push button widgets (ok/cancel) forms a dialog widget which may be used in combination with other widgets to form a program.
The widgets themselves are all simply instances of classes which define how they look and behave.
Qt comes with a lot of build in widgets for commonly used elements such as: radio buttons, check boxes, line edit boxes, etc.
Custom widgets can easily be created by subclassing an existing widget class.

For a quick and easy introduction to PyQt4 see `the PyQt4 tutorial <http://zetcode.com/tutorials/pyqt4>`_.

Creating a GUI with Qt designer
===============================

All interfaces can in principle be hand coded.
For small programs this is usually fine, but for larger projects this can be quite time consuming.
Qt comes with a program called *Qt designer* to simplify this task.
With designer you can create your interface widget using simple drag and drop.
You may add and position elements and set default values.
It is even possible to define basic interactions between elements in your widget such as: if this checkbox is clicked then hide that text input box.
This is great because it allows you to separate the design of the user interface from the program code thereby making it much easier to change something without having to do a full rewrite.
Also it allows users to test the interface before any code has been written so you can experiment with different designs.

Building the GUI tools for PySALT
=================================

In principle you should not need to rebuild the tools since the resulting Python modules are already included in the repository.
However if you are a PySALT developer and want to make changes the following procedure should be followed.

* The user interfaces for PySALT have all been created using Qt designer. These interfaces are saved as *.ui* files in the *lib/ui* directory. In order to modify the interface simply open the corresponding *.ui* file in Qt designer and make the desired modifications.
* Now that you have updated the *.ui* file you need to convert it into the corresponding Python module. This is done by invoking the *pyuic4* tool. For instance to rebuild the photometry star selection dialog run::

    pyuic4 -o ui_photometryconfigwidget.py PhotometryConfigWidget.ui

* If you want to include additional icons in your GUI you should save them as *.png* images in the *lib/ui/icons* directory and add them to the resource file *lib/ui/icons.qrc*. This file also needs to be compiled using pyrcc4. pyrcc4 is PyQt's equivalent to Qt's rcc  utility and is used in exactly the same way. pyrcc4  reads the .qrc  file, and the resource files, and generates a Python module that only needs to be import ed by the application in order for those resources to be made available just as if they were the original files.
* Now you may use the interface from your Python module by importing the generated widget class. For instance::

    # Gui library imports
    from PyQt4 import QtGui, QtCore

    class PhotometryConfigWidget(QtGui.QWidget):
        """Configure dialog for photometry
        """
        def __init__(self, parent=None):
            """Setup widget.
            *parent* parent widget.
            """

            # Import gui
            from ui.ui_photometryconfigwidget import Ui_PhotometryConfigWidget

            # Setup widget
            QtGui.QWidget.__init__(self, parent)

            # Bind gui to widget
            self.ui = Ui_PhotometryConfigWidget()
            self.ui.setupUi(self)

****
Todo
****

.. todolist::