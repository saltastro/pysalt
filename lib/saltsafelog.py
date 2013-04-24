################################# LICENSE ##################################
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.                                                     #
#                                                                          #
# Redistribution and use in source and binary forms, with or without       #
# modification, are permitted provided that the following conditions       #
# are met:                                                                 #
#                                                                          #
#     * Redistributions of source code must retain the above copyright     #
#       notice, this list of conditions and the following disclaimer.      #
#     * Redistributions in binary form must reproduce the above copyright  #
#       notice, this list of conditions and the following disclaimer       #
#       in the documentation and/or other materials provided with the      #
#       distribution.                                                      #
#     * Neither the name of the South African Astronomical Observatory     #
#       (SAAO) nor the names of its contributors may be used to endorse    #
#       or promote products derived from this software without specific    #
#       prior written permission.                                          #
#                                                                          #
# THIS SOFTWARE IS PROVIDED BY THE SAAO ''AS IS'' AND ANY EXPRESS OR       #
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED           #
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE   #
# DISCLAIMED. IN NO EVENT SHALL THE SAAO BE LIABLE FOR ANY                 #
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL       #
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS  #
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)    #
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      #
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN #
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE          #
# POSSIBILITY OF SUCH DAMAGE.                                              #
############################################################################

"""Module provides uniform logging for SALT."""

# Ensure Python 2.5 compatibility
from __future__ import with_statement

import inspect
from salterror import SaltError
from contextlib import contextmanager
from time import strftime

class SaltLog:
    """Class providing uniform logging."""

    def __init__(self,logfile,with_stdout=True):
        """Constructor, binds logfile to file object *f*.
        
	"""

        # Bind to logfile
        self.logfile=logfile


    def error(self,e,with_backtrace=True, with_stdout=True):
       """Prints error message and backtrace to logfile."""

       with open(self.logfile,'a') as f:
            import traceback

            # Get current time
            time=strftime("%Y-%m-%d %H:%M:%S")

            # Define header
            header="%s ERROR --------------------------------------------\n" % time

            # Optional write to standard output
            if with_stdout:
                # Print header
                print header

                # Print message
                print str(e)+'\n'

                if with_backtrace:
                    # Print error traceback
                    traceback.print_exc()

            # Write header to logfile
            f.write(header)

            # Write message to logfile
            f.write(str(e)+'\n')

            if with_backtrace:
                # Write error traceback to logfile
                traceback.print_exc(file=f)


    def warning(self,m, with_stdout=True):
        """Prints warning message *m* to logfile."""

        with open(self.logfile,'a') as f:
            # Get current time
            time=strftime("%Y-%m-%d %H:%M:%S")

            # Define header
            header="%s WARNING ------------------------------------------\n" % time

            # Compose final message
            log_message=header+m+'\n'

            # Optional print to standard output
            if with_stdout:
                print log_message

            # Write header+message to logfile
            f.write(log_message)

    def message(self,m,with_header=True, with_stdout=True):
        """Prints message *m* to logfile."""

        with open(self.logfile,'a') as f:
            # Get current time
            time=strftime("%Y-%m-%d %H:%M:%S")

            # Define header
            header="\n%s MESSAGE ------------------------------------------\n" % time

            # Compose final message
            if with_header:
                log_message=header+m
            else:
                log_message=m

            # Optional print to standard output
            if with_stdout:
                print log_message

            # Write header+message to logfile
            f.write(log_message+'\n')

def history(level=3, wrap=True, wrapchar=80, exclude=[]):
   """Return the history of the call.  This includes return the name of the
       current program as well as the information about all the parameters
       which are passed to it.  The return is the name of the program along 
       witha string listing all the parameters
 
       level -- the level of the frame of interest
       wrap  -- wraps characters if true
       wrapchar--number of characters to wrap at
       exclude--options to exclude 
 
       returns str, str 
   """
   frame=inspect.getouterframes(inspect.currentframe())[level][0]
   args,_,_,values=inspect.getargvalues(frame)
   fname=str(inspect.getframeinfo(frame)[2])

   #msg='%s started with the following parameters:\n\n' % fname
   msg ='%s ' % fname.upper()
   lcount=0
   for i in args:
     if  i not in exclude: 
       instr="%s=%s," % (i, values[i])
       lcount += len(instr)
       if lcount>wrapchar and wrap: 
           msg += '\n'
           lcount = len(instr)
       if i.count('pass'):
           msg+="%s=%s " % (i, '****')
       else:
           msg+="%s=%s " % (i, values[i])
   return fname, msg

@contextmanager
def logging(logfile,debug=True,with_stdout=True,with_call=True):
    """Context manager to ensure proper error handling and logging.

    Example usage::

        with logging('logfile.txt') as log:
            # User code
            log.message('Hello world!') # Writes a message to the log
            # Some more user code
            log.warning('This is a warning message') # Writes warning
            # Again some user code
            raise SaltError('Oops!') # Error message and traceback are written to log

    """

    # Create LogFile object
    log=SaltLog(logfile,with_stdout)

    # Log the call
    if with_call:
        #get the history of the program
        fname, msg=history(level=3, wrap=True, wrapchar=80)       

        log.message(msg.rstrip()+'\n', with_header=False)
        log.message("Starting %s\n" % fname)

    # Wrap code block in try, except to ensure proper error handling
    try:
        # Transfer control to wrapped code with access to the log
        yield log

        # Log completion
        if with_call:
            msg='%s completed' % fname
            log.message(msg)

    except SaltError, e:
        # Catch and log any errors that may have occured
        log.error(e,debug)

        # Log abort
        if with_call:
            msg='%s aborted\n' % fname
            log.message(msg)

        #raise the error to quit out of the program and allow any
        #wrapper program to catch the error
        raise SaltError(e)

    finally:
        # Any additional cleanup code and logging goes here
        pass

