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

"""Module for printing.

.. warning::
    This module is **deprecated** and has been superceded by :mod:`saltsafeio`.
    Please adjust your code accordingly.

"""

def log(file,message,verbose):

    if (verbose and verbose != 'no'):

# print to shell

        print (message)

# print to log file

        try:
            output = open(file,'a')
            output.write(message+'\n')
            output.close()
        except IOError:
            print "Could not write message '"+message+"' to file "+file

def err(file,message):
    """write message to log file and shell"""

    log(file,message,True)
    return 1

def time(text,file,verbose):
    """write time to log file and shell"""

    import time

    if (verbose):
        message = text + ': '+time.asctime(time.localtime())
        log(file,message,verbose)

def line(file,verbose):
    """write line descriminator to log file and shell"""

    if (verbose):
        message = '----------------------------------------------------------------------------'
        log(file,message,True)

def history(taskname, plist, slist, logfile, verbose):
    """Using the parameter list and the string list, create the output history

    return  history
    """
    status=0
    history=''
    if len(plist)!=len(slist):
        message='SALTPRINT.history--ERROR:  Input parameters not equal to string parameters'
        status=err(logfile, message)
        return history, status

    try:
        history='%s--' % taskname.upper()
        for i in range(len(plist)):
            history += '%s=%s ' % (slist[i], plist[i])
    except Exception, e:
        message='SALTPRINT.history--ERROR:  Could not create output history string because %s' % e
        status=err(logfile, message)

    return history, status
