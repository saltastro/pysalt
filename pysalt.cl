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

# Package script for the PySALT package

# IRAF patch level V2.12.2a or greater is required.
if (defpar ("release")) {
    if (release <= "2.12.1") {
        printf ("WARNING: IRAF patch level V2.12.2 or greater is required\n")
        printf ("         to run the SALT IRAF package\n")
        sleep 1
    }
} else {
    printf ("WARNING: IRAF patch level V2.12.2 or greater is required\n")
    printf ("         to run the SALT IRAF package\n")
    sleep 1
}
;

# Load necessary packages - only those that are used by all packages

# Set the imtype to fits
reset imtype = "fits"
flpr

# Set up tasks which report on PyRAF-based tasks
task pyexecute = "pysalt$pyexecute.cl"
hidetask pyexecute 

# This task simply prints a message stating that a task needs PyRAF to run
task nopyraf = "pysalt$nopyraf.cl"
hidetask nopyraf

# Add python tree to default Python path, if running PyRAF
pyexecute("pysalt$addpath.py",verbose=no)

#define the different tasks

set saltred = "pysalt$saltred/"
set slottools = "pysalt$slottools/"
set bvittools = "pysalt$bvittools/"
set saltfp = "pysalt$saltfp/"
set saltspec = "pysalt$saltspec/"
set proptools = "pysalt$proptools/"
set salthrs = "pysalt$salthrs/"

package pysalt

task saltred.pkg = saltred$saltred.cl

task slottools.pkg = slottools$slottools.cl

if (access('saltfp$saltfp.cl')) {
  task saltfp.pkg = saltfp$saltfp.cl
}

if (access('saltspec$saltspec.cl')) {
   task saltspec.pkg = saltspec$saltspec.cl
}

if (access('proptools$proptools.cl')) {
   task proptools.pkg = proptools$proptools.cl
}


if (access('salthrs$salthrs.cl')) {
   task salthrs.pkg = salthrs$salthrs.cl
}
if (access('bvittools$bvittools.cl')) {
   task bvittools.pkg = bvittools$bvittools.cl
}


task $sed = $foreign
hidetask sed

print(" ")
print("     +----------------------------------------------------+")
print("     |                  _ _                               |")
print("     |   __ _ _  __ __ | | |_    http://www.salt.ac.za    |")
print("     |  | _\ | || _| _`| |  _|   Southern African Large   |")
print("     |  |  /__ |__|\_,_|_|\__|   Telescope PyRAF Package  |")
print("     |  |_| \__'                                          |")
print("     |               Development PRERELEASE               |")
print("     |               Version 0.50  1 Oct 2014             |")
print("     |               Recommend IRAF 2.14/PyRAF 1.8.1      |")
print("     |               Bug reports: salthelp@salt.ac.za     |")
print("     +----------------------------------------------------+")
print(" ")
print("     Setting imtype=fits")
print(" ")


;
clbye()
