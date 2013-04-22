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

procedure pyexecute(filename)

string filename {prompt="Filename of Python module to be executed."}
string tasknames = "" {prompt="Name(s) of IRAF task defined by Python module."}
bool verbose=yes {prompt="Warn about existence of Python tasks."}

string mode="h"

string *list

begin

    string curpack, tmpfile, i_file
    string tasknm, i_tasks, doit
    int i, i1, i2, istr

    i_tasks = tasknames
    i_file = filename
        
	# get current package
	tmpfile = "tmp$" // mktemp("pyexecute")
	package(> tmpfile)
	list = tmpfile
	curpack = list
	list = ""
	delete(tmpfile, verify=no, go_ahead=yes)
	for (i=1; i<=strlen(curpack); i = i+1) {
		if (substr(curpack,i,i) != " ") {
			curpack = substr(curpack,i,strlen(curpack))
			break
		}
	}
    # Warn user that this package contains a PyRAF task
	if (verbose) {
		print ("Warning: package `", curpack,
			"' includes Python tasks that require PyRAF")
	}

    # Set up the hidden tasks for each named PyRAF task
    # that will redirect user calls to the task to a dummy
    # script which warns the user that it is a PyRAF task
    #
    # Count number of tasks listed
    istr = strlen(i_tasks)
    
    if (istr > 0) {
   
        # For each task listed...
        i1 = 1
        i2 = 0
        i = 1
        while (i <= istr){
            if (substr(i_tasks,i,i) == ",") {
               tasknm = substr(i_tasks,i1,i-1)
               # Try to ignore blanks
               if (tasknm != "") {
                    doit = "task $"+tasknm+" = \"pysalt$nopyraf.cl\" ; keep "
                    print (doit) | cl
                    keep
                    i2 = 0
               }
               i1 = i+1
            }
            i = i+1
            i2 = 1
        }
        if (i2 == 1)
            tasknm = substr(i_tasks,i1,istr)
            if (tasknm != "") {
                doit = "task $"+tasknm+" = \"pysalt$nopyraf.cl\" ; keep "
                print (doit) | cl
                keep
            }
    }
end
