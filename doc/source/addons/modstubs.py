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

# Ensure python 2.5 compatibility
from __future__ import with_statement

import os

"""Autogenerate module documentation in *modules* directory for all modules listed in *modules.lst*"""

def autogen(list,output_path):
    """
    Automatically generates module documentation files for modules in list file.
    *list* filename of module list.
    *output_path* output module documentation goes in this directory.
    """

    print 'Autogenerating module documentation for modules in',list

    # Read modules from list file
    with open(list) as f:
        modules=[m.strip() for m in f.readlines() if m.strip()[0]!='#']

    for m in modules:
        # Try importing the module this is needed for sphinx to function correctly
        try:
            exec('import %s' % m)
            exec('del %s' % m)
        except:
            print 'Module',m,'cannot be imported no documentation will be generated.'
            continue

        # Check if documentation file exists
        if os.path.isfile(output_path.rstrip('/')+'/'+m+'.rst'):
            print 'Module documentation for',m,'exists, skipping.'
            continue

        # Write empty configuration file
        with open(output_path.rstrip('/')+'/'+m+'.rst','w') as f:
            f.write('*'*len(m)+'\n')
            f.write(m+'\n')
            f.write('*'*len(m)+'\n')
            f.write('\n')
            f.write('.. automodule:: '+m+'\n')
            f.write('   :members:\n')
            f.write('   :undoc-members:\n')
            f.write('   :show-inheritance:\n')
            f.write('\n')
