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

""" Inserts the directory containing links to all Python-based tasks
    in the SALT tree to the default Python path.
"""
import os, sys
from pyraf import iraf 
import matplotlib
matplotlib.use('Agg')


# define and add path to top level Python directory in SALT

_path = iraf.osfn('pysalt$')
if _path not in sys.path:
    sys.path.insert(1,_path)

# define and add path to lib directory in SALT package

_path = iraf.osfn('pysalt$lib')
if _path not in sys.path:
    sys.path.insert(1,_path)

# define and add path to directories in SALT package

_path = iraf.osfn('pysalt$saltred')
if _path not in sys.path:
    sys.path.insert(1,_path)

_path = iraf.osfn('pysalt$slottools')
if _path not in sys.path:
    sys.path.insert(1,_path)

_path = iraf.osfn('pysalt$saltspec')
if _path not in sys.path:
    sys.path.insert(1,_path)

_path = iraf.osfn('pysalt$saltfp')
if _path not in sys.path:
    sys.path.insert(1,_path)

_path = iraf.osfn('pysalt$bvittools')
if _path not in sys.path:
    sys.path.insert(1,_path)

#add the proptools class
_path = iraf.osfn('pysalt$proptools')
if _path not in sys.path:
    sys.path.insert(1,_path)

#Add the plugin directory if it exists
_path = iraf.osfn('pysalt$plugins')
if _path not in sys.path and os.path.isdir(_path):
    sys.path.insert(1,_path)

#Add the plugin directory if it exists
_path = iraf.osfn('pysalt$saltfirst')
if _path not in sys.path and os.path.isdir(_path):
   sys.path.insert(1,_path)

#add the proptools class
_path = iraf.osfn('pysalt$salthrs')
if _path not in sys.path:
    sys.path.insert(1,_path)
