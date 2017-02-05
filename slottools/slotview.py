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

"""SlotView is the interactive tool for analyzing lightcurves produced by 
   slotphot for SALT SALTICAM slotmode data.

Updates:

20110516
    * Updated the code to use PyQt4 for the GUI backend
    * Updated the code to handle the new error handling


"""
# Ensure python 2.5 compatibility
from __future__ import with_statement


import time, math
import numpy as np
import scipy as sp
from pyraf import iraf
from pyraf.iraf import pysalt
import saltprint, salttime
import slottool as st

import Tkinter as Tk
from matplotlib.widgets import Cursor, SpanSelector, Slider, CheckButtons

# Gui library imports
from PyQt4 import QtGui, QtCore
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT

# Salt imports
from saltgui import ImageDisplay, MplCanvas
from salterror import SaltIOError

import saltsafeio as saltio
import saltsafekey as saltkey
from saltsafelog import logging

from SlotViewWindow import SlotViewWindow

debug=True

def slotview(newfits,indata , fileout, srcfile, fps=10.0, phottype='square', sigdet=5, contpix=10, \
    driftlimit=10, clobber=True,logfile='slotview.log',verbose=True):

    
#set up the variables
    status = 0
    entries = []
    vig_lo = {}
    vig_hi = {}
    hour = 0
    min = 0
    sec = 0.
    time0 = 0.
    nframes = 0
    sleep=0

    with logging(logfile,debug) as log:

        #enter in the input data
        saltio.fileexists(newfits)

        #set the sleep parameter
        if fps>0: sleep=1.0/(fps)

        # read in the data file
        id, time, ratio, rerr, tx, ty, tflux, terr, cx, cy, cflux, cerr=st.readlcfile(indata)

        # read extraction region defintion file
        amp, x, y, x_o, y_o, r, br1, br2=st.readsrcfile(srcfile)



        #determine the size of the data arrays
        struct = saltio.openfits(newfits)
        naxis1 = saltkey.get('NAXIS1',struct[1])
        naxis2 = saltkey.get('NAXIS2',struct[1])

        # Plot all of the data and the first image
        # Create GUI
        App = QtGui.QApplication([])
        aw=SlotViewWindow(struct, id, tflux, cflux, ratio, time, phottype, sleep,   \
                     tx, ty, cx, cy, r, br1, br2, naxis1, naxis2, sigdet, contpix, driftlimit)
        aw.show()

        # Start application event loop
        app_exit=App.exec_()
        
        # Check if GUI was executed succesfully
        if app_exit!=0:
             raise SALTError('InterIdentify GUI has unexpected exit status '+str(exit))

        ratio, tflux, cflux, gframe, newphot=aw.ratio, aw.tflux, aw.cflux, aw.goodframes, aw.newphot

        #close the input file
        saltio.closefits(struct)

        # Update the indata file if necessary
        lc=saltio.openascii(fileout,'w')

        for i in range(len(ratio)):
            x['target']=tx[i]
            x['comparison']=cx[i]
            y['target']=ty[i]
            y['comparison']=cy[i]
            reltime=False
            if gframe[i]:
                st.writedataout(lc, id[i], time[i], x, y, tflux[i], terr[i], \
                    cflux[i], cerr[i], ratio[i], rerr[i], time[0], reltime)

        saltio.closeascii(lc)

# -----------------------------------------------------------
# main code

parfile = iraf.osfn("slottools$slotview.par")
t = iraf.IrafTaskFactory(taskname="slotview",value=parfile,function=slotview, pkgname='slottools')
