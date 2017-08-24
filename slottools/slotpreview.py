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

"""Preview slotmode data and set-up config file for photometry."""

# Ensure Python 2.5 compatibility
from __future__ import with_statement

# General imports
import sys
import os
import pyfits
from pyraf import iraf
from pyraf.iraf import pysalt

# Gui library imports
from PyQt4 import QtGui, QtCore
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT

# Salt imports
import saltsafeio
from saltsafelog import logging
from saltgui import ImageDisplay
from PhotometryConfigWidget import PhotometryConfigWidget
from salterror import SaltIOError

debug=True

class ApplicationWindow(QtGui.QMainWindow):
    """Main application window."""

    def __init__(self, imlist, number, config,
                target_line_color='g', comparison_line_color='g',
                target_line_width=2, comparison_line_width=2,
                distance=5, cmap='gray', scale='zscale', contrast=1.0):
        """Default constructor."""

        # Setup widget
        QtGui.QMainWindow.__init__(self)

        # Set main widget
        self.main = QtGui.QWidget(self)

        # Set window title
        self.setWindowTitle("Slotpreview")

        # Add FITS display widget with mouse interaction and overplotting
        self.imdisplay = ImageDisplay()
        self.imdisplay.setMinimumHeight(150)
        self.imdisplay.setMinimumWidth(400)

        # Set colormap
        self.imdisplay.setColormap(cmap)

        # Set scale mode for dynamic range
        self.imdisplay.scale=scale
        self.imdisplay.contrast=contrast

        # Add navigation toolbars for each widget to enable zooming
        self.toolbar=NavigationToolbar2QT(self.imdisplay,self)

        # Add configuration widget and connect it to the display widget
        self.conf = PhotometryConfigWidget(imdisplay=self.imdisplay, imlist=imlist, number=number, config=config)

        # Set default linewidth and color
        self.conf.setSearchDistance(distance)
        self.conf.setMarkerType('square')
        self.conf.setLineColor('target', target_line_color)
        self.conf.setLineColor('comparison', comparison_line_color)
        self.conf.setLineWidth('target', target_line_width)
        self.conf.setLineWidth('comparison', comparison_line_width)

        # Layout the widgets
        l = QtGui.QVBoxLayout(self.main)
        l.addWidget(self.imdisplay)
        l.addWidget(self.toolbar)
        l.addWidget(self.conf)

        # Set focus to main widget
        self.main.setFocus()

        # Set the main widget as the central widget
        self.setCentralWidget(self.main)

        # Destroy widget on close
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)

        # Close when config dialog is closed
        self.connect(self.conf, QtCore.SIGNAL('destroyed()'), self, QtCore.SLOT('close()'))

def slotpreview(images,outfile,ampperccd=2,ignorexp=6,recenter_radius=5,
                tgt_col='b',cmp_col='g', tgt_lw=2,cmp_lw=2,cmap='gray',
                scale='zscale',contrast=0.1,clobber=True,logfile='salt.log',
                verbose=True):

    with logging(logfile,debug) as log:
        
        # is the input file specified?
        saltsafeio.filedefined('Input',images)

        # if the input file is a list, does it exist?
        if images[0] == '@':
            saltsafeio.listexists('Input',images)

        # parse list of input files
        imlist=saltsafeio.listparse('Raw image',images,'','','')

        # check input files exist
        saltsafeio.filesexist(imlist,'','r')

        # is the output file specified?
        saltsafeio.filedefined('Output',outfile)

        # check output file does not exist, optionally remove it if it does exist
        if os.path.exists(outfile) and clobber:
            os.remove(outfile)
        elif os.path.exists(outfile) and not clobber:
            raise SaltIOError('File '+outfile+' already exists, use clobber=y')

        # Get the number of ccds to calculate the number of frames to skip
        try:
            nccds=pyfits.getheader(imlist[0])['NCCDS']
        except:
            raise SaltIOError('Could not read NCCDS parameter from header of first fits file.')

        # Create GUI
        App = QtGui.QApplication(sys.argv)
        aw = ApplicationWindow(imlist=imlist,
                number=ignorexp*ampperccd*nccds+1,
                config=outfile, target_line_color=tgt_col,
                comparison_line_color=cmp_col, target_line_width=tgt_lw,
                comparison_line_width=cmp_lw, distance=recenter_radius,
                cmap=cmap, scale=scale, contrast=contrast)
        aw.show()

        # Start application event loop
        exit=App.exec_()

        # Check if GUI was executed succesfully
        if exit!=0:
            log.warning('Slotpreview GUI has unexpected exit status'+str(exit))

# -----------------------------------------------------------
# main code 

parfile=iraf.osfn("slottools$slotpreview.par") 
t=iraf.IrafTaskFactory(taskname="slotpreview",value=parfile,function=slotpreview,pkgname='slottools')

