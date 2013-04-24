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

"""
Module containing generic graphical user interface widgets.
"""

# Ensure python 2.5 compatibility
from __future__ import with_statement
import matplotlib.cm

# General imports
import pyfits
import numpy as np

# Gui library imports
from PyQt4 import QtGui, QtCore
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import CirclePolygon, Rectangle

# Salt imports
import saltsafeio
from salterror import SaltError, SaltIOError
from saltimagetools import find_object, zscale

class MplCanvas(FigureCanvas):
    """Base class for embedding a matplotlib canvas in a PyQt4 GUI.
    """

    def __init__(self):
        """Default constructor."""

        # Initialize base class
        FigureCanvas.__init__(self,Figure())

        # Set resize policy
        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)

        # Set geometry
        FigureCanvas.updateGeometry(self)

    def set_focus(self, event):
        self.setFocus()

    def connectMatplotlibKeyEvents(self):
        """Bind events to event handlers."""

        self.key_press_id = self.figure.canvas.mpl_connect('key_press_event', self.onKeyPress)
        self.key_release_id = self.figure.canvas.mpl_connect('key_release_event', self.onKeyRelease)
    
    def connectMatplotlibMouseEvents(self):
        """Bind events to event handlers."""

        self.button_press_id = self.figure.canvas.mpl_connect('button_press_event', self.onButtonPress)
        self.button_release_id = self.figure.canvas.mpl_connect('button_release_event', self.onButtonRelease)

    def connectMatplotlibMouseMotion(self):
        self.mouse_motion_id = self.figure.canvas.mpl_connect('motion_notify_event', self.set_focus)

    def disconnectMatplotlibKeyEvents(self):
        """Unbind events."""

        self.figure.canvas.mpl_disconnect(key_press_id)
        self.figure.canvas.mpl_disconnect(key_release_id)

    def disconnectMatplotlibMouseEvents(self):
        """Unbind events."""

        self.figure.canvas.mpl_disconnect(button_press_id)
        self.figure.canvas.mpl_disconnect(button_release_id)

    def onButtonPress(self, event):
        """Overload this function to implement mousebutton press events."""
        pass

    def onButtonRelease(self, event):
        """Overload this function to implement mousebutton release events."""
        pass

    def onKeyPress(self, event):
        """Overload this function to implement mousebutton press events."""
        print "I'm here:", event, dir(event)
        pass

    def onKeyRelease(self, event):
        """Overload this function to implement mousebutton release events."""
        pass

class ImageDisplay(MplCanvas):
    """Class for displaying FITS images using matplotlib imshow() embedded in a Qt 4 GUI.
    With extra methods for overplotting patches."""

    def __init__(self):
        """Default constructor."""

        # Initialize base class
        MplCanvas.__init__(self)

        # Add central axes instance
        self.axes = self.figure.add_subplot(111)

        # Connect mouse events
        self.connectMatplotlibMouseEvents()

        # Keep track of all patches
        self.patches={}
        self.zorder=10

        # Set display parameters
        self.vmin=None
        self.vmax=None
        self.interpolation=None
        self.cmap=None
        self.scale='minmax'
        self.contrast=1.0
        self.origin='lower'
        self.aspect='auto'

    def onButtonPress(self, event):
        """Emit signal on selecting valid image position."""

        if event.xdata and event.ydata:
            self.emit(QtCore.SIGNAL("positionSelected(float, float)"), 
	            float(event.xdata), float(event.ydata))

    def setColormap(self, cmap_name):
        """Set colormap based on name."""
        try:
            self.cmap=matplotlib.cm.get_cmap(cmap_name)
        except:
            raise SaltError('Cannot get colormap instance for specified cmap')

    def setScale(self):
        # Set scale parameters
        if self.scale=='minmax':
            self.vmin=np.min(self.image)
            self.vmax=np.max(self.image)
        elif self.scale=='zscale':
            self.vmin,self.vmax=zscale(self.image,self.contrast)
        else:
            self.vmin=None
            self.vmax=None

    def loadImage(self, image):
        """Load image array."""

        # Set image
        self.image=image

        # Set vmin and vmax parameters
        self.setScale()

    def drawImage(self):
        """Draw image to canvas."""

        # Display image
        self.axes.imshow(self.image, cmap=self.cmap, aspect=self.aspect, vmin=self.vmin, vmax=self.vmax,interpolation=self.interpolation, origin=self.origin)

    def addPatch(self, label, patch):
        """Add a matplotlib *patch* instance with a given *label*."""

        # There shall be one and only one patch for each label
        if label in self.patches:
            del self.patches[label]

        # Add patch to list
        self.patches[label]=patch

        self.zorder+=1
        
    def removePatch(self, label):
        """Remove patch instance referenced by *label* from figure."""

        # Remove patch if it exists
        if label in self.patches:
            del self.patches[label]

    def addCircle(self, label, x, y, r, color='y', lw=1):
        """Add circle patch at postion (*x*,*y*) with radius *r*
        using a line with color *color* and thickness *lw*."""

        circ=CirclePolygon((x,y),radius=r,ec=color,zorder=self.zorder,lw=lw,fill=False)

        # Add patch to figure
        self.addPatch(label,circ)

    def addSquare(self, label, x, y, r, color='y', lw=1):
        """Add square patch at postion (*x*,*y*) with radius *r*
        using a line with color *color* and thickness *lw*."""

        # Calculate coordinates
        xl=x-r
        yl=y-r
        w=2*r
        h=2*r

        # Create Rectangle patch instance
        rect=Rectangle((xl,yl),width=w,height=h,ec=color,zorder=self.zorder,lw=lw,fill=False)

        # Add patch to figure
        self.addPatch(label,rect)

    def addRectangle(self, label, x1, y1, x2, y2, color='y', lw=1):
        """Add rectangle patch from (*x1*,*y1*) to (*x2*,*y2*)
        using a line with color *color* and thickness *lw*."""

        # Calculate coordinates
        w=x2-x1
        h=y2-y1

        # Create Rectangle patch instance
        rect=Rectangle((x1,y1),width=w,height=h,ec=color,zorder=self.zorder,lw=lw,fill=False)

        # Add patch to figure
        self.addPatch(label,rect)

    def reset(self):
        # Delete all patches
        self.patches={}

        # Redraw canvas
        self.redraw_canvas()

    def redraw_canvas(self,keepzoom=False):
        if keepzoom:
            # Store current zoom level
            xmin, xmax = self.axes.get_xlim()
            ymin, ymax = self.axes.get_ylim()

        # Clear plot
        self.axes.clear()

        # Draw image
        self.drawImage()

        # Draw patches
        for key in self.patches.keys():
            self.axes.add_patch(self.patches[key])

        if keepzoom:
            # Restore zoom level
            self.axes.set_xlim((xmin,xmax))
            self.axes.set_ylim((ymin,ymax))

        # Force redraw
        self.draw()

