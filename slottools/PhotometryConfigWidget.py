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


class PhotometryConfigWidget(QtGui.QWidget):
    """Configure dialog for photometry.

    Has settings for:

    * target position, size
    * target background
        * type (anulus/region)
        * parameters
    * comparison position, size
    * comparison background
        * type (anulus/region)
        * parameters
    """

    def __init__(self, imdisplay, config, imlist=None, number=1, parent=None):
        """Setup widget.
        
        *imdisplay* a `FitsDisplay` derived fits display widget,
        *imlist* a list of fits image filenames,
        *config* filename used for output configuration file,
        *number* image number to load on startup,
        *parent* parent widget.
        """

        # Set default parameters
        self.imlist=imlist
        self.number=number
        self.config=config
        self.amp={'target' : 1, 'comparison' : 1 }

        # Set default marker
        self.mark_with='circle'

        # Set default search distance for recentering
        self.distance=5

        # Default line style parameters
        self.line={ 'target'     : { 'color' : 'g', 'width' : 2 },
                    'comparison' : { 'color' : 'g', 'width' : 2 }}

        # Import gui
        from ui_photometryconfigwidget import Ui_PhotometryConfigWidget

        # Setup widget
        QtGui.QWidget.__init__(self, parent)

        # Bind gui to widget
        self.ui = Ui_PhotometryConfigWidget()
        self.ui.setupUi(self)

        # Destroy widget on close
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)

        # Connect to display window
        self.imdisplay=imdisplay

        # Connect position selected signal from display to event handler
        self.connect(self.imdisplay, QtCore.SIGNAL('positionSelected(float, float)'), self.selectionHandler)

        # Set current display widget for positionSelected signal
        self.xdisplay=[]
        self.ydisplay=[]
        self.rdisplay=[]

        # Keep track of currently displayed objects
        self.display={'target'     : {'position' : False,
                                      'annulus'  : False,
                                      'region'   : False },
                      'comparison' : {'position' : False,
                                      'annulus'  : False,
                                      'region'   : False }}

        # Keep track of input widgets
        self.parameters=['x','y','r','r1','r2','x1','y1','x2','y2']

        self.input={'target'     : { 'x'  : self.ui.tgtXLineEdit,
                                     'y'  : self.ui.tgtYLineEdit,
                                     'r'  : self.ui.tgtRLineEdit,
                                     'r1' : self.ui.tgtR1LineEdit,
                                     'r2' : self.ui.tgtR2LineEdit,
                                     'x1' : self.ui.tgtX1LineEdit,
                                     'y1' : self.ui.tgtY1LineEdit,
                                     'x2' : self.ui.tgtX2LineEdit,
                                     'y2' : self.ui.tgtY2LineEdit},
                    'comparison' : { 'x'  : self.ui.cmpXLineEdit,
                                     'y'  : self.ui.cmpYLineEdit,
                                     'r'  : self.ui.cmpRLineEdit,
                                     'r1' : self.ui.cmpR1LineEdit,
                                     'r2' : self.ui.cmpR2LineEdit,
                                     'x1' : self.ui.cmpX1LineEdit,
                                     'y1' : self.ui.cmpY1LineEdit,
                                     'x2' : self.ui.cmpX2LineEdit,
                                     'y2' : self.ui.cmpY2LineEdit}}

        # Keep track of capture buttons
        self.buttons=['position','radius','annulus','region']

        self.capture={'target' \
                        : {'position' : self.ui.captureTgt,
                           'radius'   : self.ui.captureTgtRadius,
                           'annulus'  : self.ui.captureTgtAnulusBackground,
                           'region'   : self.ui.captureTgtRegionBackground},
                      'comparison' \
                        : {'position' : self.ui.captureCmp,
                           'radius'   : self.ui.captureCmpRadius,
                           'annulus'  : self.ui.captureCmpAnulusBackground,
                           'region'   : self.ui.captureCmpRegionBackground}}

        # Keep track of checkbox recenter widgets
        self.recenter={'target'     : self.ui.tgtRecenterCheckBox,
                       'comparison' : self.ui.cmpRecenterCheckBox}

        self.centered={'target'     : False,
                       'comparison' : False}

        # Enable blocking of redraws
        self.block={'target'     : { 'x'  : False,
                                     'y'  : False,
                                     'r'  : False,
                                     'r1' : False,
                                     'r2' : False,
                                     'x1' : False,
                                     'y1' : False,
                                     'x2' : False,
                                     'y2' : False},
                    'comparison' : { 'x'  : False,
                                     'y'  : False,
                                     'r'  : False,
                                     'r1' : False,
                                     'r2' : False,
                                     'x1' : False,
                                     'y1' : False,
                                     'x2' : False,
                                     'y2' : False}}

        # Set validator to ensure valid input on lineEdit input widgets
        self.validator = QtGui.QDoubleValidator(self)

        for object in ['target','comparison']:
            for key in self.parameters:
                self.input[object][key].setValidator(self.validator)

        # Set signal mapper for lineEdit updates
        self.drawMapper = QtCore.QSignalMapper(self)

        # Connect lineEdit updates to signal mapper 
        for object in ['target','comparison']:
            for key in self.parameters:
                # Add signal map entry
                self.drawMapper.setMapping(self.input[object][key],
                    QtCore.QString(object+','+key))

                # Connect to signal mapper
                self.connect(self.input[object][key], QtCore.SIGNAL('textChanged(QString)'), self.drawMapper, QtCore.SLOT('map()'))

        # Connect signal mapper to draw handler
        self.connect(self.drawMapper, QtCore.SIGNAL('mapped(QString)'),
            self.textUpdated)

        # Set signal mapper for capture buttons
        self.captureMapper = QtCore.QSignalMapper(self)

        # Connect capture button signals to signal mapper 
        for object in ['target','comparison']:
            for key in self.buttons:
                # Add signal map entry
                self.captureMapper.setMapping(self.capture[object][key],
                    QtCore.QString(object+','+key))

                # Connect to signal mapper
                self.connect(self.capture[object][key], QtCore.SIGNAL('clicked()'), self.captureMapper, QtCore.SLOT('map()'))

        # Connect signal mapper to capture handler
        self.connect(self.captureMapper, QtCore.SIGNAL('mapped(QString)'),
            self.captureHandler)

        # Connect save button
        self.connect(self.ui.saveButton, QtCore.SIGNAL('clicked()'), self.save)

        # If an image list is given
        if self.imlist is not None:
            # Connect image selection spinBox to event handlers
            self.connect(self.ui.imageSpinBox, QtCore.SIGNAL('valueChanged(int)'), self.loadImage)
            self.connect(self.ui.imageSpinBox, QtCore.SIGNAL('valueChanged(int)'), self.redraw)

            # Load first image
            self.setImageNumber(self.number)

            # Hide end selection widgets (not implemented here)
            self.ui.tgtEndPosLabel.hide()
            self.ui.tgtEndXLabel.hide()
            self.ui.tgtEndYLabel.hide()
            self.ui.cmpEndPosLabel.hide()
            self.ui.cmpEndXLabel.hide()
            self.ui.cmpEndYLabel.hide()
            self.ui.tgtXEndLineEdit.hide()
            self.ui.tgtYEndLineEdit.hide()
            self.ui.cmpXEndLineEdit.hide()
            self.ui.cmpYEndLineEdit.hide()
            self.ui.captureTgtEnd.hide()
            self.ui.captureCmpEnd.hide()

    def setImageNumber(self,number):
        """Set the image number."""

        self.ui.imageSpinBox.setValue(number)

    def loadImage(self, number):
        """Loads a new image.
        
        *number* is the image number to be loaded.

        This function uses `saltsafeio.getexposure` to get the correct
        exposure from a list of fits files containing an arbitrary number
        of extensions.
        """

        # Emit signal
        self.emit(QtCore.SIGNAL("imageNumberUpdated(int)"), number)

        # Load image from file
        self.img=saltsafeio.get_exposure(self.imlist,number)

        # Display image
        self.imdisplay.loadImage(self.img)

        # Redraw canvas
        self.imdisplay.redraw_canvas()

    def mark(self,*args,**kwargs):
        if self.mark_with=='square':
            self.imdisplay.addSquare(*args,**kwargs)
        elif self.mark_with=='circle':
            self.imdisplay.addCircle(*args,**kwargs)

    def textUpdated(self,key):
        # Get object and parameter from key
        obj,par=str(key).split(',')

        # Check block
        if self.block[obj][par]:
            return

        # Set block to prevent infinite repeat
        self.block[obj][par]=True

        # Recenter on object if requested
        if par=='x' and self.recenter[obj].isChecked() and not self.centered[obj]:
            x=float(self.input[obj]['x'].text())
            y=float(self.input[obj]['y'].text())
            r=float(self.input[obj]['r'].text())

            x,y=find_object(self.img,x,y,self.distance)
            
            self.input[obj]['x'].setText(str(x))
            self.input[obj]['y'].setText(str(y))

            self.centered[obj]=not(self.centered[obj])

        # Check if object region size locking is on
        if self.ui.lockObjectSizes.isChecked():
            if par=='r':
                r=self.input[obj]['r'].text()
                if obj=='target':
                    self.input['comparison']['r'].setText(r)
                elif obj=='comparison':
                    self.input['target']['r'].setText(r)

        # Check if background size locking is on
        if self.ui.lockBackgroundSize.isChecked():
            if par in ['r1','r2']:
                r=self.input[obj][par].text()
                if obj=='target':
                    self.ui.cmpAnulusRadioButton.setChecked(True)
                    self.input['comparison'][par].setText(r)
                elif obj=='comparison':
                    self.ui.tgtAnulusRadioButton.setChecked(True)
                    self.input['target'][par].setText(r)
            elif par in ['x1','y1','x2','y2']:
                c=self.input[obj][par].text()
                if obj=='target':
                    self.ui.cmpRegionRadioButton.setChecked(True)
                    self.input['comparison'][par].setText(c)
                elif obj=='comparison':
                    self.ui.tgtRegionRadioButton.setChecked(True)
                    self.input['target'][par].setText(c)

        # Check if background region centering
        if self.ui.allignTgtVerticalCenter.isChecked():
            if par in ['y1','y2']:
                y=float(self.input[obj][par].text())
                center=self.img.shape[0]/2.0
                height=abs(y-center)
                self.input[obj]['y1'].setText(str(center+height))
                self.input[obj]['y2'].setText(str(center-height))
            
        # Draw markers
        self.draw(key)

        # Unset block
        self.block[obj][par]=False

    def draw(self,key):
        """Draws markers for object positions, and backgrounds.

        To be called when any input widget value changes.

        *key* is given by the signal mapper and consists of a string with
        the object and parameter separated by a comma.
        """

        # Get object and parameter from key
        obj,par=str(key).split(',')

        try:
            # Set amplifier
            self.amp[obj]=self.getCurrentAmp()

            # Draw markers
            if par=='x' or par=='y' or par=='r':
                x=float(self.input[obj]['x'].text())
                y=float(self.input[obj]['y'].text())
                r=float(self.input[obj]['r'].text())

                self.display[obj]['position']=True

                self.mark(obj,x,y,r,color=self.line[obj]['color'],lw=self.line[obj]['width'])

            elif par=='r1' or par=='r2':
                # Annulus is selected so remove region marker
                self.imdisplay.removePatch(obj+'_region')

                x=float(self.input[obj]['x'].text())
                y=float(self.input[obj]['y'].text())
                r=float(self.input[obj][par].text())

                # Keep track of the selected background mode
                self.display[obj]['annulus']=True
                self.display[obj]['region']=False

                self.mark(obj+'_'+par,x,y,r,color=self.line[obj]['color'],lw=self.line[obj]['width'])

            elif par=='x1' or par=='y1' or par=='x2' or par=='y2':
                # Region is selected so remove annulus markers
                self.imdisplay.removePatch(obj+'_r1')
                self.imdisplay.removePatch(obj+'_r2')

                x1=float(self.input[obj]['x1'].text())
                y1=float(self.input[obj]['y1'].text())
                x2=float(self.input[obj]['x2'].text())
                y2=float(self.input[obj]['y2'].text())

                # Keep track of the selected background mode
                self.display[obj]['annulus']=False
                self.display[obj]['region']=True

                self.imdisplay.addRectangle(obj+'_region',x1,y1,x2,y2,
                    color=self.line[obj]['color'],lw=self.line[obj]['width'])
            
            # Redraw canvas
            self.imdisplay.redraw_canvas(keepzoom=True)

        except ValueError:
            pass

    def redraw(self, number):
        """Redraws object and background markers for all objects on the
        currently displayed amplifier *number*.
        """

        self.imdisplay.reset()

        # Find wich amplifier is currently displayed
        amp=self.getCurrentAmp()

        # (Re)draw markers
        for obj in ['target','comparison']:
            if self.amp[obj]==amp:
                if self.display[obj]['position']:
                    self.draw(obj+','+'r')
                if self.display[obj]['annulus']:
                    self.draw(obj+','+'r1')
                    self.draw(obj+','+'r2')
                if self.display[obj]['region']:
                    self.draw(obj+','+'y2')

    def getCurrentAmp(self, namps=4):
        """Returns the currently displayed amplifier.
        
        *namps* is the number of amplifiers on the CCD.
        """

        # Get exposure number
        n=int(self.ui.imageSpinBox.value())

        # Convert exposure number to current amplifier number
        amp=n%namps
        if amp==0:
            amp=namps

        return amp

    def captureHandler(self, key):
        """Called when a capture button is clicked.

        *key* is given by the signal mapper and consists of a string with
        the object and parameter separated by a comma.

        Depending on the *key* input widgets are added to the current
        display lists.
        Subsequent calls to `self.selectionHandler` get displayed in
        the listed widgets.
        """

        # Get object and parameter from key
        obj,par=str(key).split(',')
        
        # Add input widgets to lists
        if par=='position':
            self.xdisplay=[self.input[obj]['x']]
            self.ydisplay=[self.input[obj]['y']]
            self.rdisplay=[]
        elif par=='radius':
            self.xdisplay=[]
            self.ydisplay=[]
            self.x=float(self.input[obj]['x'].text())
            self.y=float(self.input[obj]['y'].text())
            self.rdisplay=[self.input[obj]['r']]
        elif par=='annulus':
            self.xdisplay=[]
            self.ydisplay=[]
            self.x=float(self.input[obj]['x'].text())
            self.y=float(self.input[obj]['y'].text())
            self.rdisplay=[self.input[obj]['r1'], self.input[obj]['r2']]
        elif par=='region':
            self.xdisplay=[self.input[obj]['x1'], self.input[obj]['x2']]
            self.ydisplay=[self.input[obj]['y1'], self.input[obj]['y2']]
            self.rdisplay=[]

    def selectionHandler(self, x, y):
        """Event handler for click in image display window.

        *x*, *y* is the position (in image pixel coordinates) of the click.
        These positions are inserted into the first input widgets in the
        display lists.

        If a radius is requested this is calculated from the position given
        in (self.x, self.y) which should be set to the current object.
        """

        if len(self.xdisplay)>0:
            display=self.xdisplay.pop(0)
            display.setText(str(x))
        if len(self.ydisplay)>0:
            display=self.ydisplay.pop(0)
            display.setText(str(y))
        if len(self.rdisplay)>0:
            r=np.sqrt((x-self.x)**2+(y-self.y)**2)
            display=self.rdisplay.pop(0)
            display.setText(str(r))

    def setSearchDistance(self, distance):
        """Set search distance used for recentering."""
        self.distance=int(distance)

    def setMarkerType(self, marker):
        """Set marker type to 'circle' or 'square'."""
        if marker in ['circle','square']:
            self.mark_with=marker
        else:
            raise SaltIOError('Unknown marker type '+str(marker))

    def setLineColor(self, object, color):
        """Changes the default line color used for marking."""

        self.line[object]['color']=color

    def setLineWidth(self, object, width):
        """Changes the default line width used for marking."""

        self.line[object]['width']=width

    def save(self):
        """Save configuration.
        
        The format is::
            For objects that use an anullus:
            object amp x y r r1 r2 
            For objects that use a region:
            object amp x y r x1 y1 x2 y2
        """
        if (self.ui.tgtAnulusRadioButton.isChecked() and self.ui.cmpRegionRadioButton.isChecked()) or \
           (self.ui.tgtRegionRadioButton.isChecked() and self.ui.cmpAnulusRadioButton.isChecked()):
               msg='SLOTPREVIEW--SLOTPHOT can not handle different background types'
               raise SaltError(msg)
        

        # Write values to file
        with open(self.config,'w') as f:
            for i,obj in enumerate(['target','comparison']):
                b_type='region'
                if obj=='target':
                    print obj, self.ui.tgtAnulusRadioButton.isChecked()
                    if self.ui.tgtAnulusRadioButton.isChecked(): b_type='annulus'
                elif obj=='comparison':
                    if self.ui.cmpAnulusRadioButton.isChecked(): b_type='annulus'

                #  If r1 is not zero, assumes annulus
                line='%i\t%i\t' % (i+1, self.amp[obj]) 
                if b_type=='annulus':
                   line+=''.join('%3.2f\t' % float(self.input[obj][key].text()) for key in ['x', 'y', 'r', 'r1', 'r2'])
                else:
                   line+=''.join('%3.2f\t' % float(self.input[obj][key].text()) for key in ['x', 'y', 'r', 'x1', 'y2', 'x2', 'y2'])


                # Write string to configfile
                f.write(line.rstrip()+'\n')

        # Exit program
        self.close()

