"""
imageWidget is a Qt4 Widget for displaying an image in the frame
"""
import os
from PyQt4 import QtGui, QtCore
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg
from saltgui import ImageDisplay, MplCanvas


class ImageWidget(QtGui.QWidget):
   def __init__(self, hdu,  hmin=150, wmin=400, cmap='gray', scale='zscale', contrast=0.1, parent=None):
        super(imageWidget, self).__init__(parent)

        self.imdisplay = ImageDisplay()
        self.imdisplay.setMinimumHeight(hmin)
        self.imdisplay.setMinimumWidth(wmin)


        # Add FITS display widget with mouse interaction and overplotting

        # Set colormap
        self.imdisplay.setColormap(cmap)

        # Set scale mode for dynamic range
        self.imdisplay.scale=scale
        self.imdisplay.contrast=contrast
        self.imdisplay.aspect='auto'

        self.loadimage(hdu)

        # Add navigation toolbars for each widget to enable zooming
        self.toolbar=NavigationToolbar2QTAgg(self.imdisplay,self)

        #set up the information panel
        self.infopanel=QtGui.QWidget()

        #add the name of the file
        self.NameLabel = QtGui.QLabel("Filename:")
        self.NameLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised )
        self.NameValueLabel = QtGui.QLabel(self.name)
        self.NameValueLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken )


        #set up the info panel layout
        infoLayout=QtGui.QGridLayout(self.infopanel)
        infoLayout.addWidget(self.NameLabel, 0, 0, 1, 1)
        infoLayout.addWidget(self.NameValueLabel, 0, 1, 1, 5)

        # Set up the layout
        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addWidget(self.imdisplay)
        mainLayout.addWidget(self.toolbar)
        mainLayout.addWidget(self.infopanel)
        self.setLayout(mainLayout)

   def loadimage(self, hdu, ext=1):
        """Load a new image"""
        if hdu is None:
           self.name=''
           return
        self.hdu=hdu
        try:
           self.imarr=self.hdu[ext].data
        except:
           raise SaltError("Could not open the data file for %s" % self.image)
        self.imdisplay.loadImage(self.imarr)
        self.imdisplay.drawImage()
        self.name=os.path.basename(hdu._HDUList__file.name)



