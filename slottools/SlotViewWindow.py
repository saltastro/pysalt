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
"""SlotViewWindow is the PyQt4 GUI for viewing and analyzing light curves
   produced by slotphot.

Updates:

20110516
    * Split this section off from the main code of slotphot.py
    * Updated it to PyQt4

Limitations:
*  Animation are currently not handled.  They need to be updated to be
handled in a PyQt4 friendly manner

"""

# Ensure Python 2.5 compatibility
from __future__ import with_statement

# General imports
import os, sys, time
import numpy as np
import pyfits
from pyraf import iraf
from pyraf.iraf import pysalt

# Gui library imports
from PyQt4 import QtGui, QtCore
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT
from matplotlib.figure import Figure


# Salt imports
import saltsafeio
from saltgui import ImageDisplay, MplCanvas
from salterror import SaltError, SaltIOError

#slottools
import slottool as st

class SlotViewWindow(QtGui.QMainWindow):
   """Main application window."""

   def __init__(self, struct, pid, tflux, cflux, ratio, time, phottype, sleep,  
                tx, ty, cx, cy, r, br1, br2, naxis1, naxis2, sigdet, contpix, driftlimit):
        """Default constructor."""
        maxcolumn=7
        self.struct = struct
        self.infile=struct._HDUList__file.name
        self.name=self.struct[0].header['OBJECT']
        self.pid=pid
        self.dtime=time.copy()
        self.tflux=tflux
        self.cflux=cflux
        self.ratio=ratio
        self.min_xlim=10
        self.radius=r['comparison']
        self.r=r
        self.br1=br1
        self.br2=br2
        self.tx=tx
        self.ty=ty
        self.cx=cx
        self.cy=cy
        self.phottype=phottype
        self.naxis1=naxis1
        self.naxis2=naxis2
        self.sigdet=sigdet
        self.contpix=contpix
        self.driftlimit=driftlimit
        self.niter=5
        self.sigback=5
        self.fft=False
        self.stopplay=False
        self.sleep=sleep
        self.zbox=[]
        self.npoint=4
        self.id=0
        self.nframes=len(self.struct)
        self.header=self.struct[int(self.pid[self.id])].header
        self.goodframes=self.dtime*0+1
        if self.phottype=='circular': self.npoint=24


        # Setup widget
        QtGui.QMainWindow.__init__(self)

        # Set main widget
        self.main = QtGui.QWidget(self)

        # Set window title
        self.setWindowTitle("SlotView %s" % self.infile)

        #set up the different pages
        self.slotPage=QtGui.QWidget()
        
        #set up the differen panels
        self.set_optionpanel()
        self.set_imagepanel()
        self.set_controlpanel()
        self.set_plotpanel()
        self.set_infopanel()

        # Set up the layout
        slotLayout = QtGui.QVBoxLayout(self.slotPage)
        slotLayout.addWidget(self.plotpanel)
        slotLayout.addWidget(self.optipanel)
        slotLayout.addWidget(self.imdisplay)
        slotLayout.addWidget(self.contpanel)
        slotLayout.addWidget(self.infopanel)


        #create the tabs
        #self.tabWidget=QtGui.QTabWidget()
        #self.tabWidget.addTab(self.slotPage, 'Slot')
 
        #layout the widgets
        mainLayout = QtGui.QVBoxLayout(self.main)
        mainLayout.addWidget(self.slotPage)
        #mainLayout.addWidget(self.tabWidget)
        #self.setLayout(mainLayout)


        # Set focus to main widget
        self.main.setFocus()

        # Set the main widget as the central widget
        self.setCentralWidget(self.main)

        # Destroy widget on close
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)

        # Close when config dialog is closed
        self.connect(self.main, QtCore.SIGNAL('keyPressEvent'), self.keyPressEvent)
        #self.connect(self.conf, QtCore.SIGNAL('destroyed()'), self, QtCore.SLOT('close()'))
        #self.connect(self.tabWidget, QtCore.SIGNAL('currentChanged(int)'), self.currentChanged)
        #self.connect(self.imagePage, QtCore.SIGNAL('regionChange(int,int)'), self.regionChange)
        #self.connect(self.imagePage, QtCore.SIGNAL('runauto(int, int, int)'), self.runauto)

   def keyPressEvent(self, event):
       #print "Key Pressed:", event.key()
       self.keyEvents(str(event.text()))

   def keyEvents(self, key, x=None, y=None):
       if key=='?':
           self.help()
       if key=='q':
           self.close()
       if key=='d':
           self.goodframes[self.id] = 0
           self.updatepage()
       if key=='u':
           self.goodframes[self.id] = 1
       if key=='p':
           self.redophot(self.id)
           self.updatedataplot()
       if key=='P':
           self.thread=QtCore.QThread()
           self.thread.run=self.newphot
           self.thread.start()
       if key=='t' and x is not None and y is not None:
           tr=self.radius
           imarr=self.struct[int(self.pid[self.id])].data
           timage, tx, ty  = st.calcdrift(imarr, x, y, tr, self.naxis1, self.naxis2)
           if tx >= 0 and ty >= 0:
                self.tx[self.id]=tx
                self.ty[self.id]=ty
           self.updatepage()
       if key=='c' and x is not None and y is not None:
           r=self.radius
           imarr=self.struct[int(self.pid[self.id])].data
           cimage, cx, cy  = st.calcdrift(imarr, x, y, r, self.naxis1, self.naxis2)
           if cx >= 0 and cy >= 0:
                self.cx[self.id]=cx
                self.cy[self.id]=cy
           self.updatepage()


   def onKeyPress(self, event):
       self.keyEvents(event.key, event.xdata, event.ydata)

   def onButtonPress(self, event):
       if event.button==2:
           self.lcpickframe(event)

   def currentChanged(self, event):
       #print event
       pass

   def redophot(self, i):
       """Redo the photometry for a single frame"""
       self.id=i
       x={}
       y={}
       x['target']=self.tx[self.id]
       y['target']=self.ty[self.id]
       x['comparison']=self.cx[self.id]
       y['comparison']=self.cy[self.id]
       image=self.struct[int(self.pid[self.id])].data

       #these will need to be changed
       gain=1
       rdnoise=1
       verbose=False

       try:
           tflux, tflux_err, cflux, cflux_err, ratio, ratio_err = \
               st.dophot(self. phottype, image, x, y, self.r, self.br1, self.br2,  \
               gain, rdnoise, self.naxis1, self.naxis2)
       except:
           msg="SLOTVIEW--ERROR:  Could not measure photometry in frame %i" % i
           raise SaltError(msg)
       
       self.tflux[self.id]=tflux
       self.cflux[self.id]=cflux
       self.ratio[self.id]=ratio

   def newphot(self):
       self.redophot(self.id)
       i=self.id+1
       self.stopplay=True
       while i < self.nframes-1 and self.stopplay:
           self.id=i
           imarr=self.struct[int(self.pid[self.id])].data
           carray, fx,fy = st.finddrift(imarr, self.cx[self.id-1], self.cy[self.id-1], self.radius,  \
                self.naxis1, self.naxis2, self.sigdet, self.contpix, self.sigback, self.driftlimit, self.niter)
           if 0 <= fx < self.naxis1 and 0 <= fy < self.naxis2:
               dx=self.cx[i-1]-fx
               dy=self.cy[i-1]-fy
               self.cx[i]=fx
               self.cy[i]=fy
               self.tx[i]=self.tx[i-1]-dx
               self.ty[i]=self.ty[i-1]-dy
           self.redophot(i)
           i = i+1
       print 'Stopped at', i
       self.updatepage()


   def changefluxplot(self, event):
       self.fluxplot=event
       self.updatedataplot(save_zoom=False)

   def changetstarplot(self, event):
       self.tstarplot=event
       self.updatedataplot(save_zoom=False)

   def changecstarplot(self, event):
       self.cstarplot=event
       self.updatedataplot(save_zoom=False)

   def plotlightcurve(self):
       """Plot the light curve"""
 
       #cut the data
       self.make_lcdata()

       #make the figure
       self.light_plot=self.lccanvas.figure.add_axes([0.10,0.15,0.8,0.80], autoscale_on=False, adjustable='datalim'  )
       self.light_plot.hold(True)

       #plot the curve
       self.lightcurve,=self.light_plot.plot(self.tarr,self.rarr,linewidth=0.5,linestyle='-',marker='',color='b')
       if self.fluxplot:
           self.lightcurve.set_visible(True)
       else:
           self.lightcurve.set_visible(False)

       #plot the flux curve for star 1
       self.tstarcurve,=self.light_plot.plot(self.tarr,self.tfarr,linewidth=0.5,linestyle='-',marker='',color='y')
       if  self.tstarplot:
           self.tstarcurve.set_visible(True)
       else:
           self.tstarcurve.set_visible(False)

       #plot the flux curve for star 1
       self.cstarcurve,=self.light_plot.plot(self.tarr,self.cfarr,linewidth=0.5,linestyle='-',marker='',color='g')
       if self.cstarplot:
           self.cstarcurve.set_visible(True)
       else:
           self.cstarcurve.set_visible(False)

       #plot a point which matches the time
       self.ptime=self.dtime[self.id]
       self.pratio=self.ratio[self.id]
       self.light_point,=self.light_plot.plot(np.asarray([self.ptime]), np.asarray([self.pratio]), linestyle='', marker='o', mec='#FF0000', mfc='#FF0000')


       self.find_lclimits()
       ylabel='Star1/Star2'
       self.light_plot.set_ylabel(ylabel)
       self.light_plot.set_xlabel('Time (s)')
 
   def make_lcdata(self):
       #cut the data
       mask = (self.goodframes>0)
       self.tarr=np.compress(mask,self.dtime)
       self.rarr=np.compress(mask,self.ratio)
       self.tfarr=np.compress(mask,self.tflux)
       self.cfarr=np.compress(mask,self.cflux)

   def find_lclimits(self, save_zoom=False):
       """Find the limits on the Light Curve plot"""
       if save_zoom: return
       self.lcx1=self.tarr.min()
       self.lcx2=self.tarr.max()
       self.light_plot.set_xlim(self.lcx1, self.lcx2)
       #determine the minimum y value based on what plots are turned on
       if self.fluxplot:
           self.lcy1=self.rarr.min()
           self.lcy2=self.rarr.max()
       if self.tstarplot:
           self.lcy1=self.tfarr.min()
           self.lcy2=self.tfarr.max()
       if self.cstarplot:
           self.lcy1=self.cfarr.min()
           self.lcy2=self.cfarr.max()
       if self.fluxplot and self.tstarplot and self.cstarplot:
           self.lcy1=min(self.rarr.min(), self.tfarr.min(), self.cfarr.min())
           self.lcy2=max(self.rarr.max(), self.tfarr.max(), self.cfarr.max())
       if self.tstarplot and self.cstarplot and not self.fluxplot:
           self.lcy1=min(self.tfarr.min(), self.cfarr.min())
           self.lcy2=max(self.tfarr.max(), self.cfarr.max())
       if self.tstarplot and not self.cstarplot and  self.fluxplot:
           self.lcy1=min(self.tfarr.min(), self.rarr.min())
           self.lcy2=max(self.tfarr.max(), self.rarr.max())
       if not self.tstarplot and self.cstarplot and self.fluxplot:
           self.lcy1=min(self.rarr.min(), self.cfarr.min())
           self.lcy2=max(self.rarr.max(), self.cfarr.max())
           
       self.light_plot.set_ylim(self.lcy1, self.lcy2)

   def lcpickframe(self, event):
       self.set_id(event.xdata)
       self.updatepage()

   def set_id(self, t):
       """Given a time, set the object id"""
       self.id=np.abs(self.dtime-t).argmin()
           

   def set_plotpanel(self, hmin=250):
       #set up the control panel
       self.plotpanel=QtGui.QWidget()

       self.lccanvas=MplCanvas()
       self.plotlightcurve()
       #add the actions
       self.lccanvas.mpl_connect('button_press_event',self.onButtonPress)

       # Add navigation toolbars for each widget to enable zooming
       self.toolbar=NavigationToolbar2QT(self.lccanvas,self)

       # Set up the layout
       plotLayout = QtGui.QVBoxLayout(self.plotpanel)
       plotLayout.addWidget(self.lccanvas)
       plotLayout.addWidget(self.toolbar)

   def set_optionpanel(self):
       #set up the control panel
       self.optipanel=QtGui.QWidget()

       #set up the options
       self.fluxplot=True
       self.tstarplot=False
       self.cstarplot=False
       self.fluxButton = QtGui.QCheckBox("Flux Ratio")
       self.fluxButton.setChecked(self.fluxplot)
       self.fluxButton.clicked.connect(self.changefluxplot)
       self.tstarButton = QtGui.QCheckBox("Target")
       self.tstarButton.clicked.connect(self.changetstarplot)
       self.cstarButton = QtGui.QCheckBox("Comparison")
       self.cstarButton.clicked.connect(self.changecstarplot)

       # Set up the layout
       optiLayout=QtGui.QGridLayout(self.optipanel)
       optiLayout.addWidget(self.fluxButton, 0, 0, 1,1)
       optiLayout.addWidget(self.tstarButton, 0, 1, 1,1)
       optiLayout.addWidget(self.cstarButton, 0, 2, 1,1)
       #optiLayout.addWidget(self.imagtoolbar)

   def set_imagepanel(self, name=None, cmap='gray', scale='zscale', contrast=0.1):
       """set up the Image Panel"""
       hmin=150
       wmin=400
       #set up the control panel
       self.imagpanel=QtGui.QWidget()

       #hmin=wmin*self.naxis2/self.naxis1
       #print self.naxis1, self.naxis2, wmin, hmin

       # Add FITS display widget with mouse interaction and overplotting
       self.imdisplay = ImageDisplay()
       #self.imdisplay.setMinimumWidth(wmin)

       # Set colormap
       self.imdisplay.setColormap(cmap)

       # Set scale mode for dynamic range
       imarr=self.struct[int(self.pid[self.id])].data
       self.imdisplay.scale=scale
       self.imdisplay.contrast=contrast
       self.imdisplay.aspect='auto'
       self.imdisplay.loadImage(imarr)
       #self.imdisplay.drawImage()
       hmin=self.imdisplay.width()*self.naxis2/self.naxis1
       self.imdisplay.setMaximumHeight(hmin)
       self.imdisplay.setMinimumHeight(hmin)

       #add the rectangles
       self.add_mark(self.tx[self.id], self.ty[self.id], 'target',color='b', lw=2)
       self.add_mark(self.cx[self.id], self.cy[self.id], 'comparison',color='g', lw=2)
       self.imdisplay.redraw_canvas()
       self.imdisplay.axes.set_xticklabels([])
       self.imdisplay.axes.set_yticklabels([])

       self.imdisplay.connectMatplotlibMouseMotion()
       self.imdisplay.mpl_connect('button_press_event', self.onButtonPress)
       self.imdisplay.mpl_connect('key_press_event', self.onKeyPress)


       # Add navigation toolbars for each widget to enable zooming
       self.imagtoolbar=NavigationToolbar2QT(self.imdisplay,self)

       # Set up the layout
       imagLayout = QtGui.QVBoxLayout(self.imagpanel)
       #imagLayout.addWidget(self.imdisplay)
       imagLayout.addWidget(MplCanvas())
       imagLayout.addWidget(self.imagtoolbar)

   def add_mark(self, x1, y1, label, color='g', lw=2):
       """Add the rectangle for the object"""
       self.imdisplay.removePatch(label)
       r1=self.r[label]
       rect=self.imdisplay.addSquare(label, x1, y1, r1, color=color, lw=lw)

   def set_controlpanel(self):
       """set up the Control Panel"""
       #set up the control panel
       self.contpanel=QtGui.QWidget()

       #set up the buttons
       self.frevButton = QtGui.QPushButton("<<")
       self.frevButton.clicked.connect(self.freverse)
       self.revButton = QtGui.QPushButton("<")
       self.revButton.clicked.connect(self.reverse)
       self.rev1Button = QtGui.QPushButton("-")
       self.rev1Button.clicked.connect(self.revone)
       self.stopButton = QtGui.QPushButton("Stop")
       self.stopButton.clicked.connect(self.stop)
       self.play1Button = QtGui.QPushButton("+")
       self.play1Button.clicked.connect(self.playone)
       self.playButton = QtGui.QPushButton(">")
       self.playButton.clicked.connect(self.play)
       self.fplayButton = QtGui.QPushButton(">>")
       self.fplayButton.clicked.connect(self.fplay)


       #set up the info panel layout
       contLayout=QtGui.QGridLayout(self.contpanel)
       #contLayout.addWidget(self.frevButton, 0, 0, 1, 1)
       #contLayout.addWidget(self.revButton,  0, 1, 1, 1)
       contLayout.addWidget(self.rev1Button, 0, 2, 1, 1)
       contLayout.addWidget(self.stopButton, 0, 3, 1, 1)
       contLayout.addWidget(self.play1Button,0, 4, 1, 1)
       #contLayout.addWidget(self.playButton, 0, 5, 1, 1)
       #contLayout.addWidget(self.fplayButton,0, 6, 1, 1)

   def set_infopanel(self):
       """Set up the information panel"""

       #set up the information panel
       self.infopanel=QtGui.QWidget()
        
       #add the name of the file
       self.IDValueLabel = QtGui.QLabel("%i" % self.pid[self.id])
       self.IDValueLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken )
       self.NameValueLabel = QtGui.QLabel("%s" % self.name)
       self.NameValueLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken )
       self.timeValueLabel = QtGui.QLabel("%s" % self.get_time())
       self.timeValueLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken )
  
       #set up the info panel layout
       infoLayout=QtGui.QGridLayout(self.infopanel)
       infoLayout.addWidget(self.IDValueLabel, 0, 0, 1, 1)
       infoLayout.addWidget(self.NameValueLabel, 0, 1, 1, 1)
       infoLayout.addWidget(self.timeValueLabel, 0, 2, 1, 1)
   
   def get_time(self):
       #set the time
       try:
            utime=self.struct[int(self.pid[self.id])].header['UTC-OBS']
       except:
            utime=''
       return utime


   def help(self):
        """Print the help message and the key-bindings available to the user"""
        helpmessage="""
    The following commands are available to the user:
    ? - Print this information           q - quit the viewer
    n - Move to the next image           b - move back an image
    D - Delete this image                u - undelete this image
    p - Perform photometry on this image
    P - Perform photometry starting at this image
    stop button-Stop photometry or display
    reset button-Reset the light curve plot
    Middle Click-Display image corresponding to this time
        """
        print helpmessage
        return

   def stop(self):
        self.stopplay=False
        print self.stopplay


   def playone(self):
        stopid = self.nframes-2
        if self.id < (stopid): self.id=self.id+1
        self.updatepage()


   def play(self):
        self.stopplay=True
        stopid = self.nframes-2
        while self.stopplay and self.id < stopid:
            self.id = self.id+1
            time.sleep(self.sleep)
            self.updatepage()

   def fplay(self):
        self.stopplay=True
        stopid = self.nframes-2
        while self.stopplay and self.id < stopid:
            self.id = self.id+1
            self.updatepage()

   def revone(self):
        if self.id > 0: self.id=self.id-1
        self.updatepage()


   def reverse(self):
        self.stopplay=True
        while self.stopplay and self.id > 0:
            self.id = self.id-1
            time.sleep(self.sleep)
            self.updatepage()

   def freverse(self):
        self.stopplay=True
        while self.stopplay and self.id > 0:
            self.id = self.id-1
            self.updatepage()

    
   def updatepage(self):
       """Update all the values on the page that need updating"""
       self.IDValueLabel.setText("%i" % self.pid[self.id])
       self.timeValueLabel.setText("%s" % self.get_time())

       #update the image
       imarr=self.struct[int(self.pid[self.id])].data
       self.imdisplay.loadImage(imarr)
       
       #update the boxes
       self.add_mark(self.tx[self.id], self.ty[self.id], 'target',color='b', lw=2)
       self.add_mark(self.cx[self.id], self.cy[self.id], 'comparison',color='g', lw=2)
       self.imdisplay.redraw_canvas()
       self.imdisplay.axes.set_xticklabels([])
       self.imdisplay.axes.set_yticklabels([])
 
       self.updatedataplot()

   def updatedataplot(self, save_zoom=True):
       """Update the data plot for changes in the options
       """
       #redraw the lines
       self.make_lcdata()
       self.lightcurve.set_xdata(self.tarr)
       self.lightcurve.set_ydata(self.rarr)
       self.tstarcurve.set_xdata(self.tarr)
       self.tstarcurve.set_ydata(self.tfarr)
       self.cstarcurve.set_xdata(self.tarr)
       self.cstarcurve.set_ydata(self.cfarr)
 
       #move the point
       self.light_point.set_xdata([self.dtime[self.id]])
       self.light_point.set_ydata([self.ratio[self.id]])

       if self.fluxplot:
           self.lightcurve.set_visible(True)
       else:
           self.lightcurve.set_visible(False)
       #plot the flux curve for star 1
       if  self.tstarplot:
           self.tstarcurve.set_visible(True)
       else:
           self.tstarcurve.set_visible(False)

       #plot the flux curve for star 1
       if self.cstarplot:
           self.cstarcurve.set_visible(True)
       else:
           self.cstarcurve.set_visible(False)
       self.find_lclimits(save_zoom=save_zoom)
       self.lccanvas.draw()
