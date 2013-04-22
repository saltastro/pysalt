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

"""slotview is an interactive tool to analyze slotmode data.
The input for the file is either the same input into slotphot or the output from slotphot
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
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

# Gui library imports
from PyQt4 import QtGui, QtCore
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg

# Salt imports
from saltgui import ImageDisplay, MplCanvas
from salterror import SaltIOError

import saltsafeio as saltio
import saltsafekey as saltkey
from saltsafelog import logging

from SlotViewWindow import SlotViewWindow

debug=True

# Make sure the plotting functions work with an older version of matplotlib
try:
    import matplotlib.pyplot as plt
except ImportError:
    import matplotlib.pylab as plt

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
# Plot the data

class makeplotdata(QtGui.QMainWindow):


    def __init__(self, struct, pid, tflux, cflux, ratio, time, phottype, sleep, vig_lo, vig_hi, \
                 tx, ty, cx, cy, r, br1, br2, naxis1, naxis2, clobber, logfile, verbose):
        """As the data is measured, plots the target and companion, the drift, both light curves and the ratio

        returns status
        """
        #set up the variables
        status=0
        maxcolumn=7
        self.struct = struct
        self.infile=struct._HDUList__file.name
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
        self.logfile=logfile
        self.clobber=clobber
        self.verbose=verbose
        self.fft=False
        self.stopplay=False
        self.sleep=sleep
        self.zbox=[]
        self.newphot=0
        self.npoint=4
        if self.phottype=='circular':
            self.npoint=24


        if status==0:
            self.id=0
            self.nframes=len(self.struct)
            self.header=self.struct[int(self.pid[self.id])].header
            self.goodframes=self.dtime*0+1

        # Setup widget
        QtGui.QMainWindow.__init__(self)

        # Set main widget
        self.main = QtGui.QWidget(self)

        # Set window title
        self.setWindowTitle("Slotview: "+self.infile)


        #self.root.bind("<Destroy>", self.destroy)
        #self.root.bind("D", self.deleteframe)
        #self.root.bind("u", self.undeleteframe)
        #self.root.bind("n", self.call_playone)
        #self.root.bind("b", self.call_revone)
        #self.root.bind("?", self.help)
        #self.root.bind("q", self.destroy)
        #self.root.bind("<Button-1>", self.callback)

        #set up the variables for which graphs to plot
        #self.ratiovar=Tk.IntVar(master=self.root, value=1)
        #self.star1var=Tk.IntVar(master=self.root, value=0)
        #self.star2var=Tk.IntVar(master=self.root, value=0)

        #self.slotfig=plt.figure(figsize=(8,1.5),dpi=72)
        #plot the  data
        #self.plotdataarray()


        #self.lcfig=plt.figure(figsize=(8,5),dpi=72)
        #plot the light curve
        #self.lcx1=self.dtime.min()
        #self.lcx2=self.dtime.max()
        #self.plotlightcurve()

        inrow=4
        lcrow=0
        pcrow=1
        darow=2
        cprow=3
        qurow=5

        #add light curve plot
        #self.lccanvas = FigureCanvasTkAgg(self.lcfig, master=self.root)
        #self.lccanvas.show()
        #self.lccanvas.get_tk_widget().grid(row = lcrow, column = 0, columnspan = maxcolumn, sticky = 'news')
        #self.lccanvas.mpl_connect('button_press_event',self.lcpickstar)
        #self.lccanvas.mpl_connect('motion_notify_event',self.lcdrawbox)
        #self.lccanvas.mpl_connect('button_release_event',self.lczoom)

        #add data array plot
        #self.canvas = FigureCanvasTkAgg(self.slotfig, master=self.root)
        #self.canvas.show()
        #self.canvas.blit()
        #self.canvas.get_tk_widget().grid(row = darow, column = 0, columnspan = maxcolumn, sticky = 'news')
        #self.canvas.mpl_connect('key_press_event',self.newphoto)

        #add the control widget
        #self.cpFrame = Tk.Frame(master=self.root)
        #self.cpFrame.grid(row=cprow, column=0, columnspan=maxcolumn, sticky='ew')

        #self.frevbutton = Tk.Button(master=self.cpFrame, text='< <', width=5, command=self.freverse)
        #self.frevbutton.grid(row=0, column=0, sticky='ew')
        #self.revbutton = Tk.Button(master=self.cpFrame, text='<',width=5,  command=self.reverse)
        #self.revbutton.grid(row=0, column=1, sticky='ew')
        #self.rev1button = Tk.Button(master=self.cpFrame, text='-',width=5,  command=self.revone)
        #self.rev1button.grid(row=0, column=2, sticky='ew')

        #self.play1button = Tk.Button(master=self.cpFrame, text='+',width=5,  command=self.playone)
        #self.play1button.grid(row=0, column=4, sticky='ew')
        #self.playbutton = Tk.Button(master=self.cpFrame, text='>',width=5,  command=self.play)
        #self.playbutton.grid(row=0, column=5, sticky='ew')
        #self.fplaybutton = Tk.Button(master=self.cpFrame, text='> >',width=5,  command=self.fplay)
        #self.fplaybutton.grid(row=0, column=6, sticky='ew')

        #self.stopbutton = Tk.Button(master=self.cpFrame, text='Stop',width=5,  command=self.stop)
        #self.stopbutton.grid(row=0, column=3, sticky='ew')

        #add the information panel
        #self.idtext= Tk.StringVar(master=self.root )
        #self.imgtext= Tk.StringVar(master=self.root )
        #self.timetext= Tk.StringVar(master=self.root )
        #self.idLabel  = Tk.Label(master=self.root, fg='#000000',textvariable=self.idtext, relief='solid')
        #self.idLabel.grid(row=inrow, column=0, sticky='ew')
        #self.imgLabel = Tk.Label(master=self.root, textvariable=self.imgtext, relief='solid')
        #self.imgLabel.grid(row=inrow, column=1, columnspan=3, sticky='ew')
        #self.timeLabel = Tk.Label(master=self.root, textvariable=self.timetext, relief='solid')
        #self.timeLabel.grid(row=inrow, column=4, columnspan=3, sticky='ew')
        #self.setinfolabels()

        #add the plot control panel
        #self.ratiobutton=Tk.Checkbutton(master=self.root, text='Flux Ratio', variable=self.ratiovar, \
        #                                command=self.calllccheck)
        #self.ratiobutton.grid(row=pcrow, column=0, sticky='ew')
        #self.star1button=Tk.Checkbutton(master=self.root, text='Star1 Flux', variable=self.star1var, \
        #                                command=self.calllccheck)
        #self.star1button.grid(row=pcrow, column=1, sticky='ew')
        #self.star2button=Tk.Checkbutton(master=self.root, text='Star2 Flux', variable=self.star2var, \
        #                                command=self.calllccheck)
        #self.star2button.grid(row=pcrow, column=2, sticky='ew')
        #self.resetbutton = Tk.Button(master=self.root, text='Reset', command=self.callreset)
        #self.resetbutton.grid(row=pcrow, column=6, sticky='ew')
        #self.savebutton = Tk.Button(master=self.root, text='save', command=self.callsave)
        #self.savebutton.grid(row=pcrow, column=5, sticky='ew')

        #add the quit button
        #self.quFrame = Tk.Frame(master=self.root)
        #self.quFrame.grid(row=qurow, column=0, columnspan=maxcolumn, sticky='ew')
        #self.exitbutton = Tk.Button(master=self.quFrame, text='Quit', command=self.exit)
        #self.exitbutton.grid(row=0, column=3, sticky='ew')

        #create the tabs
        self.tabWidget=QtGui.QTabWidget()

        #layout the widgets
        mainLayout = QtGui.QVBoxLayout(self.main)
        mainLayout.addWidget(self.tabWidget)


        # Set the main widget as the central widget
        self.setCentralWidget(self.main)

        # Destroy widget on close
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)


        return


    def runplotdata(self):

        Tk.mainloop()



    def destroy(self, e):
        self.root.quit()
        return

    def exit(self):
        self.root.quit()
        return

    def help(self, e):
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
    save button-Save the current light curve plot
    quit button-Quit the application
    Right Click-Display image corresponding to this time
    Left Click and Drag-In light curve plot, zoom in on this region
        """
        print helpmessage
        return

    def setinfolabels(self):
        """Set the text labels according to the current object displayed.
        Use the image header information if it is available.
        """
        #set the id number
        self.idtext.set(int(self.pid[self.id]))

        #set the image name
        oname=''
        try:
            oname=self.struct[int(self.pid[self.id])].header['ONAME']
            oext=self.struct[int(self.pid[self.id])].header['OEXT']
            oname=oname+'[%i]'%oext
        except Exception, e:
            try:
                oname=self.struct[0].header['OBJECT']
            except:
                pass
        oname=self.struct[0].header['OBJECT']
        self.imgtext.set(oname)

        #set the time
        try:
            utime=self.struct[int(self.pid[self.id])].header['UTC-OBS']
            self.timetext.set(utime)
        except:
            self.timetext.set('')

        return
    def calllccheck(self):
        #turn the ratio curve on and off
        if self.ratiovar.get():
            self.lightcurve.set_visible(True)
        else:
            self.lightcurve.set_visible(False)

        #turn the star1 curve on and off
        if self.star1var.get():
            self.star1curve.set_visible(True)
        else:
            self.star1curve.set_visible(False)

        #turn the star2 curve on and off
        if self.star2var.get():
            self.star2curve.set_visible(True)
        else:
            self.star2curve.set_visible(False)

        self.lcy1, self.lcy2, ylabel=self.lcylimits()
        self.light_plot.set_ylim(self.lcy1, self.lcy2)
        self.light_plot.set_ylabel(ylabel)
        self.lccanvas.draw()

    def lcylimits(self):
        """Determine the y-limts depending on what plots are selected """
        mask = (self.dtime > self.lcx1)*(self.dtime<self.lcx2)*(self.goodframes>0)
        if self.ratiovar.get():
            rarr=np.compress(mask,self.ratio)
            y1=rarr.min()
            y2=rarr.max()
            ylabel='Star1/Star2'
        else:
            if self.star2var.get() and self.star1var.get():
                cfarr=np.compress(mask,self.cflux).max()
                tfarr=np.compress(mask,self.tflux).max()
                y1=0
                y2=cfarr < tfarr and tfarr or cfarr
                ylabel='Star Flux'
            elif self.star2var.get():
                cfarr=np.compress(mask,self.cflux)
                y1=0
                y2=cfarr.max()
                ylabel='Star2 Flux'
            else:
                tfarr=np.compress(mask,self.tflux)
                y1=0
                y2=tfarr.max()
                ylabel='Star1 Flux'
        return y1, y2, ylabel

    def newphoto(self, e):
        """program to being new photometry"""
        if e.key=='c' and e.xdata and e.ydata:
            cx=e.xdata
            cy=e.ydata
            cr=self.radius
            image=self.struct[int(self.pid[self.id])].data
            cimage, cx, cy  = st.calcdrift(image, cx, cy, cr, self.naxis1, self.naxis2)
            if cx >= 0 and cy >= 0:
                self.cx[self.id]=cx
                self.cy[self.id]=cy
            self.updatedataplot()
        if e.key=='t' and e.xdata and e.ydata:
            tx=e.xdata
            ty=e.ydata
            tr=self.radius
            image=self.struct[int(self.pid[self.id])].data
            timage, tx, ty  = st.calcdrift(image, tx, ty, tr, self.naxis1, self.naxis2)
            if tx >= 0 and ty >= 0:
                self.tx[self.id]=tx
                self.ty[self.id]=ty
            self.updatedataplot()
        if e.key=='p':
            self.redophot(self.id)
            #self.updatelightcurve()
            #self.lccanvas.draw()
            self.lcfig.delaxes(self.light_plot)
            self.plotlightcurve()
            self.lccanvas.draw()
            #self.callreset()
        if e.key=='P':
            nstart=self.id+1
            nend=self.nframes-1
            self.redophot(self.id)
            self.stopplay=True

            i=nstart
            while i < nend and self.stopplay:
                image=self.struct[int(self.pid[self.id])].data
                # these may be changed
                sigdet=5
                contpix=10
                sigback=3
                driftlimit=10
                iter=3
                carray, fx,fy,status = st.finddrift(image, self.cx[i-1], self.cy[i-1], self.radius,  \
                        self.naxis1, self.naxis2, sigdet, contpix, sigback, driftlimit, iter, self.logfile)
                if fx > -1  and fy > -1:
                    if fx < self.naxis1 and fy < self.naxis2:
                        dx=self.cx[i-1]-fx
                        dy=self.cy[i-1]-fy
                        self.cx[i]=fx
                        self.cy[i]=fy
                        self.tx[i]=self.tx[i-1]-dx
                        self.ty[i]=self.ty[i-1]-dy
                    else:
                        message='Not able to perform photometry'
                        print message
                        return
                else:
                    message='Not able to perform photometry'
                    print message
                    return

                self.redophot(i)
                self.lcfig.delaxes(self.light_plot)
                self.plotlightcurve()
                self.lccanvas.draw()
                if self.dtime[i] < self.lcx1 or self.dtime[i] > self.lcx2: self.callreset()
                #self.updatelightcurve()
                #self.lccanvas.draw()
                self.root.update()
                if not self.stopplay: self.updatedataplot()
                i += 1

    def redophot(self, id):
        self.newphot=1
        self.id=id
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

        tflux, tflux_err, cflux, cflux_err, ratio, ratio_err, status = \
          st.dophot(self. phottype, image, x, y, self.r, self.br1, self.br2,  \
          gain, rdnoise, self.naxis1, self.naxis2)

        if status==0:
            self.tflux[self.id]=tflux
            self.cflux[self.id]=cflux
            self.ratio[self.id]=ratio


    def lcpickstar(self, e):
        if e.button==1 and e.xdata:
            self.id=self.findtime(e.xdata)+1
            self.updatedataplot()
        if e.button==3 and e.xdata:
            self.xt1 = e.xdata
            self.yt1 = self.lcy1

    def lcdrawbox(self, e):
        if e.button==3 and e.xdata:
            self.xt2=e.xdata
            self.yt2=self.lcy2
            xp=[self.xt1, self.xt1, self.xt2, self.xt2]
            yp=[self.yt1, self.yt2, self.yt2, self.yt1]
            if self.zbox:
                self.zbox.set_visible(False)
            self.zbox,=self.light_plot.fill(xp, yp, fc='#777777', ec='#FF0000', alpha=0.5,visible=True)
            self.lccanvas.draw()

    def lczoom(self, e):
        """Handles time axis zoom on the light curve.
        Once the 3-button is released, it will capture the new position and replot the zoomed in curve"""
        if e.button==3 and e.xdata:
            self.xt2=e.xdata
            self.yt2=self.lcy2
            if self.xt2<self.xt1:
                xtemp=self.xt1
                self.xt1=self.xt2
                self.xt2=xtemp
            self.lcx1=self.xt1
            self.lcx2=self.xt2
            self.lcy1=self.yt1
            self.lcy2=self.yt2

            if self.lcx2-self.lcx1>0:
                self.lcfig.delaxes(self.light_plot)
                self.plotlightcurve()
                if self.zbox:
                    self.zbox.set_visible(False)
                self.lccanvas.draw()

    def callsave(self):
        """Save a copy of the lc curve to a .ps file"""
        self.sroot=Tk.Tk()
        self.sroot.wm_title("Save File as:")
        TitleLabel  = Tk.Label(master=self.sroot, text='Please enter a filename for the output PS file', border=5)
        TitleLabel.grid(row=0, column=0, columnspan=2, sticky='ew')
        nameLabel  = Tk.Label(master=self.sroot, text='Filename:', relief='solid')
        nameLabel.grid(row=1, column=0, sticky='ew')
        self.nametext=Tk.StringVar(master=self.sroot)
        nameEntry = Tk.Entry(master=self.sroot, textvariable=self.nametext)
        nameEntry.grid(row=1, column=1, sticky='ew')
        nameEntry.focus_set()
        self.sroot.bind('<Return>', self._finishcallsave)

        return

    def _finishcallsave(self, e):
        status=0
        self.sroot.destroy()
        name=self.nametext.get()
        if not name: return
        if name[-3:]!='.ps': name=name+'.ps'

        #remove the file if the name already exists
        if saltio.filedoesnotexist(name,self.verbose, self.logfile):
            if self.clobber:
                os.remove(name)
            else:
                message = 'ERROR -- SALTVIEW: File ' + name + ' already exists, use clobber=y'
                status = saltprint.err(logfile,message)
                return

        #turn the red dot off in the graph
        self.light_point.set_visible(False)

        #save the figure
        self.lcfig.savefig(name)


        #turn the red dot on in the graph
        self.light_point.set_visible(True)


    def callreset(self):
        self.lcx1=self.dtime.min()
        self.lcx2=self.dtime.max()
        self.lcfig.delaxes(self.light_plot)
        self.plotlightcurve()
        self.lccanvas.draw()

    def undeleteframe(self, e):
        self.goodframes[self.id] = 1
        message='SALTPHOT:  Extension %i was undeleted' %  self.pid[self.id]
        saltprint.log(self.logfile, message, self.verbose)

    def deleteframe(self, e):
        self.newphot=1
        self.goodframes[self.id] = 0
        message='SALTPHOT:  Extension %i was deleted' % self.pid[self.id]
        saltprint.log(self.logfile, message, self.verbose)

    def callback(self, e):
        print e.x, e.y

    def stop(self):
        self.stopplay=False

    def call_playone(self, e):
        self.playone()

    def call_revone(self, e):
        self.revone()

    def playone(self):
        stopid = self.nframes-2
        if self.id < (stopid): self.id=self.id+1
        self.updatedataplot()


    def play(self):
        self.stopplay=True
        stopid = self.nframes-2
        while self.stopplay and self.id < stopid:
            self.id = self.id+1
            time.sleep(self.sleep)
            self.updatedataplot()
            self.root.update()

    def fplay(self):
        self.stopplay=True
        stopid = self.nframes-2
        while self.stopplay and self.id < stopid:
            self.id = self.id+1
            self.updatedataplot()
            self.root.update()

    def revone(self):
        if self.id > 0: self.id=self.id-1
        self.updatedataplot()


    def reverse(self):
        self.stopplay=True
        while self.stopplay and self.id > 0:
            self.id = self.id-1
            time.sleep(self.sleep)
            self.updatedataplot()
            self.root.update()

    def freverse(self):
        self.stopplay=True
        while self.stopplay and self.id > 0:
            self.id = self.id-1
            self.updatedataplot()
            self.root.update()



    def callsetfft(self, label):
        if label=='FFT':
            self.fft=(not self.fft)
            self.plotfft()

    def plotfft(self):
        fftfig=plt.figure(figsize=(8,8),dpi=72)
        axfft=fftfig.add_axes([0.10,0.10,0.8,0.50], autoscale_on=True)
        mask = (self.dtime > self.lcx1)*(self.dtime<self.lcx2)
        tarr=np.compress(mask,self.dtime)
        rarr=np.compress(mask,self.ratio)
        #ftarr=np.fft.fft(tarr)
        ftarr=np.arange(len(tarr))
        frarr=np.fft.fft(rarr)

        axfft.hold(True)
        fftcurve=axfft.plot(ftarr,frarr,linewidth=0.5,linestyle='-',marker='',color='b')
        plt.show()


    def slide_update(self, val):
        self.id=self.findtime(val)
        self.updatedataplot()





    def plotdataarray(self):
        """Plot the image array

        return axes
        """
        self.ob_plot = self.slotfig.add_axes([0.10,0.10,0.8,0.80], autoscale_on=True)
        plt.setp(plt.gca(),xticks=[],yticks=[])
        plt.jet()
        self.array=self.struct[int(self.pid[self.id])].data
        self.imarr=self.ob_plot.imshow(self.array,origin='lower')
        #Plot the apertures
        self.cbox,=self.plotbox('#00FF00',self.cx[self.id],self.cy[self.id],self.radius,self.npoint,self.naxis1, self.naxis2)
        self.tbox,=self.plotbox('#FFFF00',self.tx[self.id],self.ty[self.id],self.radius,self.npoint,self.naxis1, self.naxis2)


    def updatedataplot(self):
        """Handle updating the light curve plot and the data array plot when the
        data array image is changed
        """

        #update the information panel
        self.setinfolabels()

        self.ptime=self.dtime[self.id]
        self.pratio=self.ratio[self.id]
        #Check to make the red button hasn't moved outside of the plotting area
        if self.ptime < self.lcx1 or self.ptime > self.lcx2: self.callreset()

        #update the red piont on the light curve plot
        self.light_point.set_data(np.asarray([self.ptime]), np.asarray([self.pratio]))
        #self.lcfig.delaxes(self.light_plot)
        #self.plotlightcurve()
        self.lccanvas.draw()

        #update the data plot
        self.array=self.struct[int(self.pid[self.id])].data
        self.imarr.set_data(self.array)
        #Plot the apertures
        self.updateplotbox(self.cbox,'#00FF00',self.cx[self.id],self.cy[self.id],self.radius,self.npoint,self.naxis1, self.naxis2)
        self.updateplotbox(self.tbox,'#FFFF00',self.tx[self.id],self.ty[self.id],self.radius,self.npoint,self.naxis1, self.naxis2)
        self.canvas.draw()


    def updateplotbox(self,box, bcolor,x,y,r,npoint, nx1,nx2):
        apx,apy=self.makeplotbox(x,y,r,npoint,nx1,nx2)
        box.set_data(apx, apy)
        box.set_color(bcolor)


    def makeplotbox(self,x,y,r,npoint, nx1,nx2):
        apx=np.zeros(npoint+1)
        apy=np.zeros(npoint+1)
        for i in range(npoint+1):
            theta=math.pi/4+2*i*math.pi/npoint
            apx[i]=st.checkedge(x+math.cos(theta)*r,0,nx1)
            apy[i]=st.checkedge(y+math.sin(theta)*r,0,nx2)
        return apx, apy

    def plotbox(self,bcolor,x,y,r,npoint, nx1,nx2):
        apx,apy=self.makeplotbox(x,y,r, npoint, nx1,nx2)
        return self.ob_plot.plot(apx,apy,ls='-',color=bcolor,lw=2,marker='')

    def updatelightcurve(self):
        mask = (self.dtime > self.lcx1)*(self.dtime<self.lcx2)*(self.goodframes>0)
        tarr=np.compress(mask,self.dtime)
        rarr=np.compress(mask,self.ratio)
        self.lightcurve.set_data(tarr, rarr)
        self.ptime=self.dtime[self.id]
        self.pratio=self.ratio[self.id]
        self.light_point.set_data(np.asarray([self.ptime]), np.asarray([self.pratio]))

    def plotlightcurve(self):
        """Plot the light curve
        return ax
        """
        mask = (self.dtime > self.lcx1)*(self.dtime<self.lcx2)*(self.goodframes>0)
        tarr=np.compress(mask,self.dtime)
        rarr=np.compress(mask,self.ratio)
        tfarr=np.compress(mask,self.tflux)
        cfarr=np.compress(mask,self.cflux)

        self.lcy1,self.lcy2, ylabel=self.lcylimits()

        self.light_plot = self.lcfig.add_axes([0.10,0.10,0.8,0.80], autoscale_on=False, adjustable='datalim'  )
        self.light_plot.hold(True)

        #plot the curve
        self.lightcurve,=self.light_plot.plot(tarr,rarr,linewidth=0.5,linestyle='-',marker='',color='b')
        if self.ratiovar.get():
            self.lightcurve.set_visible(True)
        else:
            self.lightcurve.set_visible(False)

        #plot the flux curve for star 1
        self.star1curve,=self.light_plot.plot(tarr,tfarr,linewidth=0.5,linestyle='-',marker='',color='y')
        if  self.star1var.get():
            self.star1curve.set_visible(True)
        else:
            self.star1curve.set_visible(False)

        #plot the flux curve for star 1
        self.star2curve,=self.light_plot.plot(tarr,cfarr,linewidth=0.5,linestyle='-',marker='',color='g')
        if self.star2var.get():
            self.star2curve.set_visible(True)
        else:
            self.star2curve.set_visible(False)

        #plot a point which matches the time
        self.ptime=self.dtime[self.id]
        self.pratio=self.ratio[self.id]
        self.light_point,=self.light_plot.plot(np.asarray([self.ptime]), np.asarray([self.pratio]), linestyle='', marker='o', mec='#FF0000', mfc='#FF0000')


        self.light_plot.set_xlim(self.lcx1, self.lcx2)
        self.light_plot.set_ylim(self.lcy1, self.lcy2)
        self.light_plot.set_ylabel(ylabel)
        self.light_plot.set_xlabel('Time (s)')




    def findtime(self, x):
        arr=abs(self.dtime-x)
        return arr.argmin()

    def onselect(self, vmin, vmax):
        pass

    def release(self, event):
        ax=event.inaxes
        if ax==self.light_plot and abs(event.xdata-self.xt)>self.min_xlim:
            if (event.xdata > self.xt):
                self.lcx1=self.xt
                self.lcx2=event.xdata
            else:
                self.lcx2=self.xt
                self.lcx1=event.xdata

            self.plotlightcurve()
            del self.stime
            self.stime = plt.Slider(self.axslid, 'Time', self.lcx1, self.lcx2, valinit=1,  valfmt='%1.2f')
            self.stime.on_changed(self.slide_update)


    def onpress(self, event):
        if event.button!=1: return
        ax=event.inaxes
        if ax==self.light_plot:
            self.xt = event.xdata
            self.stime.set_val(self.xt)
            self.id=self.findtime(self.xt)
            self.updatedataplot()



# -----------------------------------------------------------
# main code

parfile = iraf.osfn("slottools$slotview.par")
t = iraf.IrafTaskFactory(taskname="slotview",value=parfile,function=slotview, pkgname='slottools')
