import numpy as np
import os, errno
from PyQt4 import QtGui,QtCore
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg
from pyraf import iraf
import pyfits
from iraf import pysalt
from saltgui import MplCanvas
from erc_ring import fit_rings
from saltflat import saltflat
from saltfp import saltfpprep, saltfpmask
from saltfpprep import saltfpprep
from saltfpmask import saltfpmask
#from fptools import findrings, findcenter

class FpRingWidget (QtGui.QWidget):

    def __init__(self,filenumber,flatnumber,parent=None):
        super(FpRingWidget,self).__init__(parent)
        
        self.filenumber=filenumber
        self.flatnumber=flatnumber


        #set up the ring plot
        self.ringplot=MplCanvas()

        self.axes=self.ringplot.figure.add_subplot(211)
        self.erraxes=self.ringplot.figure.add_subplot(212)

        self.ringplot.setMinimumHeight(500)
#        self.ringplot.setMinimumHeight(500)

        self.loadfpring()
        self.plotRing()

#        self.redraw_canvas()

       #set up the information panel
        self.infopanel=QtGui.QWidget()

        # add a label:
        self.NameLabel = QtGui.QLabel("File number:")

       #add the name of the file
        self.NameValueLabel = QtGui.QLineEdit(str(self.filenumber))

        # and a button to process the new ring
        self.ringButton = QtGui.QPushButton('Load new ring')
        self.ringButton.clicked.connect(self.redrawRing)

        # add a label for the flat field:
        self.flatLabel = QtGui.QLabel("Flat file number:")

       #add the name of the file
        self.flatValueLabel = QtGui.QLineEdit(str(self.flatnumber))


       #set up info panel layout

        infoLayout=QtGui.QGridLayout(self.infopanel)
        infoLayout.addWidget(self.NameLabel,0,0,1,1)
        infoLayout.addWidget(self.NameValueLabel,0,1,1,1)
        infoLayout.addWidget(self.ringButton,0,2,2,1)
        infoLayout.addWidget(self.flatLabel,0,3,1,1)
        infoLayout.addWidget(self.flatValueLabel,0,4,1,1)
#        fitLayout=QtGui.QGridLayout(self.infopanel)


        # add a panel to display the fit results

        self.fitpanel=QtGui.QWidget()

        self.fitLabel = QtGui.QLabel("Fit Results:")

        fitXresult="X:" + str(self.etalon_x)
        fitYresult="Y:" + str(self.etalon_y)
        fitZresult="Z:" + str(self.etalon_z)
        fitRresult="R: %.1f" % float(self.pars['R'][0])
        fitAmpresult="Amplitude: %.1f" % float(self.pars['Amplitude'][0])
        fitRmsresult="RMS: %.1f" % float(self.rms)
        fitGammaresult="Gamma: %.1f" % float(self.pars['Gamma'][0])
        fitFWHMresult="FWHM: %.1f" % float(self.pars['FWHM'][0])


        #add the text to the fit results panel
        self.fitX = QtGui.QLabel(fitXresult)
        self.fitY = QtGui.QLabel(fitYresult)
        self.fitZ = QtGui.QLabel(fitZresult)
        self.fitR = QtGui.QLabel(fitRresult)
        self.fitAmp = QtGui.QLabel(fitAmpresult)
        self.fitRms = QtGui.QLabel(fitRmsresult)
        self.fitGamma = QtGui.QLabel(fitGammaresult)
        self.fitFWHM = QtGui.QLabel(fitFWHMresult)

        # lay them out nicely...

        fitLayout=QtGui.QGridLayout(self.fitpanel)
        fitLayout.addWidget(self.fitLabel,0,0,1,1)
        fitLayout.addWidget(self.fitX,0,1,1,1)
        fitLayout.addWidget(self.fitY,0,2,1,1)
        fitLayout.addWidget(self.fitZ,0,3,1,1)
        fitLayout.addWidget(self.fitR,0,4,1,1)
        fitLayout.addWidget(self.fitAmp,0,5,1,2)
        fitLayout.addWidget(self.fitRms,0,7,1,1)
        fitLayout.addWidget(self.fitGamma,0,8,1,1)
        fitLayout.addWidget(self.fitFWHM,0,9,1,1)
 

       # Set up the main layout
        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addWidget(self.infopanel)
        mainLayout.addWidget(self.ringplot)
        mainLayout.addWidget(self.fitpanel)
        self.setLayout(mainLayout)

    def loadfpring(self):

        date = os.getcwd().split('/')[-1]
        
        fits = "mbxgpP%s%04d.fits" % (date, self.filenumber)
        if self.flatnumber != 0:
            flat = "mbxgpP%s%04d.fits" % (date, self.flatnumber)
        else:
            flat="/home/erc/FP_utils/DefaultNeFlat.fits"

        ffits='f'+fits
        fpfile='p'+ffits
        maskedfile='m'+fpfile
        

        if (not os.path.isfile(maskedfile)): 
             saltflat(fits,ffits,'',flat, clobber='yes',verbose='no')
             
             saltfpprep(ffits,fpfile,'',clobber='yes',verbose='no')
             
             saltfpmask(fpfile,maskedfile,'',axc=798,ayc=503,arad=450,clobber='yes', verbose='no')

        self.good, self.rsq, self.prof, self.fit, self.pars = fit_rings(maskedfile)
        self.resid = self.prof - self.fit
        self.rms = self.resid.std()

        self.saveFitToFile()
        
        return
    

    def plotRing(self):
        #set up the plots....

        self.axes.plot(self.rsq, self.prof, label="Profile")
        self.axes.plot(self.rsq, self.fit, label="Fit")
        self.axes.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=2, mode="expand", borderaxespad=0.)
        self.erraxes.plot(self.rsq, self.resid, label="Profile - Fit")
        self.erraxes.legend(loc=1)
        self.show()

        return

    def redrawRing(self):
        # read the new filenumber from the input box
        self.filenumber=int(self.NameValueLabel.text())
        # recalculate and re-draw
        self.loadfpring()
        self.axes.clear()
        self.erraxes.clear()
        self.plotRing()
        #Force a redraw!!!
        self.ringplot.draw()


    def getFitsHeader(self):

        date = os.getcwd().split('/')[-1]
        fits = "mbxgpP%s%04d.fits" % (date, self.filenumber)
        hdu = pyfits.open(fits)
        (data, header) = (hdu[0].data, hdu[0].header)
        etalon = int(header['ET-STATE'].split()[3])
        hetalon_x = "ET%dX" % etalon
        hetalon_y = "ET%dY" % etalon
        hetalon_z = "ET%dZ" % etalon

        self.etalon_x=int(header[hetalon_x])
        self.etalon_y=int(header[hetalon_y])
        self.etalon_z=int(header[hetalon_z])
  
        return

    def saveFitToFile(self):

        pars=self.pars
        self.getFitsHeader()

        if os.path.isfile('outparams'): 
            try:
                outparams=open('outparams','a')
                
            except IOError, e:
                print 'Failed to open the ring file'


        else:
            try:
                outparams=open('outparams','w+')
                
            except IOError, e:
                print 'Failed to open the ring file'

            outparams.write("File  X   Y   Z    R       Amp      RMS    Gamma FWHM \n")


        outparams.write("  %i %i %i %i %.3f %.3f %.3f %.3f %.3f \n" % \
                            (self.filenumber, self.etalon_x, self.etalon_y, self.etalon_z, pars['R'][0],
                             pars['Amplitude'][0],
                             self.rms,
                             pars['Gamma'][0],
                             pars['FWHM'][0]))

        outparams.close()

        return

