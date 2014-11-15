import numpy as np
import os, errno
from PyQt4 import QtGui,QtCore
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg
from saltgui import MplCanvas

class FpParallWidget (QtGui.QWidget):

    def __init__(self,parent=None):
        super(FpParallWidget,self).__init__(parent)

        #Load up the data:
        self.loadOutparams()


       #set up the file range panel
        self.rangepanel=QtGui.QWidget()

        # add a label:
        self.FromLabel = QtGui.QLabel("From file number:")
        self.ToLabel = QtGui.QLabel("To file number:")

       #add the name of the file
        self.FromValueLabel = QtGui.QLineEdit(str(min(self.outparams[:,0])))
        self.ToValueLabel = QtGui.QLineEdit(str(max(self.outparams[:,0])))

        # and a button to process the new range
        self.refreshButton = QtGui.QPushButton('Refresh')
        self.refreshButton.clicked.connect(self.plotOutparams)

       #set up file range panel layout

        rangeLayout=QtGui.QGridLayout(self.rangepanel)
        rangeLayout.addWidget(self.FromLabel,0,0,1,1)
        rangeLayout.addWidget(self.FromValueLabel,0,1,1,1)
        rangeLayout.addWidget(self.refreshButton,0,2,2,1)
        rangeLayout.addWidget(self.ToLabel,0,3,1,1)
        rangeLayout.addWidget(self.ToValueLabel,0,4,1,1)


        #add the radio buttons for the choice of x axis...


        self.radioFilenumber= QtGui.QRadioButton("Plot vs Filenumber")
        self.radioX= QtGui.QRadioButton("Plot vs etalon X")
        self.radioY= QtGui.QRadioButton("Plot vs etalon Y")

        #create a gropu for them:
        self.radioGroupX=QtGui.QButtonGroup()
        self.radioGroupX.addButton(self.radioFilenumber)
        self.radioGroupX.addButton(self.radioX)
        self.radioGroupX.addButton(self.radioY)

        #make sure the filenumber is the default
        self.radioFilenumber.setChecked(True)


        #create radio buttons for the choice of y axis:
        self.radioFWHM=QtGui.QRadioButton("Plots vs FWHM")
        self.radioAmp=QtGui.QRadioButton("Plots vs Amplitude")
        
        #add a group for the y axis:
        self.radioGroupY=QtGui.QButtonGroup()
        self.radioGroupY.addButton(self.radioFWHM)
        self.radioGroupY.addButton(self.radioAmp)

        #add a default:
        self.radioFWHM.setChecked(True)


        # display best fit in range:

        self.fitpanel=QtGui.QWidget()

        self.fitLabel = QtGui.QLabel("Lowest FWHM in file range:")
        self.cleanOutparams()
        self.getBestparams()

        fitFileresult="File number: %i" %int(self.bestparams[0])
        fitXresult="X: %i" % int(self.bestparams[1])
        fitYresult="Y: %i" % int(self.bestparams[2])
        fitZresult="Z: %i " % int(self.bestparams[3])
        fitRresult="R: %.1f" % float(self.bestparams[4])
        fitAmpresult="Amplitude: %.1f" % float(self.bestparams[5])
        fitRmsresult="RMS: %.3f" % float(self.bestparams[6])
        fitGammaresult="Gamma: %.2f" % float(self.bestparams[7])
        fitFWHMresult="FWHM: %.3f" % float(self.bestparams[8])


        #add the text to the fit results panel
                                    
        self.fitFile = QtGui.QLabel(fitFileresult)
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
        fitLayout.addWidget(self.fitLabel,0,0,1,4)
        fitLayout.addWidget(self.fitFile,3,0,1,1)
        fitLayout.addWidget(self.fitX,3,1,1,1)
        fitLayout.addWidget(self.fitY,3,2,1,1)
        fitLayout.addWidget(self.fitZ,3,3,1,1)
        fitLayout.addWidget(self.fitR,3,4,1,1)
        fitLayout.addWidget(self.fitAmp,3,5,1,1)
        fitLayout.addWidget(self.fitRms,3,6,1,1)
        fitLayout.addWidget(self.fitGamma,3,7,1,1)
        fitLayout.addWidget(self.fitFWHM,3,8,1,1)
 

        
        #set up the fwhm plot

        self.fwhmplot=MplCanvas()

        self.fwhmaxes=self.fwhmplot.figure.add_subplot(111)

        #connect mouse clicks
        
        self.fwhmplot.mpl_connect('button_press_event',self.onClick)



        #and now we know what the X and Y axis should be, make the fwhm/amp plot

        self.plotOutparams()

        # and check for radio button event signals!

        self.radioGroupX.buttonClicked.connect(self.plotOutparams)
        self.radioGroupY.buttonClicked.connect(self.plotOutparams)

        #Add the X radio buttons to a horizontal layout
        
        self.radiopanel= QtGui.QWidget()

        radioLayout=QtGui.QHBoxLayout(self.radiopanel)
        radioLayout.addWidget(self.radioFilenumber)
        radioLayout.addWidget(self.radioX)
        radioLayout.addWidget(self.radioY)

        #Add the Y radio buttons to a vertical layout

        self.radioYpanel=QtGui.QWidget()
        
        radioYLayout=QtGui.QVBoxLayout(self.radioYpanel)
        radioYLayout.addWidget(self.radioFWHM)
        radioYLayout.addWidget(self.radioAmp)

       # Set up the main layout
        mainLayout = QtGui.QGridLayout()
        mainLayout.addWidget(self.rangepanel,0,0,1,9)
        mainLayout.addWidget(self.fitpanel,1,0,1,9)
        mainLayout.addWidget(self.fwhmplot,2,0,1,4)
        mainLayout.addWidget(self.radioYpanel,2,5,1,1)
        mainLayout.addWidget(self.radiopanel,3,1,1,1)
        self.setLayout(mainLayout)



    def loadOutparams(self):

        self.outparams=np.genfromtxt('outparams', skip_header=1)

        return

    def cleanOutparams(self):

        minFile=float(self.FromValueLabel.text())
        maxFile=float(self.ToValueLabel.text())
#        print "reloading from %i to %i" % (minFile, maxFile)

        self.cleanarr=[]

        mask = (minFile <= self.outparams[:,0]) * (self.outparams[:,0] <= maxFile)
        self.cleanarr = self.outparams[mask]

#        print self.cleanarr[:,0]

        return

    def plotOutparams(self):
        #set up the plot....

        self.cleanOutparams()
        self.fwhmaxes.clear()

        if self.radioFilenumber.isChecked():
            x=self.cleanarr[:,0]
        elif self.radioX.isChecked():
            x=self.cleanarr[:,1]
        elif self.radioY.isChecked():
            x=self.cleanarr[:,2]

        # Work out the Y axis:

        if self.radioFWHM.isChecked():
            y=self.cleanarr[:,8]
        elif self.radioAmp.isChecked():
            y=self.cleanarr[:,5]



        self.fwhmaxes.plot(x, y, 'bo')

#        self.show()
        
        # don't forget to force a redraw!

        self.fwhmplot.draw()

        #ummmm we forgot to update the best fit..
        self.getBestparams()

        self.fitFile.setText("File number: %i" %int(self.bestparams[0]))
        self.fitX.setText("X: %i" % int(self.bestparams[1]))
        self.fitX.setText("X: %i" % int(self.bestparams[1]))
        self.fitY.setText("Y %i:" % int(self.bestparams[2]))
        self.fitZ.setText("Z: %i " % int(self.bestparams[3]))
        self.fitR.setText("R: %.1f" % float(self.bestparams[4]))
        self.fitAmp.setText("Amplitude: %.1f" % float(self.bestparams[5]))
        self.fitRms.setText("RMS: %.2f" % float(self.bestparams[6]))
        self.fitGamma.setText("Gamma: %.2f" % float(self.bestparams[7]))
        self.fitFWHM.setText("FWHM: %.3f" % float(self.bestparams[8]))

#        self.fitpanel.show()

        return

    def onClick(self,event):

#        What's on the X axis?

        if self.radioFilenumber.isChecked():
            mask = (self.cleanarr[:,0]==round(event.xdata))
        elif self.radioX.isChecked():
            mask = (self.cleanarr[:,1]==round(event.xdata))
        elif self.radioY.isChecked():
            mask = (self.cleanarr[:,2]==round(event.xdata))

        # get from the array the first row that matches the X value)
        datapoint = self.cleanarr[mask][0]

        #format it ready for the tooltip:
        text="FileNumber: %i, \nX: %i, \nY: %i, \nZ:%i, \nAmp: %.2f, \nRMS: %.2f, \nGamma: %.2f, \nFWHM: %.3f" % (int(datapoint[0]), int(datapoint[1]),int(datapoint[2]),int(datapoint[3]),datapoint[4],datapoint[6],datapoint[7],datapoint[8])

        #and plonk it on! :)
        QtGui.QToolTip.showText(QtCore.QPoint(338,314),text)

        return

    def getBestparams(self):

        if self.radioFWHM.isChecked():
            self.fitLabel.setText("Lowest FWHM in file range:")
            mask = (self.cleanarr[:,8]==min(self.cleanarr[:,8]))
            self.bestparams = self.cleanarr[mask][0]
        elif self.radioAmp.isChecked():
            self.fitLabel.setText("Highest Amplitude in file range:")
            mask = (self.cleanarr[:,5]==max(self.cleanarr[:,5]))
            self.bestparams = self.cleanarr[mask][0]


        return
                                    
                                    
