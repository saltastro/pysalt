"""
DQWidget is a Qt4 Widget for displaying information about the data quality 
of an image
"""
import numpy as np 
import pyfits

from PyQt4 import QtGui, QtCore
from ObsLogWidget import headerList, printList

from saltstat import iterstat


#import plugins
from rssinfo import rssinfo
from seeing import seeing_stats, airmass

class DQWidget(QtGui.QWidget):
   def __init__(self, name, imlist, parent=None):
       super(DQWidget, self).__init__(parent)
       self.imlist=imlist
       self.name=name
       #set up the information panel
       self.infopanel=QtGui.QWidget()
       self.infopanel.setFixedHeight(50)
 
       #set up some needed variables
       self.obsmode=self.getitem('OBSMODE')

       #add the name of the file
       self.NameLabel = QtGui.QLabel("Filename:")
       self.NameLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised )
       self.NameValueLabel = QtGui.QLabel("%s" % self.name)
       self.NameValueLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken )

       #set up the info panel layout
       infoLayout=QtGui.QGridLayout(self.infopanel)
       infoLayout.addWidget(self.NameLabel, 0, 0, 1, 1)
       infoLayout.addWidget(self.NameValueLabel, 0, 1, 1, 1)

       #set up the panel for the different modes
       if self.obsmode=='IMAGING':
          self.datapanel=self.set_imaging()
       elif self.obsmode=='SPECTROSCOPY': 
          self.datapanel=self.set_spectroscopy()
       else:
          self.datapanel=QtGui.QWidget()

       # Set up the layout
       mainLayout = QtGui.QVBoxLayout()
       mainLayout.addWidget(self.infopanel)
       mainLayout.addWidget(self.datapanel)
       self.setLayout(mainLayout)

   def updatetab(self, name, imlist):
       print "STRARTING UPDATE OF DQ"
       self.imlist=imlist
       self.name=name
       self.obsmode=self.getitem('OBSMODE')

       #add the name of the file
       self.NameValueLabel.setText(("%s" % self.name))
       self.NameValueLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken )

       #set up the panel for the different modes
       if self.obsmode=='IMAGING':
          self.datapanel=self.set_imaging()
       elif self.obsmode=='SPECTROSCOPY': 
          self.datapanel=self.set_spectroscopy()
       else:
          self.datapanel=QtGui.QWidget()

       # Set up the layout

       j=self.layout().indexOf(self.datapanel)
       print "INDEX", j
       self.layout().itemAt(1).widget().close()
       self.layout().insertWidget(1, self.datapanel)

   def set_datardx(self):
       """Set up the information from the data reduction

          number of cr clean
          bias levels removed
       """
            

   def set_spectroscopy(self):

       #set up the data panel
       datapanel=QtGui.QWidget()

       #get the infomration that you need about the image
       grating=self.getitem('GRATING').strip()
       slitname=self.getitem('MASKID').strip()
       graang=float(self.getitem('GR-ANGLE'))
       artang=float(self.getitem('AR-ANGLE'))
       xbin, ybin=self.getitem('CCDSUM').split()

       #get the information about the model
       wcen, w1, w2, res, R, slitsize=rssinfo(grating, graang, artang, slitname, xbin, ybin)
       print grating, slitname, graang, artang, xbin, ybin

       #Information to include in the data panel
       #central w, w1, w2, resolution, dw
       self.gratingLabel = QtGui.QLabel("Grating")
       self.gratingLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
       self.graangLabel = QtGui.QLabel("GR-ANGLE")
       self.graangLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
       self.artangLabel = QtGui.QLabel("AR-ANGLE")
       self.artangLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
       self.slitnameLabel = QtGui.QLabel("SLIT")
       self.slitnameLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
       self.slitsizeLabel = QtGui.QLabel("SIZE")
       self.slitsizeLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
       self.gratingValueLabel = QtGui.QLabel(grating)
       self.gratingValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)
       self.graangValueLabel = QtGui.QLabel(u"%5.3f \u00B0" % graang)
       self.graangValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)
       self.artangValueLabel = QtGui.QLabel(u"%5.3f \u00B0" % artang)
       self.artangValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)
       self.slitnameValueLabel = QtGui.QLabel(slitname)
       self.slitnameValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)
       self.slitsizeValueLabel = QtGui.QLabel("%3.2f''" % slitsize)
       self.slitsizeValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)

       self.bluewaveLabel = QtGui.QLabel("Blue Edge")
       self.bluewaveLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
       self.centwaveLabel = QtGui.QLabel("Center Wave")
       self.centwaveLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
       self.redwaveLabel = QtGui.QLabel("Red Edge")
       self.redwaveLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
       self.resolutionLabel = QtGui.QLabel("Resolution")
       self.resolutionLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
       self.reselementLabel = QtGui.QLabel("Resolution element")
       self.reselementLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
       self.bluewaveValueLabel = QtGui.QLabel(u"%7.2f \u00c5" % w1)
       self.bluewaveValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)
       self.centwaveValueLabel = QtGui.QLabel(u"%7.2f \u00c5" % wcen)
       self.centwaveValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)
       self.redwaveValueLabel = QtGui.QLabel(u"%7.2f \u00c5" % w2)
       self.redwaveValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)
       self.resolutionValueLabel = QtGui.QLabel("%5i" % R)
       self.resolutionValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)
       self.reselementValueLabel = QtGui.QLabel(u"%4.2f \u00c5" % res)
       self.reselementValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)

       #set up the signal to noise
       
       specfile='smbxp'+self.name.replace('fits', 'txt')
       try:
           w,f,s=np.loadtxt(specfile, usecols=(0,1,2), unpack=True)
           med_sn=np.median(s)
           try:
               s.sort()
               peak_sn=np.median(s[-100:])
           except:
               peak_sn=s.max()
       except:
          med_sn=0
          peak_sn=0

       self.snLabel = QtGui.QLabel("Median S/N")
       self.snLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Raised)
       self.snValueLabel = QtGui.QLabel("%5.2f" % med_sn)
       self.snValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)
       self.psnLabel = QtGui.QLabel("Median Peak S/N")
       self.psnLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Raised)
       self.psnValueLabel = QtGui.QLabel("%5.2f" % peak_sn)
       self.psnValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)

       #set the layout
       dataLayout=QtGui.QGridLayout(datapanel)
       dataLayout.addWidget(self.gratingLabel, 0, 0, 1, 1)
       dataLayout.addWidget(self.graangLabel, 0, 1, 1, 1)
       dataLayout.addWidget(self.artangLabel, 0, 2, 1, 1)
       dataLayout.addWidget(self.slitnameLabel, 0, 3, 1, 1)
       dataLayout.addWidget(self.slitsizeLabel, 0, 4, 1, 1)
       dataLayout.addWidget(self.gratingValueLabel, 1, 0, 1, 1)
       dataLayout.addWidget(self.graangValueLabel, 1, 1, 1, 1)
       dataLayout.addWidget(self.artangValueLabel, 1, 2, 1, 1)
       dataLayout.addWidget(self.slitnameValueLabel, 1, 3, 1, 1)
       dataLayout.addWidget(self.slitsizeValueLabel, 1, 4, 1, 1)
       dataLayout.addWidget(self.bluewaveLabel, 2, 0, 1, 1)
       dataLayout.addWidget(self.centwaveLabel, 2, 1, 1, 1)
       dataLayout.addWidget(self.redwaveLabel, 2, 2, 1, 1)
       dataLayout.addWidget(self.resolutionLabel, 2, 3, 1, 1)
       dataLayout.addWidget(self.reselementLabel, 2, 4, 1, 1)
       dataLayout.addWidget(self.bluewaveValueLabel, 3, 0, 1, 1)
       dataLayout.addWidget(self.centwaveValueLabel, 3, 1, 1, 1)
       dataLayout.addWidget(self.redwaveValueLabel, 3, 2, 1, 1)
       dataLayout.addWidget(self.resolutionValueLabel, 3, 3, 1, 1)
       dataLayout.addWidget(self.reselementValueLabel, 3, 4, 1, 1)
   
       dataLayout.addWidget(self.snLabel, 4, 0, 2, 3)
       dataLayout.addWidget(self.snValueLabel, 4, 3, 2, 2)
       dataLayout.addWidget(self.psnLabel, 6, 0, 2, 3)
       dataLayout.addWidget(self.psnValueLabel, 6, 3, 2, 2)

       datapanel.setFixedHeight(400)
       return datapanel

   def set_imaging(self):
       #set up the data panel
       datapanel=QtGui.QWidget()

       #get some variables from the data
       exptime=float(self.getitem('EXPTIME'))
       filtername=self.getitem('FILTER')
       telalt=float(self.getitem('TELALT'))
       ccdbin=float(self.getitem('CCDSUM').split()[0])
       print ccdbin
       pix_scale=0.14*ccdbin
       z=90-telalt 
       am=airmass(z)

       #determine the background
       try:
          bmean=float(self.getitem('BMEAN'))
          bmidpt=float(self.getitem('BMIDPT'))
          bstd=float(self.getitem('BSTD'))
       except Exception, e:
          outimg='mbxp'+self.name
          #hdu=pyfits.open(outimg)
          #bmean, bmidpt, bstd=iterstat(hdu[1].data, 5, 3)
          #hdu.close()
          bmean, bmidpt, bstd=(-1,-1,-1)

       #calculate the seeing
       try: 
          see=float(self.getitem('SEEING'))
          nsource=float(self.getitem('NSOURCES'))
       except Exception, e:
          outtxt='mbxp'+self.name.replace('fits', 'cat')
          try:
             mag_arr, fwhm_arr=np.loadtxt(outtxt, usecols=(2,10), unpack=True)
          except IOError:
             mag_arr=np.zeros([0])
             fwhm_arr=np.zeros([0])
          nsource=0 #len(mag_arr)
          mean, std, norm, peak=seeing_stats(fwhm_arr)
          see=mean*pix_scale
          print e
          seeing=-1

       #Display the filter, pix_scale, airmass, exptime
       self.filterLabel = QtGui.QLabel("FILTER")
       self.filterLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
       self.exptimeLabel = QtGui.QLabel("EXPTIME")
       self.exptimeLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
       self.pixscaleLabel = QtGui.QLabel("PIXSCALE")
       self.pixscaleLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
       self.airmassLabel = QtGui.QLabel("AIRMASS")
       self.airmassLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)

       self.filterValueLabel = QtGui.QLabel(filtername)
       self.filterValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)
       self.exptimeValueLabel = QtGui.QLabel('%6.2f s' % exptime)
       self.exptimeValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)
       self.pixscaleValueLabel = QtGui.QLabel("%3.2f ''/pix" % pix_scale)
       self.pixscaleValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)
       self.airmassValueLabel = QtGui.QLabel('%3.2f' % am)
       self.airmassValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)

       #display  the background stats
       self.bmeanLabel = QtGui.QLabel("BACKGROUND MEAN")
       self.bmeanLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
       self.bmidptLabel = QtGui.QLabel("BACKGROUND MIDPT")
       self.bmidptLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
       self.bstdLabel = QtGui.QLabel("BACKGROUND STD")
       self.bstdLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
       self.bmeanValueLabel = QtGui.QLabel('%5.2f' % bmean)
       self.bmeanValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)
       self.bmidptValueLabel = QtGui.QLabel('%5.2f' % bmidpt)
       self.bmidptValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)
       self.bstdValueLabel = QtGui.QLabel('%5.2f' % bstd)
       self.bstdValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)
       
       #display seeing, number of sources
       self.seeLabel = QtGui.QLabel("SEEING")
       self.seeLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
       self.seeValueLabel = QtGui.QLabel("%3.2f'' " % see)
       self.seeValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)
       self.sourceLabel = QtGui.QLabel("SOURCES Detected")
       self.sourceLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
       self.sourceValueLabel = QtGui.QLabel("%i" % nsource)
       self.sourceValueLabel.setFrameStyle(QtGui.QFrame.Panel |  QtGui.QFrame.Sunken)

       #add it to the panel
       dataLayout=QtGui.QGridLayout(datapanel)
       dataLayout.addWidget(self.filterLabel, 0, 0, 1, 1)
       dataLayout.addWidget(self.exptimeLabel, 0, 1, 1, 1)
       dataLayout.addWidget(self.pixscaleLabel, 0, 2, 1, 1)
       dataLayout.addWidget(self.airmassLabel, 0, 3, 1, 1)

       dataLayout.addWidget(self.filterValueLabel, 1, 0, 1, 1)
       dataLayout.addWidget(self.exptimeValueLabel, 1, 1, 1, 1)
       dataLayout.addWidget(self.pixscaleValueLabel, 1, 2, 1, 1)
       dataLayout.addWidget(self.airmassValueLabel, 1, 3, 1, 1)

       dataLayout.addWidget(self.bmeanLabel, 2, 0, 1, 1)
       dataLayout.addWidget(self.bmidptLabel, 2, 1, 1, 1)
       dataLayout.addWidget(self.bstdLabel, 2, 2, 1, 1)

       dataLayout.addWidget(self.bmeanValueLabel, 3, 0, 1, 1)
       dataLayout.addWidget(self.bmidptValueLabel, 3, 1, 1, 1)
       dataLayout.addWidget(self.bstdValueLabel, 3, 2, 1, 1)

       dataLayout.addWidget(self.seeLabel, 4, 0, 1, 1)
       dataLayout.addWidget(self.seeValueLabel, 4, 1, 1, 1)
       dataLayout.addWidget(self.sourceLabel, 5, 0, 1, 1)
       dataLayout.addWidget(self.sourceValueLabel, 5, 1, 1, 1)

       datapanel.setFixedHeight(400)

       return datapanel

   def getitem(self, key):
       i=headerList.index(key)
       try:
           value=str(self.imlist[i])
       except IndexError:
           value=''
       return value

