# -*- coding: utf-8 -*-
import os,sys
import numpy as np
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import SIGNAL, SLOT, QObject

from slitlets import Slitlets, ra_read, dec_read
from slitmask import SlitMask

from rsmt_gui import Ui_MainWindow
from infotab import InfoTab
from catalogtab import CatalogTab
from slittab import SlitTab
from reftab import RefTab
from optimizetab import OptimizeTab
from finalizetab import FinalizeTab

# added these two import to avoid a seg fault
from pyraf import iraf
from iraf import pysalt

from ImageDisplay import ImageDisplay

class SlitMaskGui(QtGui.QMainWindow, InfoTab, CatalogTab, OptimizeTab, SlitTab, RefTab, FinalizeTab):
    def __init__(self, parent=None, infile=None, inimage=None, center_ra=None, center_dec=None, position_angle=None):
        QtGui.QWidget.__init__(self, parent)
        
        #set up the main UI
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
#        OptimizeTab.__init__(self,self.ui)
        #set up the slitmask
        self.slitmask=SlitMask(center_ra=center_ra, center_dec=center_dec, position_angle=position_angle )
        self.slitlets=self.slitmask.slitlets

        #setup the image interaction
        self.imagedisplay=ImageDisplay(target='pySlitMask:5909')
        self.position_angle=position_angle 
        if inimage:
           self.loadimage(inimage)
        #set up some variables that will be needed later
        self.xmlfile=None
        self.fcfile=None

        #read in the input data if available
        self.infile=infile 
        if infile:
            self.ui.radioButtonInfo_Catalogue.setChecked(True)
            self.setmode2cat()
            self.entercatalog(infile)
            
            if self.slitmask.center_ra is None and self.slitmask.center_dec is None:
                self.slitmask.set_MaskPosition()
                self.displayfootprint()

            #if self.slitmask.center_ra and self.slitmask.center_dec:
                #self.imagedisplay.rssregion(self.slitmask.center_ra, self.slitmask.center_dec)
            
            self.ui.lineEditMain_CenRA.setText(str(self.slitmask.center_ra))
            self.ui.lineEditMain_CenDEC.setText(str(self.slitmask.center_dec))
            self.ui.lineEditMain_PA.setText(str(self.slitmask.position_angle))
            self.ui.lineEditMain_Equinox.setText(str(self.slitmask.equinox))
            
            self.ui.lineEditMain_TargetName.setText(self.slitmask.target_name)
            self.ui.lineEditMain_MaskName.setText(self.slitmask.mask_name)
            self.ui.lineEditInfo_Creator.setText(self.slitmask.creator)
            self.ui.lineEditInfo_Proposer.setText(self.slitmask.proposer)
            self.ui.lineEditInfo_ProposalCode.setText(self.slitmask.proposal_code)
        else:
            self.ui.radioButtonInfo_Catalogue.setChecked(False)
            self.ui.radioButtonInfo_Manual.setChecked(True)
            self.setmode2manual()
            self.ui.toolButtonCat_Load.setEnabled(True)



              #self.displayfootprint()

#        self.opttab = OptimizeTab()
        # setup default values for the optimizer
        self.opt_yspacing = 1.
        self.opt_iter = 10
        self.ui.lineEditOpt_Yspacing.setText(str(self.opt_yspacing))
        self.ui.lineEditOpt_Niter.setText(str(self.opt_iter))

#        self.slitmask.outFoV()
#        print self.slitlets.data['fov_flag']
#        self.updatetabs()
        #Listen to different signals

        #menu items
        QtCore.QObject.connect(self.ui.actionLoad_Catalogue, QtCore.SIGNAL("triggered()"), self.loadcatalog)
        QtCore.QObject.connect(self.ui.actionLoad_Image, QtCore.SIGNAL("triggered()"), self.loadimage)

        #main tabs
        QtCore.QObject.connect(self.ui.lineEditMain_CenRA, QtCore.SIGNAL("editingFinished()"), self.loadCenRA)
        QtCore.QObject.connect(self.ui.lineEditMain_CenDEC, QtCore.SIGNAL("editingFinished()"), self.loadCenDEC)
        QtCore.QObject.connect(self.ui.lineEditMain_PA, QtCore.SIGNAL("editingFinished()"), self.loadpositionangle)
        QtCore.QObject.connect(self.ui.lineEditMain_Equinox, QtCore.SIGNAL("editingFinished()"), self.loadequinox)
        QtCore.QObject.connect(self.ui.lineEditMain_TargetName, QtCore.SIGNAL("editingFinished()"), self.loadtargetname)
        QtCore.QObject.connect(self.ui.lineEditMain_MaskName, QtCore.SIGNAL("editingFinished()"), self.loadmaskname)
        
        #info tabs
        QtCore.QObject.connect(self.ui.lineEditInfo_ProposalCode, QtCore.SIGNAL("editingFinished()"), self.loadproposalcode)
        QtCore.QObject.connect(self.ui.lineEditInfo_Proposer, QtCore.SIGNAL("editingFinished()"), self.loadproposer)
        QtCore.QObject.connect(self.ui.lineEditInfo_Creator, QtCore.SIGNAL("editingFinished()"), self.loadcreator)

#        QtCore.QObject.connect(self.slitmask, SIGNAL('xmlloaded'), self.setcreator)

        QtCore.QObject.connect(self.ui.radioButtonInfo_Catalogue, QtCore.SIGNAL("clicked()"), self.setmode2cat)
        QtCore.QObject.connect(self.ui.radioButtonInfo_Manual, QtCore.SIGNAL("clicked()"), self.setmode2manual)
        QtCore.QObject.connect(self.ui.checkBoxInfo_CentroidOn, QtCore.SIGNAL("clicked()"), self.setmodecentroiding)

        
        
        #catalog tabs
        QtCore.QObject.connect(self.ui.toolButtonCat_Load, QtCore.SIGNAL("clicked(bool)"), self.loadcatalog)
        QtCore.QObject.connect(self.ui.pushButtonCat_AddSlits, QtCore.SIGNAL("clicked(bool)"), self.addslitfromcatalog)
        QtCore.QObject.connect(self.ui.pushButtonCat_Clear, QtCore.SIGNAL("clicked()"), self.clearContents)

        #slit tab
        QtCore.QObject.connect(self.ui.pushButtonSlit_ClearSlits, QtCore.SIGNAL("clicked()"), self.clearslittable)
        QtCore.QObject.connect(self.ui.pushButtonSlit_AddSlitImage, QtCore.SIGNAL("clicked()"), self.addslitletsfromimage)
        QtCore.QObject.connect(self.ui.pushButtonSlit_AddSlitfromCat, QtCore.SIGNAL("clicked()"), self.addslitletsfromcatalogue)
        QtCore.QObject.connect(self.ui.pushButtonSlit_AddSlit, QtCore.SIGNAL("clicked()"), self.addslitmanually)
        QtCore.QObject.connect(self.ui.pushButtonSlit_DeleteSlit, QtCore.SIGNAL("clicked()"), self.deleteslitmanually)
        QtCore.QObject.connect(self.ui.pushButtonSlit_DeleteSlitImage, QtCore.SIGNAL("clicked()"), self.deleteslitfromimage)
        QtCore.QObject.connect(self.ui.tableWidgetSlits, QtCore.SIGNAL("itemSelectionChanged()"), self.setposition)
        QtCore.QObject.connect(self.ui.tableWidgetSlits, QtCore.SIGNAL("cellChanged(int, int)"), self.slitchanged)

        #optimize tab
        QtCore.QObject.connect(self.ui.pushButtonOpt_Optimize, QtCore.SIGNAL("clicked()"), self.optimize)
        QtCore.QObject.connect(self.ui.lineEditOpt_Yspacing, QtCore.SIGNAL("editingFinished()"), self.setoptimizer_yspacing)
        QtCore.QObject.connect(self.ui.lineEditOpt_Niter, QtCore.SIGNAL("editingFinished()"), self.setoptimizer_iter)
        QtCore.QObject.connect(self.ui.checkBoxOpt_IncRefstars, QtCore.SIGNAL("stateChanged(int)"), self.includerefstars)
        QtCore.QObject.connect(self.ui.lineEditOpt_NumRefstars, QtCore.SIGNAL("editingFinished()"), self.setnumrefstars)

        #ref stars
        QtCore.QObject.connect(self.ui.pushButtonRef_ClearRefstars, QtCore.SIGNAL("clicked()"), self.clearrefstartable)
        QtCore.QObject.connect(self.ui.pushButtonRef_AddRefstarImage, QtCore.SIGNAL("clicked()"), self.addslitletsfromimage)
        QtCore.QObject.connect(self.ui.pushButtonRef_AddRefstarsfromCat, QtCore.SIGNAL("clicked()"), self.addrefstarsfromcatalogue)
        QtCore.QObject.connect(self.ui.pushButtonRef_AddRefstar, QtCore.SIGNAL("clicked()"), self.addrefstarmanually)
        QtCore.QObject.connect(self.ui.pushButtonRef_DeleteRefstar, QtCore.SIGNAL("clicked()"), self.deleterefstarmanually)
        QtCore.QObject.connect(self.ui.pushButtonRef_DeleteRefstar_2, QtCore.SIGNAL("clicked()"), self.deleteslitfromimage)
        QtCore.QObject.connect(self.ui.tableWidgetRefstars, QtCore.SIGNAL(" itemSelectionChanged()"), self.setrefposition)
        QtCore.QObject.connect(self.ui.tableWidgetRefstars, QtCore.SIGNAL("cellChanged(int, int)"), self.refchanged)

        # finalize tab
        QtCore.QObject.connect(self.ui.pushButtonFin_Validate, QtCore.SIGNAL("clicked(bool)"),self.validator)
        QtCore.QObject.connect(self.ui.pushButtonFin_WriteXML, QtCore.SIGNAL("clicked(bool)"), self.writexml)
        QtCore.QObject.connect(self.ui.toolButtonFin_WriteRSMT, QtCore.SIGNAL("clicked(bool)"), self.writersmt)
        QtCore.QObject.connect(self.ui.pushButtonFin_CreateFChart_Current, QtCore.SIGNAL("clicked(bool)"), self.writeFC_Current)
        QtCore.QObject.connect(self.ui.pushButtonFin_CreateFChart_DSS, QtCore.SIGNAL("clicked(bool)"), self.writeFC_DSS)

 
        self.ui.tabWidget.setCurrentIndex(0)

    def clearContents(self):
        self.slitlets.data = None
        self.ui.tableWidgetCat.clearContents()
        self.ui.tableWidgetCat.setRowCount(0)
        #TODO: Set the number of rows to the current data length
        #print 'nope not doing it'

    # loads mask coordinates
    def loadCenRA(self):
        self.slitmask.validated = False
        self.slitmask.add_center_ra(ra_read(self.ui.lineEditMain_CenRA.text()))
        if self.slitmask.center_ra == None:
            palette = self.setPalette('error')
            self.ui.lineEditMain_CenRA.setPalette(palette)
        else:
            palette = self.setPalette('normal')
            self.ui.lineEditMain_CenRA.setPalette(palette)
            self.ui.lineEditMain_CenRA.setText(str(self.slitmask.center_ra))
            self.slitmask.outFoV()
            self.updatetabs()

    def loadCenDEC(self):
        self.slitmask.validated = False
        self.slitmask.add_center_dec(dec_read(self.ui.lineEditMain_CenDEC.text()))
        if self.slitmask.center_dec == None:
            palette = self.setPalette('error')
            self.ui.lineEditMain_CenDEC.setPalette(palette)
        else:
            palette = self.setPalette('normal')
            self.ui.lineEditMain_CenDEC.setPalette(palette)
            self.ui.lineEditMain_CenDEC.setText(str(self.slitmask.center_dec))
            self.slitmask.outFoV()
            self.updatetabs()

    def loadpositionangle(self):
        #print self.ui.lineEditMain_PA.text()
        self.slitmask.validated = False
        self.slitmask.add_position_angle(dec_read(self.ui.lineEditMain_PA.text()))
        if self.slitmask.position_angle == None:
            palette = self.setPalette('error')
            self.ui.lineEditMain_PA.setPalette(palette)
        else:
            palette = self.setPalette('normal')
            self.ui.lineEditMain_PA.setPalette(palette)
            self.ui.lineEditMain_PA.setText(str(self.slitmask.position_angle))
            self.slitmask.outFoV()
            self.imagedisplay.rotate(self.slitmask.position_angle)
            self.updatetabs()

    def loadequinox(self):
        self.slitmask.validated = False
        self.slitmask.add_equinox(dec_read(self.ui.lineEditMain_Equinox.text()))
        if self.slitmask.equinox == None:
            palette = self.setPalette('error')
            self.ui.lineEditMain_Equinox.setPalette(palette)
        else:
            palette = self.setPalette('normal')
            self.ui.lineEditMain_Equinox.setPalette(palette)
            self.ui.lineEditMain_Equinox.setText(str(self.slitmask.equinox))
            self.slitmask.outFoV()
            
     # load info from the main window
    def loadtargetname(self):
        self.slitmask.validated = False
        self.slitmask.target_name=self.ui.lineEditMain_TargetName.text()
        if len(self.ui.lineEditMain_TargetName.text()) == 0:
            self.slitmask.target_name = None

    def loadmaskname(self):
        self.slitmask.validated = False
        self.slitmask.mask_name=self.ui.lineEditMain_MaskName.text()
        if len(self.ui.lineEditMain_MaskName.text()) == 0:
            self.slitmask.mask_name = None
   
    def loadValue(self):
        self.updatetabs()

    def updatetabs(self):
       """Update all of the information after changes to the slitlet class"""
       print 'Updating tabs'
       self.updatecatalogtable()
       self.updateslittable()
       self.updaterefstartable()
       self.imagedisplay.deleteregions()
       self.displayslits()
#       self.displayfootprint()

    def displayslits(self): 
       """Add the slits to the image """
       ids = np.where(self.slitlets.data['inmask_flag']==1)[0]
       fout = open('tmp.reg', 'w')
       fout.write('# Region file format: DS9 version 4.1\n# Filename: sgpR.fits\n')

       fout.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
       fout.write('fk5\n')

       for i in ids:
           fout.write(self.slitlets.asregion(i,self.slitmask.position_angle)+'\n')
           # show spectral footprint:
           #fout.write(self.slitlets.asregionspec(i,self.slitmask.position_angle)+'\n')
       fout.close()

       self.imagedisplay.regionfromfile('tmp.reg')


    def displayall(self): 
       """Add the slits to the image """
       ids=np.where(self.slitlets.data['inmask_flag']==1)[0]
       fout=open('oth.reg', 'w')
       fout.write('# Region file format: DS9 version 4.1\n# Filename: sgpR.fits\n')

       fout.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
       fout.write('fk5\n')

       for i in ids:
           #fout.write(self.slitlets.asregion(i,self.slitmask.position_angle)+'\n')
           # show spectral footprint:
           fout.write(self.slitlets.asregionspec(i,self.slitmask.position_angle)+'\n')
           # **** also stars, need to be careful about shape
       fout.close()

       self.imagedisplay.regionfromfile('tmp.reg')


    def displayfootprint(self):
        """Add the RSS footprint to the image"""
#       fout=open('tmpfp.reg', 'w')
#       fout.write('# Region file format: DS9 version 4.1\n# Filename: sgpR.fits\n')
#
##       fout.write('global color=yellow dashlist=8 3 width=2 font="helvetica 10 normal roman" select=0 highlite=0 dash=1 fixed=0 edit=0 move=0 delete=0 include=1 source=1\n')
#       fout.write('fk5\n')
#       fout.write('circle(%f,%f,4\') # color=yellow dashlist=8 3 width=2 select=0 highlite=0 dash=1 fixed=0 edit=0 move=1 delete=0 background include=1 source=0\n' % (self.slitmask.center_ra, self.slitmask.center_dec))
#       fout.close()
#       self.imagedisplay.regionfromfile('tmpfp.reg')
        self.imagedisplay.rssregion(self.slitmask.center_ra, self.slitmask.center_dec)

    def loadimage(self, inimage=None):
        if not inimage:
             #launch a file IO dialog
             ldir = os.getcwd()
             inimage = QtGui.QFileDialog.getOpenFileName(caption="Open Catalog", directory=ldir)

        self.inimage=str(inimage)
        self.imagedisplay.display(inimage, pa=self.position_angle)

    def deleteslitfromimage(self):
        #download the slits in the image
        newslits=self.imagedisplay.getregions()
        

        #loop through the list and see which one is missing
        try:
            index=np.where(self.slitlets.data['inmask_flag']==1)[0]
        except:
            return
        
 
        #check to see if it is in the mask
        for i in index:
            sra=self.slitlets.data['targ_ra'][i]
            sdec=self.slitlets.data['targ_dec'][i]
            found=False
            for k in newslits.keys():
              ra = float(newslits[k][0][0])
              dec = float(newslits[k][0][1])
              if abs(sra-ra) < 0.0003 and abs(sdec-dec) < 0.0003:
                 found=True
            if not found:
               self.slitlets.data['inmask_flag'][i]=0

        #update the tabs
        self.updatetabs()      


    def addslitletsfromimage(self):
        """Download the slits from the image and add them to the slitlet or catalog

          If catalog is selected, it will search the catalog for a corresponding object and center the slit on that object

          **TODO**
          If manual centroided is selected, it will centroid around that value in the image and use that position
          If manual uncentroided is selected, it will just use the slit position

        """
        #download the slits in the image
        newslits=self.imagedisplay.getregions()
        #loop through the objects--if they are already in the catalog check to see if they
        #need updating.  If not, then add them to the catalog

        #print "Sorting through regions"

        for i in newslits:
            if newslits[i][1]:
              #if it is tagged we assume it is already in the slitmask
              #print newslits[i]
              pass
            else:
              ra = float(newslits[i][0][0])
              dec = float(newslits[i][0][1])
              #print i,ra,dec
           #width=str(newslits[i][0])
           #height=str(newslits[i][0])
           #tilt=str(newslits[i][0])
              #TODO: This searches that catalog and adds the target that matches the slit drawn
              sid=self.slitlets.findtarget(ra,dec)
              self.slitlets.addtomask(sid)

         
        self.updatetabs()

    def setPalette(self,mode):
        palette = QtGui.QPalette()

        if mode == 'error':
            brush = QtGui.QBrush(QtGui.QColor(255, 148, 148))
            brush.setStyle(QtCore.Qt.SolidPattern)
            palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Base, brush)
            palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Base, brush)

            return palette

        if mode == 'normal':
            brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
            brush.setStyle(QtCore.Qt.SolidPattern)
            palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Base, brush)
            palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Base, brush)

            return palette

if __name__ == "__main__":
  infile = None
  inimage= None
  if len(sys.argv)>1:
     infile=sys.argv[1]
  if len(sys.argv)==3:
     inimage=sys.argv[2]
  app = QtGui.QApplication([])
  myapp = SlitMaskGui(infile=infile, inimage=inimage)
  myapp.show()
  sys.exit(app.exec_())

