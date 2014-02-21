# -*- coding: utf-8 -*-
import os, sys
import numpy as np
import zipfile
from xml.dom import minidom
import xml.parsers.expat


from PyQt4 import QtCore, QtGui


from slitmask import SlitMask
from slitlets import Slitlets

class ReadXMLError(Exception):
    """Class for handling XML reading errors"""
    pass

catcolumn_list=['name', 'targ_ra', 'targ_dec', 'width', 'len1', 'len2', 'tilt', 'mag', 'priority', 'flags']

class CatalogTab:
    def __init__(self, ui, infile=None):
        self.ui = ui
        self.slitlets = Slitlets()
        self.slitmask = Slimask()
        self.infile = infile

    def getxml(self,infile):
        """Return the xml dom file if infile is either an xml file or an rsmt file"""
        try:
            zip = zipfile.ZipFile(infile,'r')
            zip.extract('Slitmask.xml')
            try:
                dom = minidom.parse('Slitmask.xml')
            except xml.parsers.expat.ExpatError, e:
                raise ReadXMLError(e)

        except zipfile.BadZipfile:
            try:
                dom = minidom.parse(infile)
            except xml.parsers.expat.ExpatError,e:
                raise ReadXMLError(e)
        return dom

    def loadcatalog(self):
        """Locate a file and then enter it into the ui"""

        #launch a file IO dialog
        ldir = os.getcwd()
        infile = QtGui.QFileDialog.getOpenFileName(caption="Open Catalog", directory=ldir)

        #load that file into the catalog
        self.entercatalog(str(infile))
        # set the new mask values on the display
        self.ui.lineEditMain_CenRA.setText(str(self.slitmask.center_ra))
        self.ui.lineEditMain_CenDEC.setText(str(self.slitmask.center_dec))
        self.ui.lineEditMain_PA.setText(str(self.slitmask.position_angle))
        self.ui.lineEditMain_Equinox.setText(str(self.slitmask.equinox))
        
        if self.slitmask.target_name:
           self.ui.lineEditMain_TargetName.setText(self.slitmask.target_name)
        if self.slitmask.mask_name:
           self.ui.lineEditMain_MaskName.setText(self.slitmask.mask_name)
        if self.slitmask.creator:
           self.ui.lineEditInfo_Creator.setText(self.slitmask.creator)
        if self.slitmask.proposer:
           self.ui.lineEditInfo_Proposer.setText(self.slitmask.proposer)
        if self.slitmask.proposal_code:
           self.ui.lineEditInfo_ProposalCode.setText(self.slitmask.proposal_code)

    def entercatalog(self, infile=None, form='short'):
        """Given a catalog, enter it into table in the ui"""
   
        #double check that a file exists and if not, then load it
        self.infile = infile
        if self.infile is None: 
           self.loadcatalog()
           return

        # check whether the input file is a rsmt or xml file. if it is a xml file
        # load the xml slitmask info, else load it as an ascii file
        try:
            self.slitmask.readmaskxml(self.getxml(str(infile)))
        except ReadXMLError:
            #enter the file information into the slit_arr
            self.slitlets.readascii(self.infile, form=form)
            self.slitmask.set_MaskPosition()
        print self.slitmask.center_ra, self.slitmask.center_dec
        #check for objects outside the FoV
        self.slitmask.outFoV()
        ### TESTING THE COLLISION CHECKER
        self.slitmask.find_collisions()
        #update the table
        self.updatetabs()



    def updatecatalogtable(self):
        self.slitlets.update_flags()
        rows = self.ui.tableWidgetCat.rowCount()
        for i in range(0,rows):
            self.ui.tableWidgetCat.removeRow(i)

        self.ui.tableWidgetCat.setRowCount(0)
        #enter the information into the table
        for i in range(self.slitlets.nobjects):
            if i > self.ui.tableWidgetCat.rowCount():
                self.ui.tableWidgetCat.insertRow(i-1)
            for j in range(len(catcolumn_list)):
                f = catcolumn_list[j]
                item = self.parseItem(self.slitlets.data[catcolumn_list[j]][i])
                self.ui.tableWidgetCat.setItem(i-1,j,item)

    def updatetabs(self):
       """Task designed for overloading"""
       self.updatecatalogtable()

    def addslitfromcatalog(self):
       """Determine slits elected in the catalog and add them to the catalog"""

       #get the selected items
       sel_list = self.ui.tableWidgetCat.selectedItems()
       print self.ui.tableWidgetCat.selectedRanges()
       
       #for each item in sel_list, 
       #get the item, and determine the parameters from it 
       #and activite the object
       for selitem in sel_list:
           selitem.row() 
           i = selitem.row()
           stext = self.ui.tableWidgetCat.item(i,0).text()
           sid = self.slitlets.findslitlet(str(stext))
           self.slitlets.addtomask(sid)
       self.updatetabs()
           
        
    def parseItem(self, x):
       """Parse an object so it can be entered into the table"""
       if isinstance(x, str):
           return QtGui.QTableWidgetItem(x)
       elif isinstance(x, np.float32):
           return QtGui.QTableWidgetItem('%f' % x)
       elif isinstance(x, float):
           return QtGui.QTableWidgetItem('%f' % x)
       elif isinstance(x, int):
           return QtGui.QTableWidgetItem('%i' % x)
       return QtGui.QTableWidgetItem('')

  



