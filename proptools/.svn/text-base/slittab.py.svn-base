# -*- coding: utf-8 -*-
import os, sys
import numpy as np

from PyQt4 import QtCore, QtGui

from slitlets import Slitlets

catcolumn_list=['name', 'targ_ra', 'targ_dec', 'width', 'len1', 'len2', 'tilt', 'mag', 'priority', 'flags']

class SlitTab:
    def __init__(self, ui, infile=None):
        self.ui = ui
        self.slitlets=Slitlets()
        self.infile=infile
        self.cell = None
        self.rows = None

    def unique(self,seq):
        '''
        return a unique number of rows from the table selection
        '''
        seen = set()
        seen_add = seen.add
        return [ x for x in seq if x not in seen and not seen_add(x)]

    def setposition(self):
        '''
        determine the which rows are selected for the deletion on slits
        '''
        self.rows = []
        indexes = self.ui.tableWidgetSlits.selectedIndexes()
        for index in indexes:
            self.rows.append(index.row())
        self.rows = self.unique(self.rows)
        print self.rows


    def clearslittable(self):
        '''
        set all the in_mask flags to 0 and update the slit table
        '''
        rows = self.ui.tableWidgetSlits.rowCount()
        print rows
        self.slitlets.data['inmask_flag'] = 0
        self.updatetabs()
#        self.ui.tableWidgetSlits.clear()
#        for i in range(0,rows):
#            print i
#            self.ui.tableWidgetSlits.removeRow(i)

    def addslitmanually(self):
        '''
        add an empy row to the slit table that the user must fill in manually
        '''
#        rows = self.ui.tableWidgetSlits.rowCount()
#        print rows
#        if rows > 0:
#            self.ui.tableWidgetSlits.setRowCount(rows+1)
##            self.ui.tableWidgetSlits.insertRow(rows+1)
#            for j in range(len(catcolumn_list)):
#                item = self.parseItem('')
#                self.ui.tableWidgetSlits.setItem(rows+1,j,item)
#        else:
#            self.ui.tableWidgetSlits.setRowCount(0)
#            self.ui.tableWidgetSlits.insertRow(0)
#            for j in range(len(catcolumn_list)):
#                item = self.parseItem('')
#                self.ui.tableWidgetSlits.setItem(rows+1,j,item)

        ######### TESTING: ###########
        self.slitlets.add_slitlet()
        self.slitmask.outFoV()
        self.updatetabs()

    def deleteslitmanually(self):
        '''
        set the selected slits inmask_flag to 0, if none were slected, do nothing
        '''

        if len(self.rows) == 0:
            return
        else:
            print self.rows
            for i in self.rows:
                item = self.ui.tableWidgetSlits.item(i,0)
                name = str(item.text())
                ai = np.where(self.slitlets.data['name'] == name)
                self.slitlets.data[ai[0][0]]['inmask_flag'] = 0
        self.slitlets.update_flags()
        self.updatetabs()
        return

    def addslitletsfromcatalogue(self):
        '''
        if a slit has a priority >=0 add it to the slits table.
        * if slits have been added manually before adding from the catalogue
        the row count is set accordingly
        TODO: add in the FoV checker before loading slits and set the appropriate
        flags
        '''
        slits = np.where((self.slitlets.data['priority'] >= 0) * (self.slitlets.data['fov_flag'] == 1))
        for i in slits[0]:
            self.slitlets.data[i]['inmask_flag'] = 1
        self.slitmask.outFoV_all()
        self.slitmask.find_collisions()
        self.slitlets.update_flags()
        self.updatetabs()
        
#        rows = self.ui.tableWidgetSlits.rowCount()
#        if rows > 0:
#            self.ui.tableWidgetSlits.setRowCount(rows+1)
#            row = rows+1
#        else:
#            self.ui.tableWidgetSlits.setRowCount(0)
#            row = 0
#
#        for i in slits[0]:
#            self.ui.tableWidgetSlits.insertRow(row)
#            for j in range(len(catcolumn_list)-1):
#                item=self.parseItem(self.slitlets.data[catcolumn_list[j]][i])
#                self.ui.tableWidgetSlits.setItem(row,j,item)
#            row+=1
       


    def updateslittable(self):
        """Using the slitlet object, update the slit table"""

        # check which entries are in the mask and not reference stars
        inmask = np.where((self.slitlets.data['inmask_flag'] == 1)*(self.slitlets.data['refstar_flag'] == 0))
        
        #enter the information into the table
        self.ui.tableWidgetSlits.setRowCount(0)
        nobj = 0
        for i in inmask[0]:
            self.ui.tableWidgetSlits.insertRow(nobj)
            for j in range(len(catcolumn_list)):
                item=self.parseItem(self.slitlets.data[catcolumn_list[j]][i])
                self.ui.tableWidgetSlits.blockSignals(True)
                self.ui.tableWidgetSlits.setItem(nobj,j,item)
                self.ui.tableWidgetSlits.blockSignals(False)
            nobj += 1
#        for i in range(self.slitlets.nobjects):
#            if self.slitlets.data['inmask_flag'][i]==1:
#
#                if nobj>self.ui.tableWidgetSlits.rowCount():
#                    self.ui.tableWidgetSlits.insertRow(nobj)
#                for j in range(len(catcolumn_list)-1):
#                    f=catcolumn_list[j]
#                    item=self.parseItem(self.slitlets.data[catcolumn_list[j]][i])
#                    self.ui.tableWidgetSlits.setItem(nobj,j,item)
#            nobj+=1

    def slitchanged(self, x ,y):
        """When ever a cell is changed, updated the information about the slit"""
        #identify what slit changed and what attribute of that slit changed
        item = self.ui.tableWidgetSlits.item(x,0)
        name = str(item.text())
        ai = self.slitlets.findslitlet(name)

        #update the attribute
        item=self.ui.tableWidgetSlits.item(x,y)
        self.slitlets.updatevalue(ai, str(item.text()), catcolumn_list[y])

        #update all the tabs
        self.slitmask.outFoV_row(ai)
        self.slitmask.find_collisions()
        self.slitlets.update_flags()
        self.updatetabs()
        
   


    def updatetabs(self):
       """Task designed for overloading"""
       self.updateslittable()
       self.updatecatalogtable()
        
    def parseItem(self, x):
       """Parse an object so it can be entered into the table"""
       if isinstance(x, str):
           return QtGui.QTableWidgetItem(x)
       elif isinstance(x, numpy.float32):
           return QtGui.QTableWidgetItem('%f' % x)
       elif isinstance(x, float):
           return QtGui.QTableWidgetItem('%f' % x)
       elif isinstance(x, int):
           return QtGui.QTableWidgetItem('%i' % x)
       return QtGui.QTableWidgetItem('')

