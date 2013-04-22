# -*- coding: utf-8 -*-
import os, sys
import numpy as np

from PyQt4 import QtCore, QtGui

from slitlets import Slitlets

catcolumn_list=['name', 'targ_ra', 'targ_dec', 'width', 'len1', 'len2', 'tilt', 'mag', 'priority', 'flags']

class RefTab:
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

    def setrefposition(self):
        '''
        determine the which rows are selected for the deletion on refstars
        '''
        self.rows = []
        indexes = self.ui.tableWidgetRefstars.selectedIndexes()
        for index in indexes:
            self.rows.append(index.row())
        self.rows = self.unique(self.rows)
        print self.rows


    def clearrefstartable(self):
        '''
        set all the in_mask flags of objects with refstar = 1 to 0 and update
        the refstar table
        '''
        # get refstars in the mask
        inmask = np.where((self.slitlets.data['inmask_flag'] == 1)*(self.slitlets.data['refstar_flag'] == 1))
        
        self.slitlets.data['inmask_flag'][inmask] = 0
        self.updaterefstartable()

    def addrefstarmanually(self):
        '''
        add an empty row to the slit table that the user must fill in manually
        '''
#        rows = self.ui.tableWidgetSlits.rowCount()
#        print rows
#        if rows > 0:
#            self.ui.tableWidgetRefstars.setRowCount(rows+1)
##            self.ui.tableWidgetSlits.insertRow(rows+1)
#            for j in range(len(catcolumn_list)):
#                item = self.parseItem('')
#                self.ui.tableWidgetRefstars.setItem(rows+1,j,item)
#        else:
#            self.ui.tableWidgetRefstars.setRowCount(0)
#            self.ui.tableWidgetRefstars.insertRow(0)
#            for j in range(len(catcolumn_list)):
#                item = self.parseItem('')
#                self.ui.tableWidgetRefStars.setItem(rows+1,j,item)

        ######### TESTING: ###########
        self.slitlets.add_slitlet(1)
        self.slitmask.outFoV()
        self.updatetabs()

    def refchanged(self, x ,y):
        """When ever a cell is changed, updated the information about the slit"""
        #identify what slit changed and what attribute of that slit changed
        item = self.ui.tableWidgetRefstars.item(x,0)
        name = str(item.text())
        ai = self.slitlets.findslitlet(name)

        #update the attribute
        item=self.ui.tableWidgetRefstars.item(x,y)
        self.slitlets.updatevalue(ai, str(item.text()), catcolumn_list[y])

        #update all the tabs
        self.slitmask.outFoV_row(ai)
        self.slitmask.find_collisions()
        self.slitlets.update_flags()
        self.updatetabs()


    def deleterefstarmanually(self):
        '''
        set the selected slits inmask_flag to 0, if none were slected, do nothing
        '''

        if len(self.rows) == 0:
            return
        else:
            print self.rows
            for i in self.rows:
                item = self.ui.tableWidgetRefstars.item(i,0)
                name = str(item.text())
                ai = np.where(self.slitlets.data['name'] == name)
                self.slitlets.data[ai[0][0]]['inmask_flag'] = 0
        self.slitlets.update_flags()
        self.updatetabs()
        return

    def addrefstarsfromcatalogue(self):
        '''
        if a refstar has a priority >=0 add it to the refstar table.
        * if slits have been added manually before adding from the catalogue
        the row count is set accordingly
        TODO: add in the FoV checker before loading slits and set the appropriate
        flags
        '''
        refstars = np.where((self.slitlets.data['priority'] == -1) * (self.slitlets.data['fov_flag'] == 1))
        for i in refstars[0]:
            self.slitlets.data[i]['inmask_flag'] = 1
        self.slitlets.update_flags()
        self.updatetabs()
     
    def updaterefstartable(self):
        """Using the slitlet object, update the slit table"""

        # check which entries are in the mask
        inmask = np.where((self.slitlets.data['inmask_flag'] == 1) * (self.slitlets.data['refstar_flag'] == 1))
        
        #enter the information into the table
        self.ui.tableWidgetRefstars.setRowCount(0)
        nobj=0
        for i in inmask[0]:
            self.ui.tableWidgetRefstars.insertRow(nobj)
            for j in range(len(catcolumn_list)):
                item=self.parseItem(self.slitlets.data[catcolumn_list[j]][i])
                self.ui.tableWidgetRefstars.blockSignals(True)
                self.ui.tableWidgetRefstars.setItem(nobj,j,item)
                self.ui.tableWidgetRefstars.blockSignals(False)
            nobj+=1

    def updatetabs(self):
       """Task designed for overloading"""
       self.updaterefstartable()
        
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

