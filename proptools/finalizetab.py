# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
from PyQt4 import QtCore, QtGui
import zipfile

try:
    from PyQt4.QtCore import QString
except ImportError:
    QString = str

from finder_chart import finderchart

'''

'''


class ValidatorError(Exception):
    """Class for handling errors during the validation process"""

    def __init__(self, ui, errmesg):
        self.ui = ui
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(225, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Text, brush)
        self.ui.textEditFin_ErrorMessages.setPalette(palette)
        self.ui.textEditFin_ErrorMessages.setText(errmesg)
    pass

class FinalizeTab:
    def __init__(self, ui, infile=None):
        self.ui = ui
        self.infile = infile
        self.xmlfile = None
        self.fcfile = None

    def updatetabs(self):
       """Task designed for overloading"""
       pass
        

    def validation_success(self, msg=''):
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(0, 170, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Text, brush)
        self.ui.textEditFin_ErrorMessages.setPalette(palette)
        valmesg = msg + '\n ** Validation was successful! **\n'
        self.ui.textEditFin_ErrorMessages.setText(valmesg)


    def validator(self):
        '''
        go through the slitmask class and make sure that all the info neccessary
        for a successful mask has been provided and is correct
        * produce errors for information necessary for the gcode
        * produce warnings for info that is not needed but might have been
          missed by the user
        '''
        
        errmsg = ''
        warnmsg = ''
        err_encountered = False
        errmesg_hdr = '''\nThe following ERRORS were encountered: \n\n'''

        warn_encountered = False
        warnmesg_hdr = '''\nThe following WARNINGS were encountered: \n\n'''

        errmsg = errmsg + errmesg_hdr
        warnmsg = warnmsg + warnmesg_hdr

        # go through all the fields and check for errors:
        if self.slitmask.center_ra == None:
            errmsg = ' - Center RA not defined\n'
            err_encountered = True

        if self.slitmask.center_dec == None:
            errmsg = errmsg + ' - Center DEC not defined\n'
            err_encountered = True

        if self.slitmask.position_angle == None:
            errmsg = errmsg + ' - Position angle not defined\n'
            err_encountered = True

        if len(self.slitmask.target_name) == 0:
            errmsg = errmsg + ' - Target name not defined\n'
            err_encountered = True

        if len(self.slitmask.mask_name) == 0:
            errmsg = errmsg + ' - Mask name not defined\n'
            err_encountered = True

        if self.slitmask.equinox == None:
            errmsg = errmsg + ' - Equinox not defined\n'
            err_encountered = True

        if len(self.slitmask.proposal_code) == 0:
            errmsg = errmsg + ' - Proposal code is not defined\n'
            err_encountered = True

        if len(self.slitmask.proposer) == 0:
            errmsg = errmsg + ' - Proposer not defined\n'
            err_encountered = True

        if len(self.slitmask.creator) == 0:
            errmsg = errmsg + ' - Creator not defined\n'
            err_encountered = True

        if len(np.where((self.slitlets.data['refstar_flag'] == 1) & \
                    (self.slitlets.data['inmask_flag'] == 1))[0]) < 1:
            errmsg = errmsg + ' - Less than 3 refstar defined\n'
            err_encountered = True


        # warn the user about info that might be missing
        if len(np.where((self.slitlets.data['refstar_flag'] != 1) & \
                    (self.slitlets.data['inmask_flag'] == 1))[0]) == 0:
            warnmsg = warnmsg + ' - No slits were defined\n'
            warn_encountered = True

        self.slitmask.find_collisions()
        self.updatetabs()

        if len(np.where((self.slitlets.data['collision_flag'] == 1)*(self.slitlets.data['refstar_flag'] != 1) & \
                    (self.slitlets.data['inmask_flag'] == 1))[0]) > 0:
            warnmsg = warnmsg + ' - Unresolved slit collisions\n'
            warn_encountered = True

        if err_encountered:
            self.slitmask.validated = False
            if warn_encountered:
                msg = warnmsg + errmsg
                raise ValidatorError(self.ui, msg)
            else:
                raise ValidatorError(self.ui, errmsg)
        else:
            self.slitmask.validated = True
            if warn_encountered:
                self.validation_success(warnmsg)
            else:
                self.validation_success()




    def writersmt(self):
        '''
        opens a dialog box to obtain the name of the file to be saved
        * TODO:
        * check the validation values and produce and error if it has not been
        validated
        * add the writing of the xml file as Slitmask.xml
        * zip the rsmt file with the current jpg of the slits
        '''
        if not self.slitmask.validated:
            errmsg = '''The slitmask has not been validated, please validate it before proceeding.
            '''
            raise ValidatorError(self.ui,errmsg)

        ldir = os.getcwd()
        rsmtoutfilename = QtGui.QFileDialog.getSaveFileName(caption="Save RSMT File", directory=ldir)
        rsmtoutfilename = str(rsmtoutfilename).strip('.rsmt') + '.rsmt'

        #write the xml
        self.writexml2file('Slitmask.xml')
        # get the finder chart image to package with the rsmt zipfile.
        try:
             print 'writing with image'
             self.writeFC(image=self.inimage, outfile='Slitmask.png')
        except:
             print 'writing without image'
             self.writeFC(outfile='Slitmask.png')
        rsmtfile = zipfile.ZipFile(rsmtoutfilename, mode='w')
        try:
            rsmtfile.write('Slitmask.xml')#, rsmtoutfilename, zipfile.ZIP_DEFLATED)
            rsmtfile.write('Slitmask.png')
        finally:
            rsmtfile.close()

        print 'Finished writing rsmt file'

    def writexml(self):
        '''
        writes the current mask info to an user specified xml file
        '''
        ldir=os.getcwd()
        xmloutfilename = QtGui.QFileDialog.getSaveFileName(caption="Save Current XML file", directory=ldir)
        self.writexml2file(str(xmloutfilename))



    def writexml2file(self,outfile):
        '''
        * get the xml from the slitmask class
        * write the xml to the specified output file
        '''
        # make sure that the output filename contains the .xml extension
        outfile = outfile.strip('.xml') + '.xml'
        self.xmlfile=outfile
        if os.path.isfile(outfile): os.remove(outfile)
        xml_txt = self.slitmask.writexml()
        rsmt_file = open(outfile,'w')
        rsmt_file.write(xml_txt)
        rsmt_file.close()
 
    def writeFC_DSS(self):
        self.writeFC()

    def writeFC_Current(self):
        self.writeFC(image=self.inimage)

    def writeFC(self, image=None, outfile=None):
        """Write a finder chart for the current mask"""

        #check to make sure it is valid
        if not self.slitmask.validated:
            msg="Please validate your slitmask before creating a finding chart"
            raise ValidatorError(self.ui, msg)


        #determine the name of the output file
        if not outfile:
            ldir = os.getcwd()
            outfile = str(QtGui.QFileDialog.getSaveFileName(caption="Save Finding Chart file", directory=ldir))
            self.fcfile = outfile
        else:
            self.fcfile=outfile

        if os.path.isfile(outfile): os.remove(outfile)

        #if an xml file exists, then write it out
        if self.xmlfile:
            xmlfile = self.xmlfile
            tmpfile = None
        else:
            tmpfile = 'tmprsmttmp.xml'
            self.writexml2file(tmpfile)
            xmlfile = tmpfile

        #create the finding chart
        print 'running finderchart', xmlfile, image, self.fcfile
        finderchart(xmlfile, image=image, outfile=self.fcfile)
        
#        remove any temporary files
        #if tmpfile: os.remove(tmpfile)
