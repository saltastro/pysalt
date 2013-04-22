# -*- coding: utf-8 -*-
import sys
from PyQt4 import QtCore, QtGui
'''

'''

class InfoTab():
    def __init__(self, ui, slitmask, infile=None):
        super(QObject, self).__init__()
        self.ui = ui
        self.infile=infile
    
    # set the mode of operation:
    def setmode2cat(self):
        self.ui.lineEditMain_Mode.setText('Catalogue')
        self.ui.toolButtonCat_Load.setEnabled(True)
        self.mode = 'catalogue'
 
    def setmode2manual(self):
        self.ui.lineEditMain_Mode.setText('Manual')
        self.mode = 'manual'

    def setmodecentroiding(self):
        '''
        set the text on the centroiding display
        '''
        if self.ui.checkBoxInfo_CentroidOn.isChecked():
            self.ui.labelMain_CentroidingOnOff.setText('ON')
        else:
            self.ui.labelMain_CentroidingOnOff.setText('OFF')

    # load the mask info
    def loadcreator(self):
        self.slitmask.validated = False
        self.slitmask.creator=self.ui.lineEditInfo_Creator.text()
        if len(self.ui.lineEditInfo_Creator.text()) == 0:
            self.slitmask.creator = None

    def loadproposer(self):
        self.slitmask.validated = False
        self.slitmask.proposer=self.ui.lineEditInfo_Proposer.text()
        if len(self.ui.lineEditInfo_Proposer.text()) == 0:
            self.slitmask.proposer = None

    def loadproposalcode(self):
        self.slitmask.validated = False
        self.slitmask.proposal_code=self.ui.lineEditInfo_ProposalCode.text()
        if len(self.ui.lineEditInfo_ProposalCode.text()) == 0:
            self.slitmask.proposal_code = None


        


