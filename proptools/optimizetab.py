# -*- coding: utf-8 -*-
import os, sys
import numpy as np

from PyQt4 import QtCore, QtGui

from slitlets import Slitlets
import pyslit_optimize as opt

class OptimizeTab:
    def __init__(self, ui, default_yspacing=1., default_iter=10):
        print 'loading OPT'
        self.ui = ui
        self.slitlets=Slitlets()
        self.opt_yspacing = default_yspacing
        self.opt_niter = default_iter
        
    def setoptimizer_yspacing(self):
        self.opt_yspacing = self.checkyspacing_input(self.ui.lineEditOpt_Yspacing.text())
        
    def setoptimizer_iter(self):
        self.opt_niter = self.checkniter_input(self.ui.lineEditOpt_Niter.text())

    def includerefstars(self):
        if self.ui.checkBoxOpt_IncRefstars.isChecked():
            nrefstars = len(np.where(self.slitlets.data['priority'] == -1)[0])
            self.ui.lineEditOpt_AllRefstars.setText(str(nrefstars))
        else:
            self.ui.lineEditOpt_AllRefstars.setText('')

    def setnumrefstars(self):
        print self.ui.lineEditOpt_NumRefstars.text()

        
    def optimize(self): 
        """Run the optimizer program and optimize the slits"""
        msg = "Optimize the Slitlets"
        print msg
        cra = self.slitmask.center_ra
        cdec = self.slitmask.center_dec
        rotang = self.slitmask.position_angle
        equinox = 2000
        is_in_fov = np.where(self.slitlets.data['fov_flag'] == 1)[0]

        # jpk: this will need to be added in the next version
#        is_in_fov = np.where((self.slitlets.data['inmask_flag'] == 1) * (self.slitlets.data['fov_flag'] == 1))[0]

        ra = self.slitlets.data['targ_ra']
        dec = self.slitlets.data['targ_dec']
        pri = self.slitlets.data['priority']
        slen = self.slitlets.data['len1'] + self.slitlets.data['len2'] 
        swid = self.slitlets.data['width']
        stilt = self.slitlets.data['tilt']
        
        Nstars_req = 0. # **** Paul: I'm not quite sure where to get this from ****
        Niter=10 # **** as above ****
        # set all inmask flags to zero before running optimiser
        #self.slitlets.emptymask()
        
        # -- only run this on objects within FOV:        
        ok = is_in_fov
        if not ok.any(): 
           print "No objects in the field of view--update mask center and run again"
           return

        print ra[ok]
        tra = ra[ok]
        tdec = dec[ok]
        tpri = pri[ok]
        tslen = slen[ok]
        tswid = swid[ok]
        tstilt = stilt[ok]
        print 'running optimizer'
 
        tin_mask = opt.pyslit_optimize(cra, cdec, rotang, equinox, tra, tdec, \
                                    tpri,tslen,tswid,tstilt,\
                                    Niter,self.opt_yspacing, Nstars_req)
        # apply index numbers to full list:
        in_mask = ok[tin_mask]

        # reset all the in_mask values, otherwise the objects which should not
        # be in the optimized mask will still have a in_mask flag 
        self.slitlets.data['inmask_flag'] = 0
        self.slitlets.data['collision_flag'] = 0

        # now add the in_mask flag to the sources which was found by the
        # optimizer 
        for sid in in_mask:
            self.slitlets.addtomask(sid)
        self.updatetabs()


    def updatetabs(self):
        self.slitmask.outFoV_all()
        self.slitmask.find_collisions()
        self.slitlets.update_flags()
#        pass
    
    def checkyspacing_input(self, x):
        try:
            val = float(x)
            if val > 0:
                return val
            else:
                self.opt_yspacing = 1
                self.ui.lineEditOpt_Yspacing.setText(str(self.opt_yspacing))
        except ValueError,e:
            self.opt_yspacing = 1
            self.ui.lineEditOpt_Yspacing.setText(str(self.opt_yspacing))

    def checkniter_input(self,x):
        try:
            val = int(x)
            if val > 0:
                return val
            else:
                self.opt_niter = 10
                self.ui.lineEditOpt_Niter.setText(str(self.opt_niter))
        except ValueError,e:
            self.opt_niter = 10
            self.ui.lineEditOpt_Niter.setText(str(self.opt_niter))
