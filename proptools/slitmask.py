# -*- coding: utf-8 -*

import numpy as np
import pywcs
from slitlets import Slitlets
from xml.dom import minidom, getDOMImplementation
import xml.parsers.expat

from PyQt4.QtCore import Qt, QVariant, QObject, SIGNAL, QAbstractTableModel, QModelIndex

class SlitMaskError(Exception):
    """Base class for exceptions from the slitmaks class"""
    pass


class SlitMask(QObject):
    def __init__(self, center_ra=None, center_dec=None, position_angle=0, target_name='', mask_name='',
                 equinox=2000, proposal_code='', proposer='', creator='', validated=False):
        super(QObject, self).__init__();
        self.__center_ra = None
        self.__center_dec = None
        self.__position_angle = None
        self.__equinox = None
        self.add_center_ra(center_ra)
        self.add_center_dec(center_dec)
        self.add_position_angle(position_angle)
        self.add_equinox(equinox)
        
        self.target_name = target_name
        self.mask_name = mask_name
        self.proposal_code = proposal_code
        self.proposer = proposer
        self.creator = creator
        self.validated = validated

        #create the slitlets
        self.slitlets = Slitlets()

        #needed for writing xml
        self.impl = getDOMImplementation()

#        self.connect(self, SIGNAL('xmlloaded'), self.printcreator)
        
    @property
    def center_ra(self):
        return self.__center_ra
    @property
    def center_dec(self):
        return self.__center_dec
    @property
    def position_angle(self):
        return self.__position_angle
    @property
    def equinox(self):
        return self.__equinox

    def add_center_ra(self, new_center_ra):
        if new_center_ra is None: 
            self.__center_ra=None
            return
        if new_center_ra < 0 or new_center_ra >= 360:
            self.__center_ra=None
#            raise SlitMaskError()
            return
        changed = new_center_ra != self.__center_ra
        self.__center_ra = new_center_ra
        if changed:
            pass
            # send signal

    def add_center_dec(self, new_center_dec):

        if new_center_dec is None : 
            self.__center_dec=None
            return
        if new_center_dec < -90 or new_center_dec > 90:
            self.__center_dec=None
#            raise SlitMaskError()
            return
        changed = new_center_dec != self.__center_dec
        self.__center_dec = new_center_dec
        if changed:
            pass
            # send signal

    def add_position_angle(self, new_position_angle):
        """Update the value of the position angle"""
        #if None, keep as none
        if new_position_angle is None:
            self.__position_angle=None
            return
        #if a str is given, try to update to a float

        if new_position_angle < 0 or new_position_angle > 360:
            self.__position_angle=None
            return
        changed = new_position_angle != self.__position_angle
        self.__position_angle = new_position_angle
        if changed:
            pass
            # send signal

    def add_equinox(self, new_equinox):
            """Update the equinox value"""
            #if None, keep as none
            if new_equinox is None:
                self.__equinox=None
                return
            if new_equinox < 0:
                self.__equinox=None
                return
#                raise SlitMaskError()
            changed = new_equinox != self.__equinox
            self.__equinox = new_equinox
            if changed:
                pass
                # send signal

    def set_MaskPosition(self):
       """Sets the central mask position from the current slit positions.  The setting of the mask position works on the following
          priority: 1) If an object is givne with a priority of 2, that object sets the masks position.  2) If objects are
          pre-selected to be in the mask, they set the center position of the mask.  3) Set the center position to be at the 
          center of the catalog
       """
       #get the arrays that you might need
       priority = self.slitlets.data['priority']
       inmask = self.slitlets.data['inmask_flag']
       ra = self.slitlets.data['targ_ra']
       dec = self.slitlets.data['targ_dec']

       if priority.any() >= 2:
          cen_ra = ra[priority == 2].mean()
          cen_dec = dec[priority == 2].mean()
       elif inmask.any():
          cen_ra = ra[inmask == 1].mean()
          cen_dec = dec[inmask == 1].mean()
       else:
          cen_ra = ra.mean()
          cen_dec = dec.mean()
 
       self.add_center_ra(cen_ra)
       self.add_center_dec(cen_dec)
       self.add_position_angle(0)

    #
    # ... analogous methods for other properties ...
    #

    def loadxmlinfo(self,par):

        self.proposal_code = par['proposalcode']
        self.proposer = par['pi']
        self.creator = par['creator']
        self.mask_name = par['masknum']
        self.validated = par['validated']
        self.add_center_ra(float(par['centerra']))
        self.add_center_dec(float(par['centerdec']))
        self.add_position_angle(float(par['rotangle']))
        try:
           self.target_name = par['target']
        except:
           pass

#        c = str(self.creator)
#        self.emit(SIGNAL('xmlloaded'), self.creator)
      
        
#for testing the emit and connect
#    def printcreator(self, c):
#        print 'the creator is: ', c
        
    def readmaskxml(self, dom):
        # read all the parameters into dictionaries
        parameters = dom.getElementsByTagName('parameter')

        Param = {}
        for param in parameters:
            Param[str(param.getAttribute('name')).lower()] \
            = str(param.getAttribute('value')).lower()


        # read all the reference stars into dictionaries
        Refstars = {}
        refstars = dom.getElementsByTagName('refstar')
        for refstar in refstars:
            t = {}
            t['xce'] =  float(refstar.getAttribute('xce'))
            t['yce'] =  float(refstar.getAttribute('yce'))
            t['mag'] = float(refstar.getAttribute('mag'))
            t['id'] = str(refstar.getAttribute('id'))
            Refstars[refstar.getAttribute('id')] = t

        # read all the slits into dictionaries
        Slits = {}
        slits = dom.getElementsByTagName('slit')
        for slit in slits:
            t = {}
            t['xce'] = float(slit.getAttribute('xce'))
            t['yce'] = float(slit.getAttribute('yce'))
            t['width'] = float(slit.getAttribute('width'))
            t['length'] = float(slit.getAttribute('length'))
            t['len1'] = float(slit.getAttribute('length')) / 2.
            t['len2'] = float(slit.getAttribute('length')) / 2.
            t['priority'] = float(slit.getAttribute('priority'))
            t['mag'] = float(slit.getAttribute('mag'))
            t['id'] = str(slit.getAttribute('id'))
            Slits[slit.getAttribute('id')] = t
         
        self.loadxmlinfo(Param)
        self.slitlets.readxml(Slits, Refstars)

        # need to add an update all the param field when loading an xml

    def addxmlparameter(self,name,value):
        parameter = self.doc.createElement("parameter")
        parameter.setAttribute('name','%s'%name)
        parameter.setAttribute('value','%s'%value)
        return parameter

    def writexml(self):
        '''
        write out the slitmask and slits info to a xml file.
        '''

        print 'writing xml...'
         #create the xml documents and the main Element called slitmask
        self.doc = self.impl.createDocument(None, "slitmask", None)
        slitmask= self.doc.documentElement
        header = self.doc.createElement("header")
        slitmask.appendChild(header)

        header.appendChild(self.addxmlparameter("VERSION","1.1"))
        header.appendChild(self.addxmlparameter("PROPOSALCODE","%s"%self.proposal_code))
        header.appendChild(self.addxmlparameter("MASKNUM","%s"%self.mask_name))
        header.appendChild(self.addxmlparameter("TARGET","%s"%self.target_name))
        header.appendChild(self.addxmlparameter("PI","%s"%self.proposer))
        header.appendChild(self.addxmlparameter("CREATOR","%s"%self.creator))
        header.appendChild(self.addxmlparameter("ROTANGLE","%s"%self.position_angle))
        header.appendChild(self.addxmlparameter("CENTERRA","%f"%self.center_ra))
        header.appendChild(self.addxmlparameter("CENTERDEC","%f"%self.center_dec))
        header.appendChild(self.addxmlparameter("NSMODE","0"))
        header.appendChild(self.addxmlparameter("VALIDATED","%s"%str(self.validated)))
        header.appendChild(self.addxmlparameter("SPECLENGTH","12400"))
        header.appendChild(self.addxmlparameter("SPECOFFSET","0"))
        header.appendChild(self.addxmlparameter("SPECPOLSPLIT","0"))
        header.appendChild(self.addxmlparameter("SPECHEIGHT","0"))


        for i in range(0,len(self.slitlets.data)):
          if self.slitlets.data['inmask_flag'][i]:
            slitcard = self.slitlets.asxml(i)
            slitmask.appendChild(slitcard)

        xml = self.doc.toprettyxml(indent="  ")
        return xml

    def outFoV_row(self,i):
        '''
        checks if a single row in the slitlets.data array is outside the FoV
        input is a slitlets.data row
        '''

        # first check that the mask info has been defined, otherwise all
        # slits fall outside the FoV
        if self.center_ra == None or self.center_dec == None\
            or self.position_angle == None or self.equinox == None:
                self.slitlets.data[i]['fov_flag'] = 0
                return
        else:
            pass

        pixscale = 0.2507 / 2. # unbinned pixels

        dcr = 4. / 60. # radius of field (deg)
        # global CCD parameters:
        ccd_dx = 2034.
        ccd_xgap = 70.
        ccd_dy = 4102.

        # define centre in pixel coords
        ccd_cx = (2. * (ccd_dx + ccd_xgap) + ccd_dx) / 2.
        ccd_cy = ccd_dy / 2.

        wcs = pywcs.WCS(naxis=2)
        wcs.wcs.crpix = [ccd_cx,ccd_cy]
        wcs.wcs.cdelt = np.array([-pixscale, pixscale]) / 3600. # set in degrees
        wcs.wcs.crval = [self.center_ra, self.center_dec]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        wcs.wcs.crota = [self.position_angle, self.position_angle] # rotate SKY this amount?
        wcs.wcs.equinox = self.equinox

        # convert onsky coords to pix
        xp,yp = wcs.wcs_sky2pix(self.slitlets.data[i]['targ_ra'],\
                                self.slitlets.data[i]['targ_dec'], 1)

        # testing for all four slit corners
        # upper left corner
        ulx = xp - (self.slitlets.data[i]['slit_width'] /2. ) / pixscale
        uly = yp + self.slitlets.data[i]['len1'] / pixscale
        # upper right corner
        urx = xp + (self.slitlets.data[i]['slit_width'] / 2.)/ pixscale
        ury = yp + self.slitlets.data[i]['len1'] / pixscale
        # lower left corner
        llx = xp - (self.slitlets.data[i]['slit_width'] / 2.) / pixscale
        lly = yp - self.slitlets.data[i]['len2'] / pixscale
        # lower right corner
        lrx = xp + (self.slitlets.data[i]['slit_width'] /2. ) / pixscale
        lry = yp - self.slitlets.data[i]['len2'] / pixscale

        # determine the distances to each corner
        uldist = (uly - ccd_cy)**2 + (ulx - ccd_cx)**2
        urdist = (ury - ccd_cy)**2 + (urx - ccd_cx)**2
        lldist = (lly - ccd_cy)**2 + (llx - ccd_cx)**2
        lrdist = (lry - ccd_cy)**2 + (lrx - ccd_cx)**2

        # test if each corner lies outside the FoV, the boolean array returns
        # true if a corner lies outside the FoV
        maxdist = (dcr * 3600.0 / pixscale)**2
        ultest = uldist <= maxdist
        urtest = urdist <= maxdist
        lltest = lldist <= maxdist
        lrtest = lrdist <= maxdist

        # the slitlet lies in the FoV only if all tests were passed
        # if only one test is failed its False
    
        pass_test = ultest * urtest * lltest * lrtest

        # convert the list to a np array
        pass_test = np.array(pass_test)


        # set the FoV flag for slits inside the FoV
        self.slitlets.data[i]['fov_flag'] = 1 * pass_test[0]

        return

    def outFoV(self):
        '''
        this function goes throught slitlets.data and populates the FoV flag
        using the outFoV_row function for each entry in the array
        '''
        
        for i in range(0,len(self.slitlets.data)):
            self.outFoV_row(i)
            
        return

    def find_collisions(self):
        '''
        this function checks for slit collisions
        '''
        idebug = 1
        if idebug: print "Checking for collisions"

        # xleft,xright are spectral lengths in pixels - set to large dummy values for now

        # **** again, only does special case of PA=0/180, no slit tilt...
        #if not ( self.position_angle == 0 or self.position_angle == 180):
        #   msg='Collision checking only works at position angles of 0 or 180'
        #   print msg
        #   self.slitlets.data['collision_flag'] = 0
        #   return

        # first check that the mask info has been defined, otherwise
        # no collision checking can be done
        if self.center_ra == None or self.center_dec == None\
            or self.position_angle == None or self.equinox == None:
                self.slitlets.data['collision_flag'] = 0
                return
        else:
            pass

        print self.center_ra, self.center_dec, self.position_angle

        # setup default values:
        #slit_ra,slit_dec,slit_length,slit_width,cra,cdec,rotang,equinox,xpad,ypad,xleft,xright

        pixscale = 0.129 # unbinned pixels
#        xpad = 1
#        xright = 1
#        xleft = 1
        ypad = 1. / pixscale



        # global CCD parameters:
        ccd_dx = 2034.
        ccd_xgap = 70.
        ccd_dy = 4102.

        # define centre in pixel coords
        ccd_cx = (2. * (ccd_dx + ccd_xgap) + ccd_dx) / 2.
        ccd_cy = ccd_dy / 2.

        wcs = pywcs.WCS(naxis=2)
        wcs.wcs.crpix = [ccd_cx,ccd_cy]
        wcs.wcs.cdelt = np.array([-pixscale, pixscale]) / 3600. # set in degrees
        wcs.wcs.crval = [self.center_ra, self.center_dec]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        wcs.wcs.crota = [self.position_angle, self.position_angle] # rotate SKY this amount?
        wcs.wcs.equinox = self.equinox

        #TODO: only slits that lie in the FoV should be checked for collisions
        is_in_mask = np.where((self.slitlets.data['inmask_flag'] == 1) * (self.slitlets.data['fov_flag'] == 1))[0]
        print 'is in mask', is_in_mask
        if not is_in_mask.any():  
           msg='No objects in mask'
           print msg
           self.slitlets.data['collision_flag'] = 0
           return

        # convert onsky coords to pix
        xp, yp = wcs.wcs_sky2pix(self.slitlets.data['targ_ra'][is_in_mask],
                                 self.slitlets.data['targ_dec'][is_in_mask], 1)
       
        slit_length = (self.slitlets.data['len1'][is_in_mask] + self.slitlets.data['len2'][is_in_mask]) / pixscale
        nslits = len(is_in_mask)
        coll_flag = [False] * nslits
        coll_ids = [None] * nslits ## can assign a list of collisions to any item in this list

        # **** i think for general case i just need to replace slit width below with projected height which would be:
        #      total_height = slit_length*cos(theta) + 2.* 0.5 * slit_width*sin(theta)
        #         where theta is some combination of rotang and tilt
        #
        # need to update optimiser for more general case, too.
#        tx0 = xp - xleft - (xpad / 2.)
#        tx1 = xp + xright + (xpad / 2.) # [ careful switching between lower left and width and lower left and upper right notation ]
        ymin = yp - (slit_length / 2.) - (ypad / 2.)
        ymax = yp + (slit_length / 2.) + (ypad / 2.)

        # dumb but reliable way of checking:
        for i in range(nslits):
            self.slitlets.data['collision_flag'][is_in_mask[i]] = 0
            self.slitlets.data['collision_id'][is_in_mask[i]]= []

            tcoll_ids=[]
            for j in range(nslits):
                if i != j:
                   if (ymax[i] > ymin[j] and ymax[i] < ymax[j]) or (ymin[i] < ymax[j] and ymin[i] > ymin[j]):
                        self.slitlets.data['collision_flag'][is_in_mask[i]] = 1
                        tcoll_ids.append(self.slitlets.data[is_in_mask[j]]['name'])
                        print self.slitlets.data['name'][is_in_mask[i]], self.slitlets.data['collision_flag'][is_in_mask[i]]
                        
            if self.slitlets.data['collision_flag'][is_in_mask[i]]:
               self.slitlets.data['collision_id'][is_in_mask[i]] = " ".join(['%s' % x for x in tcoll_ids])
               print self.slitlets.data['collision_id'][is_in_mask[i]]

        print 'Finished Checking for collisions'
        self.slitlets.update_flags()
        
    def update_fov_slitlets(self):
        changed = False
        if changed:
            # send signal
            pass

    def add_slitlet(self, slitlet, in_mask=False):
        self.__all_slitlets[slitlet.id] = slitlet
        # connect to slit
        # send signal
        if in_mask:
            self.__mask_slitlets.add(slitlet.id)
            # send signal

    def remove_slitlet(self, slitlet):
        del(self.__all_slitlets[slitlet.id])
        # disconnect from slitlet
        # send signal
        if slitlet.id in self.__mask_slitlets:
            self.__mask_slitlets.remove(slitlet.id)
            # send signal
        if slitlet.id in self.__fov_slitlets:
            self.__fov_slitlets.remove(slitlet.id)
            # send signal

    def add_to_mask(self, slitlet):
        if not slitlet.id in self.__all_slitlets.keys():
            raise SlitError()
        if slitlet.id in self.__mask_slitlets:
            return
        self.__mask_slitlets.add(slitlet.id)
        # send signal
        # update FOV list?
        # update collisions

    def remove_from_mask(self, slitlet):
        if not slitlet.id in self.__mask_slitlets:
            return
        self.__mask_slitlets.remove(slitlet.id)
        # send signal
        # update FOV list?
        # update collisions

    def slitlet_changed(self, slitlet):
        """callback method for handling slit changes
        """
        # update FOV list?
        # update collisions
        pass

    @staticmethod
    def is_in_fov(slitlet):
        pass

    @staticmethod
    def read_from_file():
        pass

#    def FOVTest(cra,cdec,equinox,rotang,slitra,slitdec,slit_length,tilt):
    def outFoV_all(self):

        # first check that the mask info has been defined, otherwise all
        # slits fall outside the FoV
        if self.center_ra == None or self.center_dec == None\
            or self.position_angle == None or self.equinox == None:
                self.slitlets.data['fov_flag'] = 0
                return
        else:
            pass

        pixscale = 0.2507 / 2. # unbinned pixels

        dcr = 4. / 60. # radius of field (deg)
        # global CCD parameters:
        ccd_dx = 2034.
        ccd_xgap = 70.
        ccd_dy = 4102.

        # define centre in pixel coords
        ccd_cx = (2.*(ccd_dx + ccd_xgap) + ccd_dx) / 2.
        ccd_cy = ccd_dy / 2.

        # setup the field WCS coords.
        wcs = pywcs.WCS(naxis=2)
        wcs.wcs.crpix = [ccd_cx,ccd_cy]
        wcs.wcs.cdelt = np.array([-pixscale, pixscale]) / 3600. # set in degrees
        wcs.wcs.crval = [self.center_ra, self.center_dec]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        wcs.wcs.crota = [self.position_angle, self.position_angle] # rotate SKY this amount?
        wcs.wcs.equinox = self.equinox

        # convert onsky coords to pix
        xp, yp = wcs.wcs_sky2pix(self.slitlets.data['targ_ra'],
                                self.slitlets.data['targ_dec'], 1)

        # testing for all four slit corners
        # upper left corner
        ulx = xp - (self.slitlets.data['slit_width'] /2. ) / pixscale
        uly = yp + self.slitlets.data['len1'] / pixscale
        # upper right corner
        urx = xp + (self.slitlets.data['slit_width'] / 2.) / pixscale
        ury = yp + self.slitlets.data['len1'] / pixscale
        # lower left corner
        llx = xp - (self.slitlets.data['slit_width'] / 2.) / pixscale
        lly = yp - self.slitlets.data['len2'] / pixscale
        # lower right corner
        lrx = xp + (self.slitlets.data['slit_width'] /2. ) / pixscale
        lry = yp - self.slitlets.data['len2'] / pixscale

        # determine the distances to each corner
        uldist = (uly - ccd_cy)**2 + (ulx - ccd_cx)**2
        urdist = (ury - ccd_cy)**2 + (urx - ccd_cx)**2
        lldist = (lly - ccd_cy)**2 + (llx - ccd_cx)**2
        lrdist = (lry - ccd_cy)**2 + (lrx - ccd_cx)**2

        # test if each corner lies outside the FoV, the boolean array returns
        # true if a corner lies outside the FoV
        maxdist = (dcr * 3600.0 / pixscale)**2
        ultest = uldist <= maxdist
        urtest = urdist <= maxdist
        lltest = lldist <= maxdist
        lrtest = lrdist <= maxdist

        # the slitlet lies in the FoV only if all tests were passed
        # if only one test is failed its False
        pass_test = ultest * urtest * lltest * lrtest

        # convert the list to a np array
        pass_test = np.array(pass_test)

        # set the FoV flag for slits inside the FoV

        self.slitlets.data['fov_flag'] = 1 * pass_test

        return


if __name__=='__main__':
    import sys
   
