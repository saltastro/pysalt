
import numpy as np
import pywcs

# **** check that rotang we are using agrees with telescope definition! ****


# -- set geometry for RSS (and write region file)

# should probably have some smarter way of storing these global parameters

pxscale=0.2507/2. # unbinned

dcr=4./60. # radius of field (deg)
#dcr=3.9/60. # radius of field (deg)
# shrink radius slightly to avoid part of slit being drawn outside of mask circle?


# global CCD parameters:
ccd_dx=2034.
ccd_xgap=70.
ccd_dy=4102.

# define centre in pixel coords
ccd_cx=(2.*(ccd_dx+ccd_xgap)+ccd_dx)/2.
ccd_cy=ccd_dy/2.


def RSSskyfromPix(cra,cdec,rotang,equinox,ccd1,ccd2,ccd3,circ_fov):
    # -- Make a WCS to convert to RSS pixel coords
    # Create a new WCS object.  The number of axes must be set
    # from the start
    wcs = pywcs.WCS(naxis=2)
    wcs.wcs.crpix = [ccd_cx,ccd_cy] 
#    wcs.wcs.crpix = [0,0] # define centre relative to zero

    #wcs.wcs.cdelt = np.array([-pxscale, pxscale])
    wcs.wcs.cdelt = np.array([-pxscale, pxscale])/3600.
    wcs.wcs.crval = [cra, cdec]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crota = [-rotang, -rotang] # rotate SKY this amount? 
    # [This is consistent with original RSMT angle definition]
    wcs.wcs.equinox=equinox

#    wcs.wcs.print_contents()
    #ra,dec = wcs.wcs_pix2sky(xpix,ypix, 1)  

    # -- convert central coords to Ra, dec and linear sizes to angular sizes:     
    ccd1_sky=[wcs.wcs_pix2sky(np.reshape(np.array(ccd1[0]),(1,2)), 1),ccd1[1]*pxscale,ccd1[2]*pxscale ]
    ccd2_sky=[wcs.wcs_pix2sky(np.reshape(np.array(ccd2[0]),(1,2)), 1),ccd2[1]*pxscale,ccd2[2]*pxscale ]
    ccd3_sky=[wcs.wcs_pix2sky(np.reshape(np.array(ccd3[0]),(1,2)), 1),ccd3[1]*pxscale,ccd3[2]*pxscale ]
    circ_fov_sky=[wcs.wcs_pix2sky(np.reshape(np.array(circ_fov[0]),(1,2)), 1),circ_fov[1]*pxscale]
    
    
    return ccd1_sky,ccd2_sky,ccd3_sky,circ_fov_sky


def DefineCCDsPix(cra,cdec,rotang,equinox):
    #rect=Rectangle((0,0),ccd_dx,ccd_dy,color='none',ec='y')
    #ax.add_patch(rect)
    #rect=Rectangle((0+ccd_dx+ccd_xgap,0),ccd_dx,ccd_dy,color='none',ec='y')
    #ax.add_patch(rect)
    #rect=Rectangle((0+2.*(ccd_dx+ccd_xgap),0),ccd_dx,ccd_dy,color='none',ec='y')
    #ax.add_patch(rect)
    
    # Define lower left corner, + width and height
    
    # CCD 1:
#    offset=[ccd_dx/2.,ccd_dy/2.] # bottom left --> centre of chip
#    ccd1=[ [0,0],ccd_dx,ccd_dy ]
#    ccd2=[ [0+ccd_dx+ccd_xgap,0],ccd_dx,ccd_dy]
#    ccd3=[[0+2.*(ccd_dx+ccd_xgap),0],ccd_dx,ccd_dy]   
    ccd1=[ [0+ccd_dx/2.,0+ccd_dy/2.],ccd_dx,ccd_dy ]
    ccd2=[ [0+ccd_dx+ccd_xgap+ccd_dx/2.,0+ccd_dy/2.],ccd_dx,ccd_dy]
    ccd3=[[0+2.*(ccd_dx+ccd_xgap)+ccd_dx/2.,0+ccd_dy/2.],ccd_dx,ccd_dy]   
#    ccd1=[ (-(ccd_dx+ccd_xgap),0),ccd_dx,ccd_dy ]
#    ccd2=[ (0.,0.),ccd_dx,ccd_dy]
#    ccd3=[((ccd_dx+ccd_xgap),0),ccd_dx,ccd_dy]   
    circ_fov=[(ccd_cx,ccd_cy),dcr*3600./pxscale] # everything in pix
#    circ_fov=[(0.,0.),dcr*3600./pxscale] # everything in pix

    return ccd1,ccd2,ccd3,circ_fov


def toRegionSky(sccd1,sccd2,sccd3,scirc_fov,rotang,regfile):
    regstr='FK5;'    
#    print np.shape(sccd1[0])
#    print (sccd3[0][0,0], sccd3[0][0,1], sccd3[1], sccd3[2], rotang)
    xxx=open(regfile,'w')
    regstr+='box(%s, %s, %s\", %s\", %s)\n' % \
            (sccd1[0][0,0], sccd1[0][0,1], sccd1[1], sccd1[2], rotang)
    xxx.write(regstr)
    regstr='FK5;'    
    regstr+='box(%s, %s, %s\", %s\", %s)\n' % \
            (sccd2[0][0,0], sccd2[0][0,1], sccd2[1], sccd2[2], rotang)
    xxx.write(regstr)
    regstr='FK5;'    
    regstr+='box(%s, %s, %s\", %s\", %s)\n' % \
            (sccd3[0][0,0], sccd3[0][0,1], sccd3[1], sccd3[2], rotang)
    xxx.write(regstr)   
    regstr='FK5;'    
    regstr+='circle(%s, %s, %s\")\n' % (scirc_fov[0][0,0],scirc_fov[0][0,1], scirc_fov[1]) 
    xxx.write(regstr)
    xxx.close()
    #return regstr

if __name__ == '__main__':
    cra=53.934130928509504
    cdec=-27.04987517144887
    rotang=70.
    equinox=2000.0
    
    # -- detector footprints, pix.:
    # [This rotang is consistent with original RSMT file angle definition.]
    det1,det2,det3,detc=DefineCCDsPix(cra,cdec,rotang,equinox)
    #print det1,det2,det3,detc
    # -- detector footprints, sky coor.:
    skyd1,skyd2,skyd3,skycirc = RSSskyfromPix(cra,cdec,rotang,equinox,det1,det2,det3,detc)
    #print skyd1,skyd2,skyd3,skycirc
    
    # i think the -ve sign here is right. It needs to be opposite to the angle in the CD matrix, anyway
    # **** CHECK this agrees with ROTANG defined for telescope ****
    #print 
    
    # -- write region file
    regfile='test.reg'
    toRegionSky(skyd1,skyd2,skyd3,skycirc,rotang,regfile)
    #fk5;circle(35.458437,-3.7601347,5.3749084")

# ============================================================================================================================

def FOVTest(cra,cdec,equinox,rotang,slitra,slitdec,slit_length,tilt):
    # ---- Test if slit lies entirely within RSS FOV
    
    # If the top or bottom edge of the slit lies outside the nominal FOV, flag as outside
    
    # -- convert sky coords to pix:    
    wcs = pywcs.WCS(naxis=2)
    wcs.wcs.crpix = [ccd_cx,ccd_cy] 
#    wcs.wcs.crpix = [0,0] # define centre relative to zero

    #wcs.wcs.cdelt = np.array([-pxscale, pxscale])
    wcs.wcs.cdelt = np.array([-pxscale, pxscale])/3600.
    wcs.wcs.crval = [cra, cdec]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crota = [rotang, rotang] # rotate SKY this amount? 
    wcs.wcs.equinox=equinox

    xp,yp = wcs.wcs_sky2pix(slitra,slitdec, 1)
    
    # since we are only supporting rotang = 0/180 so far:
    ytop=yp+(slit_length/2.0/pxscale) 
    ybot=yp-(slit_length/2.0/pxscale) 
    # assume slit width is negligible. could add small safety padding to dcr
    
    # --test perimeter of field:
    fov_flag = np.zeros(len(cra))
    # 0 means in FOV, 1 means ootside:
    maxdis2 = dcr*3600.0/pxscale
    topdis2 = (ytop-ccd_cy)**2 + (xp-ccd_cx)**2
    botdis2 = (ybot-ccd_cy)**2 + (xp-ccd_cx)**2
    
    #oof=np.reshape((topdis2>=maxdis2).nonzero(),-1)
    oof=np.reshape((topdis2>=maxdis2).nonzero(),-1)
    fov_flag[oof]=1
    oof=np.reshape((botdis2>=maxdis2).nonzero(),-1)
    fov_flag[oof]=1
         
    return fov_flag

# =============================================================================================================================

def ConvertSlitLen(objra,objdec,slit_tilt,rotang,len1,len2):

    # ---- Convert user-defined object position and (possible) asymmetric slit lengths to a symmetrized slit centre:
    
    # take object coords, and length1, length2
    # **** only support rotang=0./180. and tilt=0. for now ****
    
#    from wiki"
#    Len1 -- (real) requested length above object (in direction of PA) in arcsec; default set by parameters
#     Len2 -- (real) requested length below object (opposite PA) in arcsec; default set by parameters
    slit_length=len1+len2 # symmetrized
    # calculate shift from midpoint (object coords):
    ddec = slit_length/2.0 - len1
    slitdec = objdec+ddec # ****
    slitra = objra # need to update for arbitrary angles later ****
     
    return slitra,slitdec

# =============================================================================================================================
    
def FindCollisions(slit_ra,slit_dec,slit_length,slit_width,cra,cdec,rotang,equinox,xpad,ypad,xleft,xright):

    idebug=0
    
    # xleft,xright are spectral lengths in pixels - set to large dummy values for now

# **** again, only does special case of PA=0/180, no slit tilt...

    from pyslit_optimize import RSSpixfromSky
    
    nslits=len(slit_ra)
    # -- identify collisions between slits and tag ID of colliding objects:
    coll_flag=['false']*nslits
    coll_ids=['none']*nslits ## can assign a list of collisions to any item in this list
    
    xp,yp=RSSpixfromSky(cra,cdec,rotang,slit_ra,slit_dec,equinox)
    
    # **** i think for general case i just need to replace slit width below with projected height which would be:
    #      total_height = slit_length*cos(theta) + 2.* 0.5 * slit_width*sin(theta)
    #         where theta is some combination of rotang and tilt 
    #
    # need to update optimiser for more general case, too.
    tx0=xp-xleft-xpad/2.
    tx1=xp+xright+xpad/2. # [ careful switching between lower left and width and lower left and upper right notation ]
    ty0=yp-(slit_length/2.)-ypad/2.
    ty1=yp+-(slit_length/2.)+ypad/2.

    # dumb but reliable way of checking:
    for ii in range(np.size(tx0)):

        if idebug: print ii
        
        tcoll_ids=[]
        for jj in range(np.size(tx0)):
            # -- check if this rectangle overlaps with any currently in mask:
            if ii == jj: continue # don't compare slit with itself
            if idebug: print 'comparing slit ',ii,' with ',jj
# http://tech-read.com/2009/02/06/program-to-check-rectangle-overlapping/
            r1x1=tx0[ii]
            r1x2=tx1[ii]
            r1y1=ty0[ii]
            r1y2=ty1[ii]
            r2x1=tx0[jj]
            r2x2=tx1[jj]
            r2y1=ty0[jj]
            r2y2=ty1[jj]

            # isOVerlap= ((r1x2 >= r2x1) &&
            #      (r1y2 >= r2y1) &&
            #      (r1x1 <= r2x2) &&
            #    (r1y1 <= r2y2));

            #print np.shape(r1x2),np.shape(r2x1)

            if ((r1x2 >= r2x1) and \
                    (r1y2 >= r2y1) and \
                    (r1x1 <= r2x2) and \
                    (r1y1 <= r2y2)) : coll_flag[ii]='true'
#            else: olap=0
            if idebug: print r1y1,r1y2,r2y1,r2y2
            
            #if (r1y2 >= r2y1) and (r1y1 <= r2y2) : olap=1

            tcoll_ids.append(jj)
    
        coll_ids[ii]=tcoll_ids
        
    # **** might want slit name rather than ID number ****
    return coll_flag,tcoll_ids
    