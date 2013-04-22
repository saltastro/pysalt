################################# LICENSE ##################################
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.                                                     #
#                                                                          #
# Redistribution and use in source and binary forms, with or without       #
# modification, are permitted provided that the following conditions       #
# are met:                                                                 #
#                                                                          #
#     * Redistributions of source code must retain the above copyright     #
#       notice, this list of conditions and the following disclaimer.      #
#     * Redistributions in binary form must reproduce the above copyright  #
#       notice, this list of conditions and the following disclaimer       #
#       in the documentation and/or other materials provided with the      #
#       distribution.                                                      #
#     * Neither the name of the South African Astronomical Observatory     #
#       (SAAO) nor the names of its contributors may be used to endorse    #
#       or promote products derived from this software without specific    #
#       prior written permission.                                          #
#                                                                          #
# THIS SOFTWARE IS PROVIDED BY THE SAAO ''AS IS'' AND ANY EXPRESS OR       #
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED           #
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE   #
# DISCLAIMED. IN NO EVENT SHALL THE SAAO BE LIABLE FOR ANY                 #
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL       #
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS  #
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)    #
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      #
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN #
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE          #
# POSSIBILITY OF SUCH DAMAGE.                                              #
############################################################################

"""Slotphottools contains the algorithms to perform photometry on slotmode data.

TODO--
* Change how the region file is read in so you can use different region formats
"""

# Ensure python 2.5 compatibility
from __future__ import with_statement

from pyraf import iraf
import saltsafekey as saltkey
import saltsafeio as saltio
import salttime
import saltstat
from salterror import SaltError, SaltIOError

import numpy as np

def calcstat(image, sig, iter):
    """Calculate the image statistics"""
    imsize=np.size(image)
    if imsize<=0: return 0

    mean=np.average(image)
    stddev=np.sqrt(np.sum((image-mean)**2)/imsize)
    for i in range(iter):
        mean=saltstat.mean2dclip(image,mean, stddev, sig)
        stddev=saltstat.std2dclip(image,mean, stddev, sig)

    return mean, stddev


def checkedge(x, x0, x1):
    """Determeine if the value x is over the edge
    of an image
    """
    x=int(x)
    if x < x0: return x0
    if x > x1: return x1
    return x

def findcentroid(array):
    """Find the centroid of an array where

        .. math:: cx=\\frac{\sum f\cdot x}{\sum f}

        returns :math:`cx,cy`.
    """

    cx=-1
    cy=-1

    try:
        i=len(array)
        j=len(array[0])
    except:
        return cx, cy

    y,x=np.indices((i,j))
    tot=array.sum()
    if tot!=0:
        cx=(array*x).sum()/tot
        cy=(array*y).sum()/tot

    return cx,cy

def objdetect(image, sigdet, stddev):
    """Detect all significant pixels in an image
    returns count
    """
    mask = (image > sigdet*stddev)
    return mask.sum()

def calcdrift(image, x, y, r, naxis1, naxis2):
    """Given a central x-y, calculate the drift"""
    x1=checkedge(x-r,0,naxis1)
    x2=checkedge(x+r,0,naxis1)
    y1=checkedge(y-r,0,naxis2)
    y2=checkedge(y+r,0,naxis2)
    cimage=image[y1:y2,x1:x2]
    fx,fy=findcentroid(cimage)
    if fx >= 0 and fy >= 0:
        cx=x1+fx
        cy=y1+fy
    else:
        cx=-1
        cy=-1
    return cimage, cx, cy

def findneareststar(array, xc, yc, rc, std, dsig, cpix, naxis1, naxis2):
    """Search an image for all stars and find the nearest one to the
    given coordinates

    returns an x,y
    """
    import scipy.ndimage as nd
    fx=-1
    fy=-1

    #detect objects and find their positions
    mask = (array > dsig*std)
    obj_labels, obj_num = nd.label(mask)
    if obj_num<1: return fx,fy

    coolist=nd.center_of_mass(array,obj_labels,range(obj_num))

    #find the closest coordinate to the previous coordinate
    if obj_num>1:
        best_rad=10000
        for y,x in coolist:
            rad=np.sqrt((x-xc)**2+(y-yc)**2)
            if rad < best_rad:
                best_rad=rad
                fx=x
                fy=y
    elif obj_num==1:
       try:
           fy=coolist[0]
           fx=coolist[1]
       except IndexError:
           #fix because of changed output in center of mass software
           fy=coolist[0][0]
           fx=coolist[0][1]
    else:
       return fx, fy


    #Confirm the source is a good source
    if fx > -1 and fy > -1:
        carray, fx, fy = calcdrift(array, fx, fy, rc, naxis1, naxis2)
        cimsize=np.size(carray)
        if cimsize < cpix:
            fx=-1
            fy=-1

    return fx,fy

def finddrift(image, cx, cy, cr, naxis1, naxis2, sigdet, cpix,  sigback, driftlimit, iter):
    """Starting with the initial x-y of the comparison source, search for
    the comparison star and see how much it has drifted.
    Expanding the search if no comparison star is found

    returns dx,dy
    """
    fx=-1
    fy=-1

    # set the image statistics for the image
    mean,stddev=calcstat(image,sigback,iter)

    # First to the fast search for the object, assuming it is roughly in the same position
    if driftlimit>0:
        dr=int(driftlimit)
    else:
        dr=cr
    cimage, cx2, cy2 = calcdrift(image, cx, cy, dr, naxis1, naxis2)
    cimsize=np.size(cimage)
    nobj=objdetect(cimage, sigdet, stddev)
    if cx2 >= 0 and cy2 >= 0 and nobj > cpix:
        fx=cx2
        fy=cy2
        # repeat several times to grab the best x-y
        for i in range(iter):
            cimage, fx, fy = calcdrift(image, fx, fy, dr, naxis1, naxis2)
        return cimage, fx, fy

    #set up the image to be searched
    if driftlimit>0:
        x1=checkedge(cx-dr,0,naxis1)
        x2=checkedge(cx+dr,0,naxis1)
        y1=checkedge(cy-dr,0,naxis2)
        y2=checkedge(cy+dr,0,naxis2)
        cimage=image[y1:y2,x1:x2]
    else:
        cimage=image.copy()

    cx2,cy2=findneareststar(cimage,cx,cy,dr,stddev, sigdet, cpix, naxis1,naxis2)
    if cx2 >= 0 and cy2 >= 0:
        fx=cx2+x1
        fy=cy2+y1
        #repeat several times to grab the best x-y
        for i in range(iter):
            cimage, fx, fy = calcdrift(image, fx, fy, cr, naxis1, naxis2)
        return cimage, fx, fy

    #if nothing was found, return values of -1
    return cimage, fx,fy

def getobstime(struct, infile):
    """Return the time since noon the previous day in seconds
    """
    time=-1
    hour,min,sec= saltkey.timeobs(struct,infile)
    time = salttime.utc2sec(hour,min,sec)
    if (hour > 12):
        time -= 43200.
    else:
        time += 43200.

    return time


def calcstarbackground(array, x, y, br1, br2, gain, rdnoise, naxis1, naxis2):
    """Determines the sky background around a single star.
    Using square annulli.
    Returns the background flux in terms of e/pix.

    In the case br1 and br2 are floats, x and y are the center of the
    annulli and br1 and br2 are floats which give the half-width of the
    square annulli. If br1 and br2 are lists, then the background rectangles
    are defined by br1 and br2.

    returns bflux, npix
    """

    # check for region/annulus
    try:
        if len(br1)==2 and len(br2)==2:
            mode='region'
        else:
            raise SaltIOError('Expected scalar or list of length 2, got lists of length '+str(len(br1))+' '+str(len(br2)))
    except TypeError:
        # Object is scalar as expected for annulus mode
        mode='annulus'
        
    if mode=='annulus':
        # calculate counts from first annullus
        bx1=checkedge(x-br1,0,naxis1)
        bx2=checkedge(x+br1,0,naxis1)
        by1=checkedge(y-br1,0,naxis2)
        by2=checkedge(y+br1,0,naxis2)
        bflux1=array[by1:by2,bx1:bx2].sum()
        barea1=np.size(array[by1:by2,bx1:bx2])

        # calculate counts from second annullus
        bx1=checkedge(x-br2,0,naxis1)
        bx2=checkedge(x+br2,0,naxis1)
        by1=checkedge(y-br2,0,naxis2)
        by2=checkedge(y+br2,0,naxis2)
        bflux2=array[by1:by2,bx1:bx2].sum()
        barea2=np.size(array[by1:by2,bx1:bx2])
        npix=barea2-barea1
    elif mode=='region':
        x1=checkedge(x-br1[0],0,naxis1)
        y1=checkedge(y-br1[1],0,naxis2)
        x2=checkedge(x+br2[0],0,naxis1)
        y2=checkedge(y+br2[1],0,naxis2)
        flux=array[y1:y2,x1:x2].sum()
        npix=np.size(array[y1:y2,x1:x2])
    else:
        raise SaltError('Selected sky background mode '+str(mode)+' not implemented.')

    # if enough pixels, calculate the sky background
    if mode=='annulus' and npix>0:
        bflux=(bflux2-bflux1)/(barea2-barea1)
        berr=np.sqrt(abs(bflux)+npix*rdnoise**2)/npix
    elif mode=='region' and npix>0:
        bflux=flux/npix
        berr=np.sqrt(abs(bflux)+npix*rdnoise**2)/npix
    else:
        bflux=0
        berr=0
        raise SaltError('Poor Background area specified, npix is '+str(npix))

    return bflux, berr, npix

def sqapphot(array, xc, yc, rc, br1, br2,gain, rdnoise,  naxis1, naxis2):
    """Perform square aperture photometry

    returns flux, flux_err
    """

    #first calculate the flux
    x1=checkedge(xc-rc,0,naxis1)
    x2=checkedge(xc+rc,0,naxis1)
    y1=checkedge(yc-rc,0,naxis2)
    y2=checkedge(yc+rc,0,naxis2)
    area=np.size(array[y1:y2,x1:x2])
    if area==0:
        flux=0
        flux_err=0
        raise SaltError('No data in the specified array')

    flux=array[y1:y2,x1:x2].sum()
    # calculate the background flux
    bflux, berr, npix=calcstarbackground(array, xc, yc, br1, br2, gain, rdnoise, naxis1, naxis2)

    # calculate the flux
    flux=flux-area*bflux
    flux_err=np.sqrt(abs(flux)+area*rdnoise**2+area*berr**2+area**2*berr**2/npix)

    return  flux, flux_err

def capphot(array, xc, yc, rc, br1, br2,gain, rdnoise,  naxis1, naxis2):
    """Perform circular aperture photometry

    return flux, flux_err
    """

    # calculate the flux
    try:
        i=len(array)
        j=len(array[0])
        y,x=np.indices((i,j))
        r=np.sqrt((x-xc)**2+(y-yc)**2)
        mask=(r<rc)
        area=mask.sum()
    except:
        area=0

    if area==0:
        flux=0
        flux_err=0
        raise SaltError('No data in the specified array')

    flux=(array*mask).sum()
    # calculate the background flux
    bflux, berr, npix=calcstarbackground(array, xc, yc, br1, br2, gain, rdnoise, naxis1, naxis2)

    # calculate the flux
    flux=flux-area*bflux
    flux_err=np.sqrt(abs(flux)+area*rdnoise**2+area*berr**2+area**2*berr**2/npix)

    return  flux, flux_err

def cogphot(array, xc, yc, rc, br1, br2,gain, rdnoise,  naxis1, naxis2):
    """Perform curve of growth photometry.
    Will use <rc> square apertures between 2 pix
    to the br2 radius.
    rc now gives the number of apertures to use.
    br1 gives the annulli to use for sky determination.
    The flux returned is the flux in the aperture with the highest signal to noise

    return flux, flux_err
    """

    min_rad=2
    rad=np.zeros(rc,dtype='float')
    flux=np.zeros(rc,dtype='float')
    flux_err=np.zeros(rc,dtype='float')
    rstep=(br2-min_rad)/rc

    best_sn=0
    best_i=-1
    for i in range(rc):
        rad[i]=min_rad+i*rstep
        flux[i], flux_err[i]=sqapphot(array, xc, yc, rad[i], br1, br2, gain, rdnoise, naxis1, naxis2)
        s_n=flux[i]/flux_err[i]
        if s_n > best_sn:
            best_sn=s_n
            best_i=i

    # if no value is found, use the middle aperture
    if best_i==-1: best_i=int(len(rad)/2)

    return rad[best_i], flux[best_i], flux_err[best_i]

def calc_optimal(array, xt, yt, xc, yc, rc, br1, br2,gain, rdnoise,  naxis1, naxis2):
    """Calculate the optimal ratio between the two sources based :math:`\\chi^{2}` pixel fit of the comparison star to the object star.
    For both stars the background will be removed and then the comparison star will be used as a 'model' psf for the :math:`\\chi^{2}` fit.
    To simplify the error calculation, we are assume the errors on the comparison star are much smaller than the errors on the target star

    return ratio, ratio_err
    """
    ratio=0
    ratio_err=0

    # find the background for both images
    btflux, bterr, ntpix=calcstarbackground(array, xt, yt, br1, br2, gain, rdnoise, naxis1, naxis2)
    bcflux, bcerr, ncpix=calcstarbackground(array, xc, yc, br1, br2, gain, rdnoise, naxis1, naxis2)

    # create the target array
    x1=checkedge(xt-rc,0,naxis1)
    x2=checkedge(xt+rc,0,naxis1)
    y1=checkedge(yt-rc,0,naxis2)
    y2=checkedge(yt+rc,0,naxis2)
    tarray=array[y1:y2,x1:x2]
    tarea=np.size(tarray)
    tarray=tarray-btflux

    # create the comparison array
    x1=checkedge(xc-rc,0,naxis1)
    x2=checkedge(xc+rc,0,naxis1)
    y1=checkedge(yc-rc,0,naxis2)
    y2=checkedge(yc+rc,0,naxis2)
    carray=(array[y1:y2,x1:x2])
    carea=np.size(carray)
    carray=carray-bcflux

    # check to make sure the objects are the same size
    if np.size(tarray)>np.size(carray):
        x1=checkedge(xc-rc,0,naxis1)
        x2=checkedge(xc+rc,0,naxis1)
        y1=checkedge(yc-rc,0,naxis2)
        y2=checkedge(yc+rc,0,naxis2)
        dx1=int(xc-x1)
        dx2=int(x2-xc)
        dy1=int(yc-y1)
        dy2=int(y2-yc)
        x1=checkedge(xt-dx1,0,naxis1)
        x2=checkedge(xt+dx2,0,naxis1)
        y1=checkedge(yt-dy1,0,naxis2)
        y2=checkedge(yt+dy2,0,naxis2)
        tarray=array[y1:y2,x1:x2]
        tarea=np.size(tarray)
        tarray=tarray-btflux
    if np.size(tarray)<np.size(carray):
        x1=checkedge(xt-rc,0,naxis1)
        x2=checkedge(xt+rc,0,naxis1)
        y1=checkedge(yt-rc,0,naxis2)
        y2=checkedge(yt+rc,0,naxis2)
        dx1=int(xt-x1)
        dx2=int(x2-xt)
        dy1=int(yt-y1)
        dy2=int(y2-yt)
        x1=checkedge(xc-dx1,0,naxis1)
        x2=checkedge(xc+dx2,0,naxis1)
        y1=checkedge(yc-dy1,0,naxis2)
        y2=checkedge(yc+dy2,0,naxis2)
        carray=array[y1:y2,x1:x2]
        carea=np.size(carray)
        carray=carray-bcflux

    # check to make sure they are the same size
    # and fail if they are not
    if carea!=tarea:
        raise SaltError('Area do not match.')

    # find the ratio between the two
    psum=(carray*carray).sum()
    if psum == 0:
        ratio = 0
        ratio_err=1
    else:
        ratio = (tarray*carray).sum()/psum
        tflux,terr=sqapphot(array, xt, yt, rc, br1, br2,gain, rdnoise,  naxis1, naxis2)
        cflux,cerr=sqapphot(array, xc, yc, rc, br1, br2,gain, rdnoise,  naxis1, naxis2)
        if not tflux==0 and not cflux==0:
            ratio_err=ratio*np.sqrt((terr/tflux)**2+(cerr/cflux)**2)
        else:
            ratio_err=0

    return ratio, ratio_err

def dophot(phottype, array, x, y, r, br1, br2, gain, rdnoise, naxis1, naxis2):
    """Routine to perform photometry on the image.
    It will perform the requested type of photometry (square aperture, circular aperture, curve of growth, psf fit, or optimal).
    See the help file for descriptions of each.

    returns list of photometry parameters for target and companion
    """
    time=1

    if phottype=='square':
        target_flux, target_flux_err=sqapphot(array, x['target'], y['target'], r['target'],br1['target'], br2['target'], gain, rdnoise,  naxis1, naxis2)

        compar_flux, compar_flux_err=sqapphot(array, x['comparison'], y['comparison'], r['comparison'], br1['comparison'], br2['comparison'], gain, rdnoise,  naxis1, naxis2)

        ratio, ratio_err=calcratio(target_flux, target_flux_err, compar_flux, compar_flux_err)

    elif phottype=='circular' :
        target_flux, target_flux_err=capphot(array, x['target'], y['target'], r['target'], br1['target'], br2['target'], gain, rdnoise,  naxis1, naxis2)

        compar_flux, compar_flux_err=capphot(array, x['comparison'], y['comparison'], r['comparison'], br1['comparison'], br2['comparison'], gain, rdnoise,  naxis1, naxis2)

        ratio, ratio_err = calcratio(target_flux, target_flux_err, compar_flux, compar_flux_err)

    elif phottype=='cog':
        r_cog, compar_flux, compar_flux_err=cogphot(array, x['comparison'], y['comparison'], r['comparison'], br1['comparison'], br2['comparison'], gain, rdnoise,  naxis1, naxis2)

        target_flux, target_flux_err=sqapphot(array, x['target'], y['target'], r_cog, br1['target'], br2['target'], gain, rdnoise,  naxis1, naxis2)

        ratio, ratio_err = calcratio(target_flux, target_flux_err, compar_flux, compar_flux_err)

    elif phottype=='psf':
        compar_flux, compar_flux_err=sqapphot(array, x['comparison'], y['comparison'], r['comparison'], br1['comparison'], br2['comparison'], gain, rdnoise,  naxis1, naxis2)

        target_flux, target_flux_err=sqapphot(array, x['target'], y['target'], r['target'], br1['target'], br2['target'], gain, rdnoise,  naxis1, naxis2)
    elif phottype=='optimal':
        compar_flux, compar_flux_err=sqapphot(array, x['comparison'], y['comparison'], r['comparison'], br1['comparison'], br2['comparison'], gain, rdnoise,  naxis1, naxis2)

        ratio, ratio_err=calc_optimal(array, x['target'], y['target'], x['comparison'], y['comparison'], r['comparison'], br1['comparison'], br2['comparison'], gain, rdnoise,  naxis1, naxis2)

        target_flux=ratio*compar_flux

        if not ratio==0 and not compar_flux:
            target_flux_err=target_flux*np.sqrt((compar_flux_err/compar_flux)**2+(ratio_err/ratio)**2)
        else:
            target_flux_err=compar_flux_err
    else:
        raise SaltError('Photometry Type does not exist')

        return 0,0,0,0,0,0

    # place the appropriate information into pinfo
    return target_flux, target_flux_err, compar_flux, compar_flux_err, ratio, ratio_err

def calcratio(target_flux, target_flux_err, compar_flux, compar_flux_err):
    """Calculate the flux ratio

    return ratio, ratio_err
    """
    if compar_flux>0:
        ratio=target_flux/compar_flux
    else:
        ratio=0

    if compar_flux > 0 and target_flux > 0:
        ratio_err=ratio*np.sqrt((target_flux_err/target_flux)**2+(compar_flux_err/compar_flux)**2)
    else:
        ratio_err=ratio*np.sqrt(target_flux_err**2+compar_flux_err**2)

    return ratio, ratio_err

def writedataout(lout, index, time, x, y, target_flux, target_flux_err, compar_flux, compar_flux_err, ratio, ratio_err, time0, reltime):
    """Write the data out to a file"""

    #calculate the ratio of the target to comparison star
    target_x=x['target']
    target_y=y['target']
    compar_x=x['comparison']
    compar_y=y['comparison']

    #if reltime is selected, only write out the time since the first image
    if reltime:
        wtime=time-time0
    else:
        wtime=time

    #write out the statement
    statement = '%5i %10.3f %10.6f %7.6f %5.2f %5.2f %10.3f %7.3f %5.2f %5.2f %10.3f %7.3f  \n' % \
         (index, wtime, ratio, ratio_err, target_x, target_y, target_flux, target_flux_err, compar_x, compar_y, compar_flux, compar_flux_err)
    try:
        lout.write(statement)
    except:
        raise SaltIOError('Unable to write to the output file')

def readsrcfile(srcfile): 
    """Read in the src file and return the parameters extracted from that file

    """
    #set up some variables
    amp = {}
    x = {}
    y = {}
    x_o = {}
    y_o = {}
    r = {}
    br1 = {}
    br2 = {}
    tgt_btype='region' 
    cmp_btype='region'
 
    # check extraction region defintion file exists
    srcfile = srcfile.strip()
    saltio.fileexists(srcfile)
    
    #try loading the region file
    try:
        region=np.round(np.loadtxt(srcfile))
    except Exception, e:
        msg='SLOTTOOLs--ERROR:  Unable to load array from %s' % srcfile
        raise SaltIOError(msg)

    # read extraction region defintion file
    try:

        # Check if region definition file has an acceptable format
        if region.shape[1] not in [7,9]:
            raise SaltIOError('Region definition file has invalid number of columns')
         

        # Check if all object numbers in the region file are unique
        if len(set(region[:,0]))!=len(region[:,0]):
            raise SaltIOError('Object numbers in region definition file not unique')

        # Find target and comparison
        tindex=None
        cindex=None
        for i in range(region.shape[0]):
            if region[i][0]==1:
                tindex=i
            elif region[i][0]==2:
                cindex=i
            else:
                msg='Object %d not used' % region[i][0]
                log.message(msg)
                continue
        if tindex is None or cindex is None:
            raise SaltIOError('Could not find target and comparison in region definition file.')

        #set the background type
        if len(region[tindex])==7:  tgt_btype='annulus'
        if len(region[cindex])==7:  cmp_btype='annulus'
        print tgt_btype, cmp_btype
        # Check if selected background mode matches defined regions
        valid=True
        if (tgt_btype=='region' or cmp_btype=='region') and region.shape[1]<9:
            valid=False
        if tgt_btype=='annulus' and not np.any(region[tindex][5:7]):
            valid=False
        elif tgt_btype=='region' and not np.any(region[tindex][5:9]):
            valid=False
        if cmp_btype=='annulus' and not np.any(region[cindex][5:7]):
            valid=False
        elif cmp_btype=='region' and not np.any(region[cindex][5:91]):
            valid=False
        if not valid:
            raise SaltIOError('Selected background types do not match regions specified in region definition file.')

        # Read the parameters
        for object in [('target',tindex),('comparison',cindex)]:
            amp[object[0]]=int(region[object[1]][1])
            x[object[0]]=int(region[object[1]][2])
            y[object[0]]=int(region[object[1]][3])
            x_o[object[0]]=int(region[object[1]][2])
            y_o[object[0]]=int(region[object[1]][3])
            r[object[0]]=int(region[object[1]][4])

        if tgt_btype=='annulus':
            # Ensure br2>=br1
            br1['target']=np.min((int(region[tindex][5]),int(region[tindex][6])))
            br2['target']=np.max((int(region[tindex][5]),int(region[tindex][6])))
        elif tgt_btype=='region':
            i=tindex
            # Ensure vertices are in proper order
            br1['target']=[np.min((int(region[i][5]),int(region[i][7]))),
                           np.min((int(region[i][6]),int(region[i][8])))]
            br2['target']=[np.max((int(region[i][5]),int(region[i][7]))),
                           np.max((int(region[i][6]),int(region[i][8])))]
        if cmp_btype=='annulus':
            # Ensure br2>=br1
            br1['comparison']=np.min((int(region[cindex][5]),int(region[cindex][6])))
            br2['comparison']=np.max((int(region[cindex][5]),int(region[cindex][6])))
        elif cmp_btype=='region':
            i=cindex
            # Ensure vertices are in proper order
            br1['comparison']=[np.min((int(region[i][7]),int(region[i][5]))),
                           np.min((int(region[i][8]),int(region[i][6])))]
            br2['comparison']=[np.max((int(region[i][7]),int(region[i][5]))),
                           np.max((int(region[i][8]),int(region[i][6])))]

    except Exception, e:
        msg='Error processing region definition file because %s'%e
        raise SaltIOError(msg)

    if (len(amp) != 2 or len(x) != 2 or len(y) != 2 or len(r) != 2 or len(br1) != 2 or len(br2) != 2):
            message = 'ERROR: SALTVIEW -- extraction region definition file '+srcfile
            message += 'does not have an acceptable format.'
            raise SaltIOError(msg)


    return amp, x, y, x_o, y_o, r, br1, br2

def readlcfile(lcfile):
    """Read in the lightcurve file

    """
    saltio.fileexists(lcfile)
    
    try:
        id, time, ratio, rerr, tx, ty, tflux, terr, cx, cy, cflux, cerr=np.loadtxt(lcfile, unpack=True)
    except Exception, e:
        msg='Error processing light curver file because %s' % e
        raise SaltIOError(msg)

    return id, time, ratio, rerr, tx, ty, tflux, terr, cx, cy, cflux, cerr
