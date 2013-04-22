
#!/usr/bin/python


# # ^^^^^^ -- indicates parameters to change

import numpy as np
import pylab as pl
import pywcs
import pickle

##
#from matplotlib.patches import Rectangle,Circle
##

pxscale = 0.2507 / 2. # unbinned

dcr = 4. / 60. # radius of field (deg)
#dcr=3.9/60. # radius of field (deg)
# shrink radius slightly to avoid part of slit being drawn outside of mask circle?

iplot = 0 #1
idebug = 0

# ^^^^^^
#niter=10


# global CCD parameters:
ccd_dx = 2034.
ccd_xgap = 70.
ccd_dy = 4102.

ccd_cx = (2. * (ccd_dx + ccd_xgap) +ccd_dx) / 2.
ccd_cy = ccd_dy / 2.



def RSSpixfromSky(cra,cdec,rotang,ra,dec,equinox):
    # -- Make a WCS to convert to RSS pixel coords
    # Create a new WCS object.  The number of axes must be set
    # from the start
    wcs = pywcs.WCS(naxis=2)
#    wcs.wcs.crpix = [3000.,2000.] # made up !! ****
    wcs.wcs.crpix = [ccd_cx,ccd_cy] 

    #wcs.wcs.cdelt = np.array([-pxscale, pxscale])
    wcs.wcs.cdelt = np.array([-pxscale, pxscale])/3600.
    wcs.wcs.crval = [cra, cdec]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
#    wcs.wcs.crota = [rotang, rotang] #** CHECK **
    wcs.wcs.crota = [-rotang, -rotang] # rotate SKY this amount? 
    # [This is consistent with original RSMT angle definition]
    wcs.wcs.equinox=equinox

#    wcs.wcs.print_contents()

    xpix,ypix = wcs.wcs_sky2pix(ra,dec, 1)  
    return xpix,ypix

def defineSpec(xobj, yobj, grating, filt, slen):
    # -- Use optical model to obtain wavelength limits (converted to pixels)
    # based on grating and filter selected. For now just make something up
    #xcen=3000. #**** made-up!
    xspec = xobj #(xobj-xcen)+(-0.2)*(xobj-xcen)+xcen # make up anamorphic factor
    yspec = yobj # make up anamorphic factor

    #^^^^^^^^
    # change by hand for grating+filter - will be set automatically eventually
    #xlo=400.  # spectral length to left of slit (unbinned pix)
    #xhi=150.  # ---"" ----          right ---""--
    xlo = 3000.
    xhi = 3000.



    # -- start off with ylo/yhi half slit length
    yhi = slen / 2.
    ylo = slen / 2.

    return xspec,yspec,xlo,xhi,ylo,yhi


def MonteRealise(x0, y0, x1, y1, pri):
    # ---- Monte Carlo assignment of weights, one realisation
    nums = np.random.random(np.shape(x0))
    keep = np.zeros(np.shape(x0)).astype('int')
    ok = (nums < pri).nonzero()
    keep[ok] = 1
    return keep


def ResolveCollisions(xspe, yspe ,x0, y0, x1, y1, tkeep, cx, cy, xpad, ypad, \
                        allocated, prisort):
    # ---- Make a first pass mask starting with object closest to central 
    #      and resolving collisions

    #**** check aboot skipping pri=-1. (setup star) obejcts ****

    #^^^^^^^^
    #idebug=0#1
    ivertdis=1 # use vertical distance (best for one tier of spectra)
    # otherwise use radial distance when distance sorting (best for multi-tier?)

#    if prisort: print "::::::::",pri
#    if prisort: print "::::::::",tkeep
        
    if prisort: # keep is priority instead of keepflag, sort by this
        tpri = tkeep
        ss = (np.argsort(tpri))[::-1] # reverse sort

        if idebug: print tpri[ss]
        tkeep = np.ones(np.shape(tpri)).astype('int')
        rej = np.reshape((tpri == 0.).nonzero(), -1)
        tkeep[rej] = 0
    else:
        if ivertdis:
            dis2 = (yspe-cy)**2
            ss = np.argsort(dis2)
        else:
            dis2 = (xspe-cx)**2 + (yspe-cy)**2
            ss = np.argsort(dis2)

    
    tx0 = x0[ss] - xpad / 2.
    tx1 = x0[ss] + x1[ss] + xpad / 2. # [ careful switching between lower left and width and lower left and upper right notation ]
    ty0 = y0[ss] - ypad / 2.
    ty1 = y0[ss] + y1[ss] + ypad / 2.
    skeep =tkeep[ss].copy()

    inmask = np.zeros(np.shape(x0)).astype('int')
    
    # -- add in already allocated slits
    ok = np.reshape((allocated[ss] == 1).nonzero(), -1)
    inmask[ok] = 1
    # --

    for ii in range(np.size(x0)):
        if skeep[ii] == 0:
            olap = 0
            continue

        if idebug: print ii

        # keep is either priority or always 1, depending on mode
#        if skeep[ii] <= 0: continue # don't care about pri=0. objects or collisions of ref stars with others
        if prisort:
            if tpri[ss[ii]] <= 0:
                continue # don't care about pri=0. objects or collisions of ref stars with others

        #--
#        if inmask[ii]==0 and keep[ii]<0.: continue # don't add a new setup star 
        #--


        if inmask[ii] == 1:
            if idebug: print ii,' already in mask'
            continue # no need to check if this slit collides
        used=np.reshape((inmask == 1).nonzero(), -1)
        if idebug:
            print np.size(used),' slits currently in mask'
        olap = 0
        for jj in range(np.size(used)):
            # -- check if this rectangle overlaps with any currently in mask:
            if idebug:
                print 'comparing slit ',ii,' with ',used[jj]
# http://tech-read.com/2009/02/06/program-to-check-rectangle-overlapping/
            r1x1 = tx0[ii]
            r1x2 = tx1[ii]
            r1y1 = ty0[ii]
            r1y2 = ty1[ii]
            r2x1 = tx0[used[jj]]
            r2x2 = tx1[used[jj]]
            r2y1 = ty0[used[jj]]
            r2y2 = ty1[used[jj]]

            # isOVerlap= ((r1x2 >= r2x1) &&
            #      (r1y2 >= r2y1) &&
            #      (r1x1 <= r2x2) &&
            #    (r1y1 <= r2y2));

            #print np.shape(r1x2),np.shape(r2x1)

            if ((r1x2 >= r2x1) and \
                    (r1y2 >= r2y1) and \
                    (r1x1 <= r2x2) and \
                    (r1y1 <= r2y2)) : olap=1
#            else: olap=0
            if idebug:
                print r1y1,r1y2,r2y1,r2y2
            
            #if (r1y2 >= r2y1) and (r1y1 <= r2y2) : olap=1


            if idebug:
                print 'olap: ',olap
# as soon as this slit overlaps with any one already in mask, we can reject it
            if olap == 1:
                break # hopefully just exits this inner loop
        # if we checked against all slits on mask and no collisions, then keep:
        if olap==0: 
            inmask[ii]=1
#            if prisort: 
#                print "adding  new slit:"
#                print ii,skeep[ii],ty1[ii]-ty0[ii]
#                print ii,tpri[ss[ii]],ty1[ii]-ty0[ii]
#                print

    #*** careful at end. inmask index is aligned with sorted-distance coords ***
    keepflag = np.zeros(np.shape(x0))
    ok = np.reshape((inmask == 1).nonzero(),-1)
    keepflag[ss[ok]] = 1
    return keepflag
    
def firstpass(Nstars_req, xspe, yspe, x0, y0, x1, y1, keep, cx, cy, xpad, \
                ypad, pri):
    allocated=np.zeros(np.shape(x0)).astype('int') # none already allocated
    # -- Now add setup *before* any science slits, otherwise each star wipes out 
    #     a vertical strip equivalent to the length of a science slit!
    allocated=setupStars(Nstars_req, xspe, yspe, x0, y0, x1, y1, pri, cx, cy, \
                xpad, ypad, allocated)
    res=ResolveCollisions(xspe, yspe, x0, y0, x1, y1, keep, cx, cy, xpad, \
                ypad, allocated,0)
    return res

def addmore(xspe, yspe, x0, y0, x1, y1, pri, cx, cy, xpad, ypad, in1stmask):
    # -- add more slits from full list (pri>0) after firstpass:
#    keep=np.zeros(np.size(pri))
#    ok=np.reshape((pri>0.).nonzero(),-1)
#    keep[ok]=1
    # sort by priority, so when turning off sorting, highest priorities are assigned first
    res=ResolveCollisions(xspe, yspe, x0, y0, x1, y1, pri, cx, cy, xpad, ypad, \
                in1stmask, 1)
    return res

def tweakslits(xspe, yspe, x0, y0, x1, y1, pri, cx, cy, xpad, ypad, inprevmask):
    # ******** TODO ********
    # -- add more slits from full list (pri>0) after firstpass:
#    keep=np.zeros(np.size(pri))
#    ok=np.reshape((pri>0.).nonzero(),-1)
#    keep[ok]=1
    # sort by priority, so when turning off sorting, highest priorities are assigned first
    y1sh = y1 * 0.8
    y0sh = y0 + 0.1 * y1 # make 80% slit length, but keep same centre
    res = ResolveCollisions(xspe, yspe, x0, y0sh ,x1, y1sh, pri, cx, cy, xpad, \
            ypad, inprevmask, 1)
    return res


def bestStars(sx,sy):
    # choose "optimal" stars from a larger list.
    # take 4-6 with best distribution across field
    
    # assume all have already been pre-selected to lie in a suitable mag range
    #**** dummy - RANDOM ****
    inds = np.random.random(np.shape(sx))
    outind = np.argsort(inds)[0:6]
#    print outind

    return outind


def setupStars(Nstars_req,xspe,yspe,x0,y0,x1,y1,pri,cx,cy,xpad,ypad,inprevmask):
    # for now just do this in a dumb way. Design the science mask first, 
    # then throw in setup stars, removing science slits where collisions occur. 
    # Probably best to do this after "firstmask" and before "addmore"

    # **** NOTE: need to deal with "musthave" objects correctly here and in 
    # firstmask, etc. ****

    # -- This should run like a simplified version of ResolveCollisions

    # check initial star list. This should be objects with priority=-1
    #****
    # (i think this will just work naturally in the science slit allocation)
    #****

#    print np.sum(inprevmask).astype('int'),"slits before setup stars added"

    stars = np.reshape((pri==-1.).nonzero(),-1)
    nstars = np.size(stars)
#    print 'NSTARS = ',nstars
#    if nstars < 4:
#        print "NO SETUP STARS IN CATALOGUE"
#        print "You must select these manually"
#        return inprevmask
    #if nstars > 6:
    if nstars > Nstars_req:
        # select best stars:
        tusestars = bestStars(x0[stars],y0[stars])
        usestars = stars[tusestars]
        nstars = np.size(usestars)
    else:
        usestars = stars
#    print 'N_USESTARS = ',np.size(usestars)
#    print usestars


    # == Setup stars use 1" diameter holes
    # reset width and length of these slit
    #x1[usestars]=1./pxscale # spectral length not slit width??
    #y1[usestars]=1./pxscale
            # ** set in input catalogue now **

    skeep = inprevmask

    skeep[usestars] = 1 # add stars to mask
    #-- increase priority so stars override science objects
#    print pri
#    print (pri==99).nonzero()    
    pri[usestars] = 10.
#    print (pri==10).nonzero()    
    #--


    tx0 = x0 - xpad / 2.
    tx1 = x0 + x1 + xpad / 2. # [ careful switching between lower left and width and lower left and upper right notation ]
    ty0 = y0 - ypad / 2.
    ty1 = y0 + y1 + ypad / 2.

    # for each science slit already in mask, check if it hits a setup star
    # if so, remove the former:
    for ii in range(np.size(x0)):
        if inprevmask[ii] != 1: # only care about collisions with allocated objects
            continue

        if pri[ii] == 10:
            continue # this is a star itself!

        for jj in range(nstars):
            # -- check if this rectangle overlaps with any currently in mask:
            if idebug: print 'comparing slit ',ii,' with ',usestars[jj]
            r1x1 = tx0[ii]
            r1x2 = tx1[ii]
            r1y1 = ty0[ii]
            r1y2 = ty1[ii]
            r2x1 = tx0[usestars[jj]]
            r2x2 = tx1[usestars[jj]]
            r2y1 = ty0[usestars[jj]]
            r2y2 = ty1[usestars[jj]]

            if ((r1x2 >= r2x1) and \
                    (r1y2 >= r2y1) and \
                    (r1x1 <= r2x2) and \
                    (r1y1 <= r2y2)) : 
                skeep[ii] = 0 #olap=1
                continue # once collided, no need to consider further
#            else: olap=0
            if idebug: print r1y1,r1y2,r2y1,r2y2
            if idebug: print "removing science slit due to setup star collision"
            if idebug: 
                if skeep[ii] == 0: print "removing science slit due to setup star collision"
    # summarise results
#    print np.sum(skeep).astype('int'),"slits afterwards, with ",nstars,"setup stars"
    
#    print (pri==10).nonzero()
#    print skeep[np.reshape((pri==10).nonzero(),-1)]
    pri[np.reshape((pri == 10).nonzero(), -1)] = -1. # this is a global variable, so need to rest after each iter!

    return skeep



# ==============================================================================
# ---- design an RSS mask

# added option to select Nstars_req from pri=-1. objects.
# needs to be some external check that this is not larger than number in catalogue! 
def pyslit_optimize(cra, cdec, rotang, equinox, ra, dec, pri, slen_as, swid_as, \
            stilt, Niter, opt_ypad, Nstars_req=0):

    slen = slen_as / pxscale
    swid = swid_as / pxscale


    # ---- Convert sky coords to RSS pixel coords (create dummy WCS for now)
    # These are coordinates for the object. Need to keep separate positions for the ref. wavelength under slit position due to anamorphic magnification
    xobj, yobj = RSSpixfromSky(cra, cdec, rotang, ra, dec, equinox)

    # ---- set-up properties of spectra in pixel coords
    # need lengths of spectra left and right, above and below object
    grating = 'fred'
    filt = 'harry'
    xspe, yspe, xle, xri, yab, ybe = defineSpec(xobj, yobj, grating, filt, slen)



    # **** Need to add a check for desired wavelength range staying on detector ****


    x0 = xspe - xle
    x1 = np.zeros(np.size(x0)) + xle + xri
    y0 = yspe - ybe
    y1 = ybe + yab
    #for kk in np.arange(np.size(x0)):
    #    if pri[kk]<1.0: col='r'
    #    if pri[kk]<0.7: col='b'
    #    if pri[kk]<0.5: col='g'
    #    
    #    rect=Rectangle((x0[kk],y0[kk]),x1[kk],y1[kk],color='none',ec=col)
    #    ax.add_patch(rect) 
    #pl.draw()


    # /////
    # pri=-1 -- setup star
    # specify nsetupstars. if there are more p=-1 objs than this, choose in some optimum way 
    # 
    # pri=9 is a must-have object which specifies point in mask to sort in distance from when optimising (not nec. cra,cdec)
    # if not set, sort from centre of mask
    # other must-have objects just p=1.
    #////

    #***need to set safety margins too. just extend slit length and shrink back later?

    # -- construct first pass realisation, just keep subsample of full object list
    #    based on priorities. No checking for slit collisions yet.


    maxwt = 0

    for mm in range(Niter):
#        print 'iter:',(mm+1)
        keep = MonteRealise(x0, y0, x1, y1, pri)


        # -- Now, order these by distance from desired object and remove colliding slits

        # -- centre of mask if not set to specific object
        #cx = 3000.
        #cy = 2000.
        cx = ccd_cx
        cy = ccd_cy

        xpad = 1. / pxscale # arcsec
        ypad = opt_ypad / pxscale
        in1mask = firstpass(Nstars_req, xspe, yspe, x0, y0, x1, y1, keep, cx, \
                cy, xpad, ypad, pri)


        ok = np.reshape((pri > 0.).nonzero(), -1)
#        print "First pass allocated",np.sum(in1mask[ok]).astype('int'),"science slits"\
#            ' with a total weight of %.1f'%np.sum(in1mask[ok]*pri[ok])
    #        ' with a total weight of %.1f'%np.sum(in1mask*pri)

        #--
    ##    inwstars=setupStars(xspe,yspe,x0,y0,x1,y1,pri,cx,cy,xpad,ypad,in1mask)
    ##    in1mask=inwstars
    #    sss=np.reshape( ((inwstars == 1) & (pri==-1)).nonzero() ,-1)
    #    print '>>> after setupStars ',np.size(sss)

        #--

        # -- add more slits from the full list (pri>0.) before 1st MC weights were applied
        bssss = np.reshape( ((in1mask == 1) & (pri==-1)).nonzero() ,-1)
    #    print '<<< before admore ',np.size(bssss)
        cumul = 0.

        inmask1 = in1mask
    #    print pri
        cumul = addmore(xspe,yspe,x0,y0,x1,y1,pri,cx,cy,xpad,ypad,inmask1)
    ## **** don't try 2nd iteration for now. for some reason this adds lots of extra setup stars! ****
    #    cumul=in1mask
    ##

        ssss = np.reshape( ((cumul == 1) & (pri==-1)).nonzero() ,-1)
    #    print '>>> after admore ',np.size(ssss)

        ok = np.reshape((pri > 0.).nonzero(), -1)
#        print "After second pass: ",np.sum(cumul[ok]).astype('int'),"science slits"\
#            ' with a total weight of %.1f'%np.sum(cumul[ok]*pri[ok])



        totwt = np.sum(cumul * pri)
        if totwt > maxwt: 
            maxwt = totwt
            # -- write results for this realisation
            # only need to write indices of slits in mask:
            output = open('data.pkl', 'wb')
            pickle.dump(cumul,output)
            output.close()

        #      reset
        cumul[:] = 0.
        in1mask[:] = 0.
    #    inwstars[:]=0.

#    print
#    print 'maxwt = %.1f'%maxwt
#    print


    # -- Read in indices of best result
    pkl_file = open('data.pkl', 'rb')
    inds = pickle.load(pkl_file)

#    print np.sum(inds),' slits in mask'

#    drawspectra(xspe,yspe,x0,y0,x1,y1,pri,cx,cy,xpad,ypad,inds,swid,xobj)

    use = np.reshape((inds == 1.0).nonzero(),-1)
##    np.savetxt('maskout.txt',np.transpose((ra[use],dec[use],pri[use])))


    # try adding some shorter slits in to identify where it might be possible to shift slits around

    ##cumul2 = tweakslits(xspe,yspe,x0,y0,x1,y1,pri,cx,cy,xpad,ypad,cumul)

    # **** this is indices of the used flags. prob want running_index[use] as return argument
    print 'these are the returned indexes', use
    print 'Optimizer completed...'
    return use
