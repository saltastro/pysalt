"""A general module with tools for use with the saltfp package"""
import math
import numpy as np
import scipy.ndimage as nd

from saltfit import interfit
from salterror import SaltError, SaltIOError
from FPRing import FPRing, ringfit


def fpfunc(z, r, t, coef=None):
   """A functional form fitting the Fabry Perot parameterization.  The 
      FP parameterization is given by:

      $\lambda = \\frac{A+Bz+Cz^2+Dz^3}{(1+r/F)^{0.5}}+Et$

      Parameters
      ----------
      z: float or ndarray
         z position of the etalon
      r: float or ndarray
         r position in the image
      t: float or ndarray
         r position in the image
      coef: list or ndarray
         An array with the coefficients for the FP equations.  The
         coefficients should be given as [A B C D E F] 

      Returns
      -------
      w:  float or ndarray
         wavelength at position r
   """
   if len(coef)!=6: raise Exception('Not enough FP Coefficients')
   if coef[5]==0:  raise Exception('F must not be zero')

   #calcuate the value according to the above equation
   w=coef[0]+coef[1]*z+coef[2]*z**2+coef[3]*z**3+coef[4]*t
   w=w/(1+(r/coef[5])**2)**0.5
   return w

def findrings(data, thresh=5, niter=5, minsize=10, axc=None, ayc=None):
    """findrings makes a rough calculation for the parameters of the rings
       based on single line cuts through the data.  It returns a list of rings
    """
    ring_list=[]

    #first guess the middle is in the middle of the data
    if axc is None:   
        xc=int(0.5*len(data[0]))
    else:
        xc=axc
   
    if ayc is None: 
        yc=int(0.5*len(data))
    else: 
        yc=ayc

    #take a look at the y cut through the data
    xdata=data[yc,:]
    #take a look through the xdata.  check for the same thing and make sure they are consistent
    ydata=data[:,xc]

    #get rid of all the lower points
    #find the peaks in the data

    ypeak_list=findpeaks(ydata, 0.4, minsize)
    xpeak_list=findpeaks(xdata, 0.4, minsize)

    if abs(len(ypeak_list)-len(xpeak_list))>1:
       msg="Non-symmetrically rings in the image"
       #raise SaltError(msg)

    nrings=int(max(len(ypeak_list)/2, len(xpeak_list)/2))

    #throw an error if no rings are detected
    if nrings<1:
       msg="No rings detected in image"
       raise SaltError(msg)
    #loop through the image and determine parameters of rings
    for i in range(0,nrings,2):
       #determine the y-center
       try:
           y1,y2=ypeak_list[i]
           yarr=np.arange(y1,y2)
           ypa=y1+ydata[y1:y2].argmax()
           ysiga=(abs(np.sum((yarr-ypa)**2*ydata[y1:y2])/ydata[y1:y2].sum()))**0.5

           y1,y2=ypeak_list[i+1]
           yarr=np.arange(y1,y2)
           ypb=y1+ydata[y1:y2].argmax()
           ysigb=(abs(np.sum((yarr-ypb)**2*ydata[y1:y2])/ydata[y1:y2].sum()))**0.5

           if ayc is None:  
               yc=0.5*(ypa+ypb)
           else:
              yc=ayc
           ymax=max(ydata[ypa], ydata[ypb])
           yrad=0.5*abs(ypb-ypa)
           ysig=0.5*(ysiga+ysigb)
       except Exception, e:
           yc=yc
           yrad=0
           ysig=0
           ymax=ydata.max()

       #determine the x-center
       try:
           x1,x2=xpeak_list[i]
           xarr=np.arange(x1,x2)
           xpa=x1+xdata[x1:x2].argmax()
           xsiga=(abs(np.sum((xarr-xpa)**2*xdata[x1:x2])/xdata[x1:x2].sum()))**0.5

           x1,x2=xpeak_list[i+1]
           xpb=x1+xdata[x1:x2].argmax()
           xarr=np.arange(x1,x2)
           xsigb=(abs(np.sum((xarr-xpb)**2*xdata[x1:x2])/xdata[x1:x2].sum()))**0.5

           if axc is None:  
               xc=0.5*(xpa+xpb)
           else:
               xc=axc
           xmax=max(xdata[xpa], xdata[xpb])
           xsig=0.5*(xsiga+xsigb)
           xrad=0.5*abs(xpa-xpb)
       except:
           xc=xc
           xrad=0
           xsig=0
           xmax=xdata.max()
       prad_err=max(1.0, 0.5*abs(yrad-xrad))
       ring_list.append(FPRing(xc, yc, max(yrad,xrad), max(xmax,ymax), max(xsig,ysig), prad_err=prad_err))

    return ring_list


def findcenter(data, ring, method, niter=5, conv=0.05):
   method=method.upper()
   if method == 'FIT':
      ring=ringfit(data, fpring=ring)
   elif method == 'MAX':
      i=0
      c=conv+1
      while i < niter and c>conv:
         xc,yc=maxflux_center(data, ring.xc, ring.yc, ring.prad, 10, maxiter=20)
         c=((ring.xc-xc)**2+(ring.yc-yc)**2)**0.5
         i+=1
         rad, rad_err=findradius(data, xc, yc, ring.prad, 10)
         ring.xc=xc
         ring.yc=yc
         ring.prad=rad
         ring.prad_err=rad_err
   elif method == 'CENTER':
         xc,yc,rad, rad_err=centerring(data, ring.xc, ring.yc, radmax=ring.prad, radstep=ring.sigma, nbins=8)
         c=((ring.xc-xc)**2+(ring.yc-yc)**2)**0.5
         ring.xc=xc
         ring.yc=yc
         ring.prad=rad
         ring.prad_err=rad_err
   elif method == 'MOMENT':
      pass
   else:
      raise SaltError('%s is not a valid method' % method)
   return ring

def maxflux_center(data, axc=None, ayc=None,radmax=450, radstep=5, maxiter=100):
   """Find the center of the data by trying to maximize the flux in the radial distribution
   """

   ylen,xlen=data.shape
   if axc is None: axc=0.5*xlen
   if ayc is None: ayc=0.5*ylen

   mflux=0
   niter=0
   found=True
   bxc=axc
   byc=ayc
   while found and niter<maxiter:
    niter+=1
    found=False
    for i in [-1,0,1]:
     axc=bxc+i
     for j in [-1,0,1]:
      ayc=byc+j
      flux=calcflux(data, axc, ayc, radmax, radstep)
      if mflux< flux:
         bxc=axc
         byc=ayc
         mflux=flux
         found=True
         #print bxc, byc, flux
         continue
   return bxc, byc

def findradius(data, axc=None, ayc=None,radmax=450, radstep=5, maxiter=100, rstep=0.25):
   """Find the radius of the ring by trying to maximum the value"""
   ylen,xlen=data.shape
   if axc is None: axc=0.5*xlen
   if ayc is None: ayc=0.5*ylen

   mflux=calcflux(data, axc, ayc, radmax, radstep)
   niter=0
   found=True
   brad=radmax
   brad_err=max(brad*(mflux/mflux**2)**0.5,1.0)
   while found and niter<100:
    niter+=1
    found=False
    for i in [-rstep,rstep]:
      rad=brad+i
      flux=calcflux(data, axc, ayc, rad, radstep)
      if mflux< flux:
         brad=rad
         brad_err=(rad*(flux/flux**2)**0.5)
         mflux=flux
         found=True
         continue
   return brad, brad_err

def calcflux(data, axc, ayc, radmax, radstep):
   """Calculate the total flux between radmax-radstep and radmax+radstep"""
   y,x=np.indices(data.shape)
   r=((x-axc)**2+(y-ayc)**2)**0.5
   mask=(r>radmax-radstep)*(r<radmax+radstep)
   return data[mask].sum()


def centerring(data, axc, ayc, radmax=450, radstep=50, nbins=8):
   """Calculate the center of the ring by determining the radius of the 
      line in several bins
   """
   rad=radmax
   rad_err=1.0

   #set up the radius and theta
   y,x=np.indices(data.shape)
   r=((x-axc)**2+(y-ayc)**2)**0.5
   theta=np.arctan((y-ayc)/(x-axc))
   theta[(x-axc<0)]+=math.pi
   theta += 0.5*math.pi

   #calculate the centroid in each bin
   nsteps=2*math.pi/nbins
   rad_arr=np.zeros(nbins)
   theta_arr=np.zeros(nbins)
   for i in range(nbins):
       t1=i*nsteps
       t2=t1+nsteps
       mask=(theta>t1)*(theta<t2)*(abs(r-radmax)<radstep)
       theta_arr[i]=0.5*(t1+t2)
       try:
           rad_arr[i]=fitradius(r[mask], data[mask])
       except:
           pass
       #r[mask][j]
   x_arr=rad_arr*np.cos(theta_arr-0.5*math.pi)
   y_arr=rad_arr*np.sin(theta_arr-0.5*math.pi)
   return axc+x_arr.mean(), ayc+y_arr.mean(), rad_arr.mean(), rad_arr.std()

def fitradius(radius, data):
    """Fit a line to the data"""
    it=interfit(radius, data, function='polynomial', order=3)
    it.interfit()
    d=it(radius)
    return radius[d.argmax()]


def findpeaks(data, fpeak=0.8,minsize=10):
    """Find peakes median filters an image and finds the peaks in the image
    """
    #median filter the image
    mdata=nd.filters.median_filter(data, size=minsize)

    #take the 80% points and find the peaks
    mask=(mdata>fpeak*mdata.max())

    #find all the objects
    obj_arr, obj_num=nd.label(mask)

    #ypeaks
    peaks=[]
    for i in range(obj_num):
        pid=np.where(obj_arr==i+1)[0]
        if len(pid)>=minsize:
           peaks.append((pid.min(), pid.max()))

    return peaks
   

