
import math
import numpy as np
from saltfit import Parameter, fit
from scipy import optimize

class FPRing:
    """This is a general class describing a ring.   As of now, we will define
       a ring to be defined as being described as an ellipse with an Gaussian
       profile around the central line

       A ring is described by the following equation:
       F=A*exp(-(r-r_o)**2/sigma) 
       where:
       A--Peak intensity
       r_o--peak radius
       sigma--width of the line
       And the radius is given by:
       r**2=(x-xc/A)**2+(y-yc/B)**2
 
       For the time being, we will assume that A=B=1

    """
    def __init__(self, xc=0.0, yc=0.0, prad=100.0, norm=1,  sigma=10, ell=0, prad_err=1.0):
        self.xc=xc
        self.yc=yc
        self.norm=norm
        self.prad=prad
        self.prad_err=prad_err
        self.sigma=sigma
        self.ell=ell

    def __call__(self, x, y):
        r=((x-self.xc)**2+(y-self.yc)**2)**0.5
        if isinstance(x, np.ndarray):
           return self.norm*np.exp(-((r-self.prad)/self.sigma)**2)
        else:
           return self.norm*math.exp(-((r-self.prad)/self.sigma)**2)


def ringfit(data, mask=None, fpring=None, full_output=1):
   """Ring fit takes an array of data and then fits an fpring to it.  
      
      Parameters
      ----------
    
      data: ndarray
        Array of values 

      mask: ndarray
        Optional mask for the data

      fpring: FPRing
        Option initial guess for the FPRing

      Output
      ------
   
      ring: FPRing
        Returns an FPRing that is the best fit to the data
   """

   if mask is None: 
      mask=(data>-999)

   if fpring is None:
      ylen,xlen=data.shape
      fpring=FPRing(0.5*xlen, 0.5*ylen, 383)

   #setup the indices
   y,x=np.indices(data.shape)
 
   #setup the parameters

   param=[fpring.xc, fpring.yc, fpring.prad, fpring.norm, fpring.sigma, 0]
   #set up the call
   def call(p, x, y):
       fpr=FPRing(p[0], p[1], p[2], p[3],p[4])
       return fpr(x,y)+p[5]

   #setup the erf
   def erf(p, x, y):
      return data[mask]-call(p,x,y)[mask] #call(x,y).ravel()

   results=optimize.leastsq(erf, param, args=(x,y)) #, full_output=full_output)
   fpring=FPRing(results[0][0], results[0][1],results[0][2],results[0][3],results[0][4])
   return fpring
