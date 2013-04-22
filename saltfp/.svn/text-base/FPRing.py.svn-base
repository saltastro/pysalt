
import math
import numpy as np
from saltfit import Parameter, fit

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
    def __init__(self, xc=0.0, yc=0.0, prad=100.0, norm=1,  sigma=10, ell=0):
        self.xc=xc
        self.yc=yc
        self.norm=norm
        self.prad=prad
        self.sigma=sigma
        self.ell=ell

    def __call__(self, x, y):
        r=((x-self.xc)**2+(y-self.yc)**2)**0.5
        if isinstance(x, np.array):
           return self.norm*np.exp(-((r-self.prad)/self.sigma)**2)
        else:
           return self.norm*math.exp(-((r-self.prad)/self.sigma)**2)


