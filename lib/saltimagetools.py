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

"""Module for working with image data."""

import numpy as np

def find_centroid(image):
    """Find the centroid *(cx,cy)* of *image* where:

    .. math::

        cx=\sum f\cdot x/\sum f

        cy=\sum f\cdot y/\sum f

    """

    # Ensure that input is a numpy array
    image=np.asarray(image)

    # Setup indices
    x,y=np.indices(image.shape)

    # Calculate total flux
    tot=image.sum()

    cx=(image*x).sum()/tot
    cy=(image*y).sum()/tot

    return cx,cy

def find_object(image,x,y,distance=5):
    """Returns image pixel coordinates of centroid in box of size::

        2*distance+1

    around coordinates *(x,y)* in *image*.
    """

    # Ensure input image is numpy array
    image=np.asarray(image).transpose()

    # Round to nearest integer for search box
    x=int(round(x))
    y=int(round(y))
    distance=int(round(distance))

    #print image

    # Set range and check for regions outside boundary
    xstart=x-distance
    if xstart<0:
        xstart=0

    ystart=y-distance
    if ystart<0:
        ystart=0

    xend=x+distance+1
    yend=y+distance+1

    #print xstart,xend,ystart,yend

    #print image.shape

    section=image[xstart:xend,ystart:yend]

    #print section
    
    cx,cy=find_centroid(section)

    return cx+xstart,cy+ystart

def zscale(image, contrast=1.0):
    """Implementation of the IRAF zscale algorithm to find vmin and vmax parameters for the dynamic range of a display. It finds the image values near the median image value without the time consuming process of computing a full image histogram."""

    from scipy import optimize
    #import matplotlib.pyplot as plt

    # Get ordered list of points
    I=np.sort(image.flatten())

    # Get number of points
    npoints=len(I)

    # Find the midpoint (median)
    midpoint=(npoints-1)/2

    # Fit a linear function
    # I(i) = intercept + slope * (i - midpoint)

    fitfunc = lambda p, x: p[0]*x+p[1]
    errfunc = lambda p, x, y: fitfunc(p, x) - y

    # Initial guess for the parameters
    p0 = [(I[-1]-I[0])/npoints,I[midpoint]] 

    # Fit
    i=np.arange(len(I))
    p1, success = optimize.leastsq(errfunc, p0[:], args=(i, I))

#    plt.plot(i,I,'r+')
#    plt.plot(i,fitfunc(p1,i))
#    plt.show()

    if success in [1,2,3,4]:
        slope=p1[0]
        z1=I[midpoint]+(slope/contrast)*(1-midpoint)
        z2=I[midpoint]+(slope/contrast)*(npoints-midpoint)
    else:
        z1=np.min(image)
        z2=np.max(image)
    
    return z1, z2
