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

"""Contains the algorithms to perform photometry on slotmode data."""

# Ensure python 2.5 compatibility
from __future__ import with_statement

from pyraf import iraf
import saltprint, saltsafeio, saltstat
from salterror import SaltError
from saltstat import median_image
import numpy as np

def subbackground(image, sig, nbin, order, iter, type):
    """Calculate and subtract a global background from an image
    using a number of different methods.  This is assuming that
    slotmode data is being used.

    return data array
    """

    # return the same image if no background subtraction is required
    if type=='none': return image

    # calculate the image statistics
    if np.size(image)>0:
        mean, median, stddev=saltstat.iterstat(image,sig,iter)
    if stddev==0:
        raise SaltError('Standard deviation = 0')

    # Create a row median image. This will be useful for a number of purposes
    image_mrow=med_row_image(image)

    # Replace bright objects in the frame with the median value of that row
    mask=((image-mean) < sig*stddev)*(image-mean > -sig*stddev)
    image_sig=image*mask+(1-mask)*image_mrow

    if type == 'median-row':
        image_back=image*0.0+image_mrow

    # Median smooth the image
    if type == 'median':
        image_back=median_image(image_sig,nbin)

    # Make a fit to the surface
    if type=='surf':
        image_back=image.copy()

    if type=='both':
        image_surf=image.copy()
        for x in range(len(image[0])):
            y=np.arange(len(image))
            my=((image[:,x]-mean) < sig*stddev)*(image[:,x] -mean > -sig*stddev)
            cy=np.compress(my,y)
            cim=np.compress(my,image_med[:,x])
            coef=np.polyfit(cy,cim,order)
            image_surf[:,x]=np.polyval(coef,y)

    # subtract the fit
    image=image-image_back

    # return the image
    return image

def med_row_image(arr):
    """Calculate the median value for each row in the array

    returns arr
    """
    nrows=len(arr)

    try:
        marr=np.zeros((nrows,1),dtype=float)
        for i in range(nrows):
            marr[i]=np.median(arr[i])
    except Exception, e:
        raise SaltError('Could not calculate median row values ' + e)

    return marr

