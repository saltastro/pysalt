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

"""saltstat contains statistical functions"""

import numpy as np

from salterror import SaltError

def mean(list):
    """calculate mean of numeric list"""

    total = 0
    for item in list:
        total += item
    mean = total / len(list)
    return mean


def median(x,logfile=None):
   """calculate median of numeric list

       logfile--depreciated variable
   """
   try:
       return np.median(x)
   except Exception, e:
       message = 'Cannont calculate median because %s' % e
       raise SaltError(message)

def mad(x):
   """Calculated the Median Absolute Deviation defined as:
      MAD=median(|x - median(x)|)
   """
   return np.median(abs(x-np.median(x)))

def median2d(arrays,logfile=None):
    """calculate median of 2d array

       logfile--depreciated variable
    """

    try:
        arrays = arrays.ravel()
        median = np.median(arrays)
    except Exception, e:
        median=None
        message = 'ERROR -- SALTSTAT.MEDIAN2D: Cannot median image arrays because %s' % e
        raise SaltError(message)

    return median

def mean2d(arrays):
    """calculate mean of 2d array"""

    mean = arrays[0]
    for image in arrays[1:]:
        mean += image
    mean /= len(arrays)

    return mean

def std2dclip(arrays, mean, std, sig):
    """calculate clipped std of 2d array"""

    if np.size(arrays)==0: return 0
    mask=(abs(arrays-mean) < sig*std)
    nsize=np.sum(mask)
    if nsize > 0:
       stddev=arrays[mask].std()
    else:
        return 0
    return stddev


def mean2dclip(arrays, mean, std, sig):
    """calculate the sigma clipped mean of 2d array"""

    if np.size(arrays)==0: return 0
    mask=(abs(arrays-mean) < sig*std)
    if np.sum(mask) > 0:
        mean=arrays[mask].mean()
    else:
        return 0
    return mean

def median2dclip(arr, mean, std, sig):
    """calculate the sigma clipped median of 2d array"""
    if np.size(arr)==0: return 0
    try:
        arr = arr.ravel()
        mask=(abs(arr-mean) < sig*std)
        median = np.median(arr[mask])
    except Exception, e:
        median=-1
    return median

def iterstat(arr, sig, niter, verbose=False):
    """iterstas calculates an arrays statistics using
    a sigma clipped values
    """
    mean=arr.mean()
    std=arr.std()
    median=np.median(arr)
    if verbose:  print mean, median, std
    for i in range(niter):
        mask=(abs(arr-mean)<sig*std)
        mean=arr[mask].mean()
        std=arr[mask].std()
        median=np.median(arr[mask])
        if verbose:  print i,mask.sum(), mean, median, std
    return mean, median, std

def median_combine(arrays, logfile=None, axis=0):
    """Median combine a set of arrays
  
      logfile--depreciated variable
    """

    status = 0
    try:
        median = np.median(arrays, axis=axis)
    except Exception, e:
        median=None
        message = 'ERROR -- SALTSTAT.MEDIAN_COMBINE: Cannot median combine arrays because %s' % e
        raise SaltError(message)

    return median, status

def median_image(arr, nbin):
    """Median smooth an image with a filter size set by bin

    returns arr
    """
    from scipy.ndimage.filters import median_filter
    try:
        arr=median_filter(arr,size=(nbin,nbin))
    except Exception, e:
        raise SaltError('Could not median filter image because %s' % e)

    return arr

def median_absolute_deviation(a, axis=None):
    """Compute the median absolute deviation

    Returns the median absolute deviation  of the array elements.  The MAD is
    defined as median(|a-median(a)|).

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    axis : int, optional
        Axis along which the medians are computed. The default (axis=None)
        is to compute the median along a flattened version of the array.

    Returns
    -------
    median_absolute_deviation : ndarray
        A new array holding the result. If the input contains
        integers, or floats of smaller precision than 64, then the output
    Examples
    --------

    This will generate random variates from a Gaussian distribution and return
    the median absolute deviation for that distribution::

        >>> from astropy.stats import median_aboslute_deviation
        >>> from numpy.random import randn
        >>> randvar = randn(10000)
        >>> mad = median_absolute_deviation(randvar)

    See Also
    --------
    median

    """

    a = np.array(a, copy=False)
    a_median = np.median(a, axis=axis)

    #re-broadcast the output median array to subtract it
    if axis is not None:
        shape = list(a_median.shape)
        shape.append(1)
        a_median = a_median.reshape(shape)

    #calculated the median average deviation
    return np.median(np.abs(a - a_median), axis=axis)


