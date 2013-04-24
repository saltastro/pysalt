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

# salttran contains tasks involving the transformation of images

import saltprint
import numpy, pyfits

from salterror import SaltError

# -----------------------------------------------------------
# Embed an array in a larger array

def embed(arr, x0, y0, emarr):
    """Directly embed a 2-D array into a larger array"""
    status=0
    try:
        ylen=len(arr)
        xlen=len(arr[0])
        section='0:%i,0:%i' % (ylen, xlen)
        x1=int(x0)
        x2=int(x1+xlen)
        y1=int(y0)
        y2=int(y1+ylen)
        emarr[y1:y2,x1:x2]=arr
    except Exception, e:
        message = 'ERROR: SALTTRAN.EMBED--Could not embed array because %s ' % e
        raise SaltError(message)

    return emarr

#-----------------------------------------------------------
# shift the data in x-y

def xyshift(arr, xshift, yshift):
    """Shift the data by an integer x and y value"""
    narr=arr
    status=0
    x0=0
    y0=0
    try:
        ylen=len(arr)
        xlen=len(arr[0])
        xlen=xlen+abs(xshift)
        ylen=ylen+abs(yshift)
        if xshift > 0: x0=xshift
        if yshift > 0: y0=yshift
        narr=numpy.zeros((ylen,xlen), numpy.float64)
        narr=embed(arr, x0, y0, narr)
    except:
        message = 'ERROR: SALTTRAN.XYSHIFT -- Could not shift array'
        raise SaltError(message)

    return narr

#-----------------------------------------------------------
# rebin the data

def rebin_factor(a, newshape):
    """Rebin the data to a newshape"""

    status=0
    try:
        assert len(a.shape) == len(newshape)

        slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
        coordinates = numpy.mgrid[slices]
        indices = coordinates.astype('i')   #choose the biggest smaller integer index
        b=a[tuple(indices)]
    except Exception, e:
        b=a
        message='SALTTRAN--ERROR:  Unable to rebin the data because %s' %e
        raise SaltError(message)
    return b

#-----------------------------------------------------------
# block average the data

def blockave (a, newshape):
    """block average the image 
        
       From scipy cookbook
    """
    try:
       shape=a.shape
       lenShape = len(shape)
       factor = numpy.asarray(shape)/numpy.asarray(newshape)
       #copied this from M. Tewes
       evList = ['a.reshape('] + \
             ['newshape[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.mean(%d)'%(i+1) for i in range(lenShape)]
       #print ''.join(evList)
       return eval(''.join(evList))
    except Exception, e:
        b=a
        message='Unable to block average array becase %s ' % e
        raise SaltError(message)
