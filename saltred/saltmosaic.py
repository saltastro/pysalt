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

#!/usr/bin/env python
"""
SALTMOSAIC is a task to apply the CCD geometric corrections to MEF style SALT data.

Author                 Version      Date
-----------------------------------------------
Martin Still (SAAO)    0.1          16 Oct 2006
SM Crawford (SAAO)      0.2           19 Mar 2006

Updates
--------------------
20120319   - Update to new error handling
           - Changed the mosaic to use the whole frame and not trim some data off


"""

import os, time
import numpy, pyfits
from pyraf import iraf

from math import cos, sin, pi
from scipy.ndimage import geometric_transform

import saltsafekey as saltkey
import saltsafeio as saltio
import saltsafestring as saltstring
from saltsafelog import logging, history

from salterror import SaltError

debug=True



# -----------------------------------------------------------
# core routine

def saltmosaic(images,outimages,outpref,geomfile,interp='linear',geotran=True, cleanup=True,clobber=False,logfile=None,verbose=True):

   #Start the logging
   with logging(logfile,debug) as log:

       # Check the input images 
       infiles = saltio.argunpack ('Input',images)

       # create list of output files 
       outfiles=saltio.listparse('Outfile', outimages, outpref,infiles,'')

       #verify that the input and output lists are the same length
       saltio.comparelists(infiles,outfiles,'Input','output')



       # does CCD geometry definition file exist
       geomfilefile = geomfile.strip()
       saltio.fileexists(geomfile)

       gap = 0
       xshift = [0, 0]
       yshift = [0, 0]
       rotation = [0, 0]
       gap, xshift, yshift, rotation=saltio.readccdgeom(geomfile)

       # open each raw image file and apply the transformation to it
       for img, oimg in zip(infiles, outfiles):

           #open the structure
           struct = saltio.openfits(img)

           #create the mosaic
           ostruct=make_mosaic(struct, gap, xshift, yshift, rotation, interp_type=interp, geotran=geotran, cleanup=cleanup, log=log, verbose=verbose)

           #update the header information
           # housekeeping keywords
           fname, hist=history(level=1, wrap=False, exclude=['images', 'outimages', 'outpref'])
           saltkey.housekeeping(ostruct[0], 'SMOSAIC', 'Images have been mosaicked ', hist)

           #write the image out
           saltio.writefits(ostruct,oimg, clobber=clobber)
 
           #close the files
           saltio.closefits(struct)
           saltio.closefits(ostruct)



def make_mosaic(struct, gap, xshift, yshift, rotation, interp_type='linear', boundary='constant', constant=0, geotran=True, cleanup=True, log=None, verbose=False):
   """Given a SALT image struct, combine each of the individual amplifiers and apply the geometric CCD transformations to the image"""

   #get the name of the file
   infile=saltkey.getimagename(struct[0], base=True)
   outpath='./'

   # identify instrument
   instrume,keyprep,keygain,keybias,keyxtalk,keyslot = \
            saltkey.instrumid(struct)

   # how many amplifiers?
   nccds = saltkey.get('NCCDS',struct[0])
   amplifiers = nccds * 2

   # CCD geometry coefficients
   if (instrume == 'RSS' or instrume == 'PFIS'):
        xsh = [0., xshift[0], 0., xshift[1]]
        ysh = [0., yshift[0], 0., yshift[1]]
        rot = [0., rotation[0], 0., rotation[1]]
   elif instrume == 'SALTICAM':
        xsh = [0., xshift[0], 0.]
        ysh = [0., yshift[0], 0.]
        rot = [0., rotation[0], 0]

   # how many extensions?
   nextend = saltkey.get('NEXTEND',struct[0])

   # CCD on-chip binning
   xbin, ybin = saltkey.ccdbin(struct[0])

   # create temporary primary extension
   outstruct = []
   outstruct.append(struct[0])
   # define temporary FITS file store tiled CCDs

   tilefile = saltio.tmpfile(outpath)
   tilefile += 'tile.fits'
   tilehdu = [None] * int(nextend/2+1)
   tilehdu[0] = pyfits.PrimaryHDU()
   tilehdu[0].header = struct[0].header

   if log: log.message('', with_stdout=verbose)


   # iterate over amplifiers, stich them to produce file of CCD images
   for i in range(int(nextend/2)):
       hdu = i * 2 + 1
       amplifier = hdu%amplifiers
       if (amplifier == 0): amplifier = amplifiers

       # read DATASEC keywords
       datasec1 = saltkey.get('DATASEC',struct[hdu])
       datasec2 = saltkey.get('DATASEC',struct[hdu+1])
       xdsec1, ydsec1 = saltstring.secsplit(datasec1)
       xdsec2, ydsec2 = saltstring.secsplit(datasec2)

       # read images
       imdata1 = saltio.readimage(struct,hdu )
       imdata2 = saltio.readimage(struct,hdu+1)

       # tile 2n amplifiers to yield n CCD images
       outdata = numpy.zeros((ydsec1[1]+abs(ysh[i+1]/ybin),xdsec1[1]+xdsec2[1]+abs(xsh[i+1]/xbin)),numpy.float32)

       x1=xdsec1[0]-1
       if x1!=0:
          message='The data in %s have not been trimmed prior to mosaicking.' % infile
          log.error(message) 
       if xsh[i+1]<0: x1+=abs(xsh[i+1]/xbin)
       x2=x1+xdsec1[1]
       y1=ydsec1[0]-1
       if ysh[i+1]<0: y1+=abs(ysh[i+1]/ybin)
       y2=y1+ydsec1[1]
       outdata[y1:y2,x1:x2]=\
           imdata1[ydsec1[0]-1:ydsec1[1],xdsec1[0]-1:xdsec1[1]]

       x1=x2
       x2=x1+xdsec2[1]
       y1=ydsec2[0]-1
       if ysh[i+1]<0: y1+=abs(ysh[i+1]/ybin)
       y2=y1+ydsec2[1]
       outdata[y1:y2,x1:x2]=\
           imdata2[ydsec1[0]-1:ydsec1[1],xdsec1[0]-1:xdsec1[1]]

       # size of new image
       naxis1 = str(xdsec1[1] + xdsec2[1])
       naxis2 = str(ydsec1[1])

       # add image and keywords to HDU list
       tilehdu[i+1] = pyfits.ImageHDU(outdata)
       tilehdu[i+1].header = struct[hdu].header
       tilehdu[i+1].header['DATASEC'] = '[1:' + naxis1 + ',1:' + naxis2 + ']'

       # image tile log message #1
       if log:
           message =  os.path.basename(infile) + '[' + str(hdu) + ']['
           message += str(xdsec1[0]) + ':' + str(xdsec1[1]) + ','
           message += str(ydsec1[0]) + ':' + str(ydsec1[1]) + '] --> '
           message += os.path.basename(tilefile) + '[' + str(i+1) + ']['
           message += str(xdsec1[0]) + ':' + str(xdsec1[1]) + ','
           message += str(ydsec1[0]) + ':' + str(ydsec1[1]) + ']'
           log.message(message,with_stdout=verbose, with_header=False)
           message =  os.path.basename(infile) + '[' + str(hdu+1) + ']['
           message += str(xdsec1[0]) + ':' + str(xdsec1[1]) + ','
           message += str(ydsec1[0]) + ':' + str(ydsec1[1]) + '] --> '
           message += os.path.basename(tilefile) + '[' + str(i+1) + ']['
           message += str(xdsec1[1]+1) + ':' + str(xdsec1[1]+xdsec2[1]) + ','
           message += str(ydsec2[0]) + ':' + str(ydsec2[1]) + ']'
           log.message(message,with_stdout=verbose, with_header=False)

   # write temporary file of tiled CCDs
   hdulist = pyfits.HDUList(tilehdu)
   hdulist.writeto(tilefile)

   # iterate over CCDs, transform and rotate images
   yrot = [None] * 4
   xrot = [None] * 4

   tranfile = [' ']
   tranhdu = [0]

   #this is hardwired for SALT where the second CCD is considered the fiducial
   for hdu in range(1,int(nextend/2+1)):
       tranfile.append(' ')
       tranhdu.append(0)
       tranfile[hdu] = saltio.tmpfile(outpath)
       tranfile[hdu] += 'tran.fits'
       ccd = hdu%nccds
       if (ccd == 0): ccd = nccds

       # correct rotation for CCD binning
       yrot[ccd] = rot[ccd] * ybin / xbin
       xrot[ccd] = rot[ccd] * xbin / ybin
       dxshift = xbin * int(float(int(gap) / xbin) + 0.5) - gap

       # transformation using geotran IRAF task
       #if (ccd == 1):
       if (ccd != 2):

           if geotran:
               message  = '\nSALTMOSAIC -- geotran ' + tilefile + '[' + str(ccd) + '] ' + tranfile[hdu]
               message += ' \"\" \"\" xshift=' + str((xsh[ccd]+(2-ccd)*dxshift)/xbin) + ' '
               message += 'yshift=' + str(ysh[ccd]/ybin) + ' xrotation=' + str(xrot[ccd]) + ' '
               message += 'yrotation=' + str(yrot[ccd]) + ' xmag=1 ymag=1 xmin=\'INDEF\''
               message += 'xmax=\'INDEF\' ymin=\'INDEF\' ymax=\'INDEF\' ncols=\'INDEF\' '
               message += 'nlines=\'INDEF\' verbose=\'no\' fluxconserve=\'yes\' nxblock=2048 '
               message += 'nyblock=2048 interpolant=\'' + interp_type + '\' boundary=\'constant\' constant=0'
               log.message(message, with_stdout=verbose)

               yd,xd=tilehdu[ccd].data.shape
               ncols='INDEF' #ncols=xd+abs(xsh[ccd]/xbin)
               nlines='INDEF' #nlines=yd+abs(ysh[ccd]/ybin)
    
               iraf.images.immatch.geotran(tilefile+"["+str(ccd)+"]",
                                                tranfile[hdu],
                                                "",
                                                "",
                                                xshift=(xsh[ccd]+(2-ccd)*dxshift)/xbin,
                                                yshift=ysh[ccd]/ybin,
                                                xrotation=xrot[ccd],
                                                yrotation=yrot[ccd],
                                                xmag=1,ymag=1,xmin='INDEF',xmax='INDEF',ymin='INDEF',
                                                ymax='INDEF',ncols=ncols,nlines=nlines,verbose='no',
                                                fluxconserve='yes',nxblock=2048,nyblock=2048,
                                                interpolant="linear",boundary="constant",constant=0)
               #open the file and copy the data to tranhdu
               tstruct=pyfits.open(tranfile[hdu])
               tranhdu[hdu]=tstruct[0].data
               tstruct.close()

           else:
               log.message("Transform CCD #%i using dx=%s, dy=%s, rot=%s" % (ccd, xsh[ccd]/2.0, ysh[ccd]/2.0, xrot[ccd]), with_stdout=verbose, with_header=False)
               tranhdu[hdu]=geometric_transform(tilehdu[ccd].data, tran_func, prefilter=False, order=1, extra_arguments=(xsh[ccd]/2, ysh[ccd]/2, 1, 1, xrot[ccd], yrot[ccd]))
               tstruct=pyfits.PrimaryHDU(tranhdu[hdu])
               tstruct.writeto(tranfile[hdu])
               
               
       else:
           log.message("Transform CCD #%i using dx=%s, dy=%s, rot=%s" % (ccd, 0, 0, 0), with_stdout=verbose, with_header=False)
           tranhdu[hdu] = tilehdu[ccd].data

   # open outfile
   outlist = 2*[None] 
   #outstruct[0] = pyfits.PrimaryHDU()
   outlist[0] = struct[0].copy()
   naxis1 = int( gap / xbin * (nccds - 1))
   naxis2 = 0
   for i in range(1,nccds+1):
       yw,xw=tranhdu[i].shape
       naxis1 += xw + int(abs(xsh[ccd]/xbin))+1
       naxis2 = max(naxis2, yw)
   outdata = numpy.zeros((naxis2,naxis1),numpy.float32)
   outdata.shape = naxis2, naxis1

   # iterate over CCDs, stich them to produce a full image
   hdu = 0
   totxshift=0
   for hdu in range(1,nccds+1):

       # read DATASEC keywords
       ydsec, xdsec = tranhdu[hdu].shape

       # define size and shape of final image
       # tile CCDs to yield mosaiced image
       x1 = int((hdu-1)*(xdsec+gap/xbin))+int(totxshift)
       x2 = xdsec+x1
       y1 = int(0)
       y2 = int(ydsec)
       outdata[y1:y2,x1:x2] = tranhdu[hdu]
       totxshift+=int(abs(xsh[hdu]/xbin))+1

   #add to the file
   outlist[1] = pyfits.ImageHDU(outdata)

   #create the image structure
   outstruct = pyfits.HDUList(outlist)

   #update the head informaation
   # housekeeping keywords
   saltkey.put('NEXTEND',2, outstruct[0])
   saltkey.new('EXTNAME','SCI','Extension name', outstruct[1])
   saltkey.new('EXTVER',1,'Extension number',outstruct[1])

   try:
      saltkey.copy(struct[1], outstruct[1], 'CCDSUM')
   except:
      pass

   #Add keywords associated with geometry
   gstr='%i %f %f %f %f %f %f' % (gap, xshift[0], yshift[0], rotation[0], xshift[1], yshift[1], rotation[1])
   saltkey.new('SALTGEOM', gstr, 'SALT geometry coefficients', outstruct[0])

   # WCS keywords
   saltkey.new('CRPIX1',0,'WCS: X reference pixel', outstruct[1])
   saltkey.new('CRPIX2',0,'WCS: Y reference pixel', outstruct[1])
   saltkey.new('CRVAL1',float(xbin), 'WCS: X reference coordinate value', outstruct[1])
   saltkey.new('CRVAL2',float(ybin), 'WCS: Y reference coordinate value', outstruct[1])
   saltkey.new('CDELT1',float(xbin),'WCS: X pixel size', outstruct[1])
   saltkey.new('CDELT2',float(ybin),'WCS: Y pixel size', outstruct[1])
   saltkey.new('CTYPE1','pixel','X type', outstruct[1])
   saltkey.new('CTYPE2','pixel','Y type', outstruct[1])

   # cleanup temporary files
   if cleanup:
       for tfile in tranfile:
           if os.path.isfile(tfile):
               saltio.delete(tfile)
       if os.path.isfile(tilefile):
           status = saltio.delete(tilefile)
 

   #return the file
   return outstruct

def tran_func(a, xshift, yshift, xmag, ymag, xrot, yrot):
    return (ymag*a[0]*cos(yrot*pi/180.0)-xmag*a[1]*sin(xrot*pi/180)-yshift, ymag*a[0]*sin(yrot*pi/180.0)+xmag*a[1]*cos(xrot*pi/180)-xshift)


# -----------------------------------------------------------
# main code

parfile = iraf.osfn("saltred$saltmosaic.par")
t = iraf.IrafTaskFactory(taskname="saltmosaic",value=parfile,function=saltmosaic, pkgname='saltred')
