################################# LICENSE ##################################
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.  See LICENSE file for more details                  #
#                                                                          #
############################################################################


#!/usr/bin/env python

"""
QUICKSPEC -- QUICKSPEC  provides a plugin for saltfirst that provides quick 
reductions for spectroscopic data

Author                 Version      Date
-----------------------------------------------
S M Crawford (SAAO)    0.1           5 Jun 2011

"""
import os 
import pyfits
import numpy as np
import matplotlib
matplotlib.use('GTkAgg')

#pysalt imports
from pyraf import iraf
from pyraf.iraf import pysalt

from scipy.ndimage.filters import median_filter

from specrectify import specrectify
from specextract import extract, write_extract
from specsky import specsky
from specslitnormalize import create_response 


from PySpectrograph.Models import RSSModel

import spectools as st
import saltsafekey as saltkey
import saltsafeio

import pylab as pl 

sky_lines = [4358.34, 5577.338, 6300.304,6363.8000,7715.0116]

def quickspec(profile, lampid=None, findobj=False, objsection=None, skysection=None, clobber=False, logfile='saltclean.log', verbose=True):
   """From mosaicked data, produce wavelength calibrated files"""
   profile = os.path.basename(profile)

   #fill in the mosaic 
   fillgaps(profile)

   #specrectify
   specrectify(profile, outimages='', outpref='s', solfile=None, caltype='rss',
                   function='legendre',  order=3, inttype='interp', w1=None, w2=None, dw=None, nw=None,
                   blank=0.0, clobber=True, logfile=logfile, verbose=True)
   y1,y2=quickap('s' + profile, lampid, findobj, objsection, skysection, clobber, logfile, verbose)

   if skysection is None:
      ylen=100
      skysection='[%i:%i]' % (y2+0.1*ylen,y2+0.2*ylen)
   
   specsky('s'+profile, outimages='s'+profile, outpref='', method='normal', section=skysection, function='polynomial', order=3, clobber=clobber, logfile=logfile, verbose=verbose)

   return y1,y2

def fillgaps(profile):
   """Fill in the gaps"""
   hdu = pyfits.open(profile)

   #get the binning
   xbin, ybin = saltkey.ccdbin( hdu[0], '')

   #fill in the second gap
   x1 = 2044/xbin
   x2 = 2152/xbin
   d1=np.min(hdu[1].data[:,x1-10:x1], axis=1)
   d2=np.min(hdu[1].data[:,x2:x2+10], axis=1)
   d1 = 0.5 * (d1+d2)
   hdu[1].data[:,x1:x2] = d1.reshape(len(d1),1)


   #fill in the second gap
   x1 = 4096/xbin
   x2 = 4308/xbin
   d1=np.min(hdu[1].data[:,x1-10:x1], axis=1)
   d2=np.min(hdu[1].data[:,x2:x2+10], axis=1)
   d1 = 0.5 * (d1+d2)
   hdu[1].data[:,x1:x2] = d1.reshape(len(d1),1)
   
   if os.path.isfile(profile): os.remove(profile)
   hdu.writeto(profile)


def quickap(profile, lampid=None, findobj=False, objsection=None, skysection=None, clobber=False, logfile='saltclean.log', verbose=True):
   profile = os.path.basename(profile)
   outfile = profile.replace('fits', 'txt')
   #open the file
   struct=pyfits.open(profile)
   data=struct[1].data
   xlen=len(struct[1].data[0])
   ylen=len(struct[1].data)

   #determine the resolution element
   dres = 1.0*calc_resolution(struct)

   #correct for the response function
   response = create_response(data, spatial_axis=1, order=3, conv=1e-2, niter=10)
   data = data / response


   if objsection is None:
      y1=0.5*ylen-0.05*ylen
      y2=0.5*ylen+0.05*ylen
      objsection='[%i:%i]' % (y1,y2)
   else:
      y1=int(objsection[1:].split(':')[0])
      y2=int(objsection[:-1].split(':')[1])
 

   #extract object
   minsize=5
   thresh=5
   ap_spec=extract(struct, method='normal', section=[(y1,y2)], minsize=minsize, thresh=thresh, convert=True)[0]

   #if it is a lamp, do not sky subtract it
   if saltsafeio.checkfornone(lampid):
        write_extract(outfile, [ap_spec], outformat='ascii', clobber=clobber)
        return y1,y2

   if skysection is None:
      skysection='%i:%i' % (y2+0.1*ylen,y2+0.2*ylen)
   #sky spectrum
   sy1, sy2 = skysection.split(':')
   sy1=int(sy1)
   sy2=int(sy2)

   sk_spec=extract(struct, method='normal', section=[(sy1,sy2)], minsize=minsize, thresh=thresh, convert=True)[0] 
   sk_spec.ldata - sk_spec.ldata
   
   w1=ap_spec.wave.min()
   w2=ap_spec.wave.max()
   nsect=10
   dw=(w2-w1)/nsect
   for w in np.arange(w1,w2,dw):
      mask=(ap_spec.wave>w)*(ap_spec.wave<w+dw)
      sk_spec=quickcross(ap_spec, sk_spec, mask=mask, dx=0.02, nstep=500)
      ap_spec.ldata[mask]=ap_spec.ldata[mask]-float(y2-y1)/(sy2-sy1)*sk_spec.ldata[mask]

   #set up masks for sky subtraction
   amask = (ap_spec.ldata==0)
   smask = (sk_spec.ldata==0)


   #find offset between two spectra
   d1 = ap_spec.ldata
   d2 = sk_spec.ldata*float(y2-y1)/(sy2-sy1) 
   y=d1[d2>0]/d2[d2>0]
   y=np.median(y)

   #subtract the skyspectrum

   clean_spectra(ap_spec, dres=2*dres, grow=dres)
   #pl.plot(ap_spec.wave, ap_spec.ldata, ls='', marker='o', ms=2)
   #pl.show()
   #write it out and return
   write_extract(outfile, [ap_spec], outformat='ascii', clobber=clobber)

   #clean up the data 
   if clobber: 
      print 'Replacing pixels'
      os.remove(profile)
      struct[1].data[data==0]=np.median(data)
      struct.writeto(profile)

   return  y1,y2

def clean_spectra(ap_spec, dres=2, grow=0):
    """Clean a spectrum and make it more presenatable"""
    xmax=len(ap_spec.ldata)
    for w in sky_lines:
        mask=(abs(ap_spec.wave-w)<dres)
        ap_spec.ldata[mask] = 0

    ap_spec.ldata[0] = 0
    ap_spec.ldata[-1] = 0


    #grow the spectra
    if grow:
      nid = np.where(ap_spec.ldata==0)[0]
      for i in nid:
          x1=max(0,int(i-grow))
          x2=min(int(i+grow),xmax-1)
          ap_spec.ldata[x1:x2]=0
 
    #remove the bad pixels 
    mask = (ap_spec.ldata!=0)
    ap_spec.ldata=ap_spec.ldata[mask]
    ap_spec.wave=ap_spec.wave[mask]

    return ap_spec

def quickcross(spec1, spec2, mask=None, dx=1, nstep=20):
    "Calculate teh offset and return matched spectrum"
    if mask is None: mask=(spec1.wave>0)
    cc_arr=np.zeros(nstep)
    start_x = - 0.5*nstep*dx
    for i in range(nstep):
        x=start_x+i*dx
        data= np.interp(spec1.wave[mask], spec2.wave[mask]+x, spec2.ldata[mask])
        cc_arr[i]=st.ncor(data, spec1.ldata[mask])
    x=start_x+cc_arr.argmax()*dx
    spec2.wave[mask] = spec2.wave[mask]+x #np.interp(spec1.wave, spec2.wave+x, spec2.ldata)
    return spec2

def calc_resolution(hdu):
    """Calculate the resolution for a setup"""
    instrume=saltkey.get('INSTRUME', hdu[0]).strip()
    grating=saltkey.get('GRATING', hdu[0]).strip()
    grang=saltkey.get('GR-ANGLE', hdu[0])
    grasteps=saltkey.get('GRTILT', hdu[0])
    arang=saltkey.get('AR-ANGLE', hdu[0])
    arsteps=saltkey.get('CAMANG', hdu[0])
    rssfilter=saltkey.get('FILTER', hdu[0])
    specmode=saltkey.get('OBSMODE', hdu[0])
    masktype=saltkey.get('MASKTYP', hdu[0]).strip().upper()
    slitname=saltkey.get('MASKID', hdu[0])
    xbin, ybin = saltkey.ccdbin( hdu[0], '')
    slit=st.getslitsize(slitname)
    
    #create RSS Model
    rss=RSSModel.RSSModel(grating_name=grating.strip(), gratang=grang, \
                          camang=arang,slit=slit, xbin=xbin, ybin=ybin, \
                          )
    res=1e7*rss.calc_resolelement(rss.alpha(), -rss.beta())
    return res 

    

