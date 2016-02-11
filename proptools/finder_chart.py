#!/opt/local/bin/python

import cStringIO
import os
import sys
import xml
import base64
import urllib2
import xml.dom.minidom 
import pyfits
import ephem
import numpy as np
import aplpy


# grab MOS xml definition from WM given account and barcode
def get_slitmask_xml(username, password, barcode):
    """
    Return the slit mask XML as a DOM document.
    """
    encoded_username = base64.encodestring(username).strip()
    encoded_password = base64.encodestring(password).strip()

    mask_url = 'https://www.salt.ac.za/wm/downloads/SlitmaskXml.php'

    # We pass the parameters in a GET request.
    url = mask_url + '?username=%s&password=%s&Barcode=%s' % (encoded_username, encoded_password, barcode)
    response = urllib2.urlopen(url)
    dom = xml.dom.minidom.parse(response)

    # Handle the case that the request wasn't successful.
    if dom.documentElement.tagName == 'Invalid':
        raise Exception('You are not allowed to view the slit mask XML.')

    return dom

# grab 10' x 10' image from STScI server and pull it into pyfits
def get_dss(imserver, ra, dec):
    url = "http://archive.stsci.edu/cgi-bin/dss_search?v=%s&r=%f&d=%f&e=J2000&h=10.0&w=10.0&f=fits&c=none" % (imserver, ra, dec)
    print url
    fitsData = cStringIO.StringIO()
    data = urllib2.urlopen(url).read()
    fitsData.write(data)
    fitsData.seek(0)
    return pyfits.open(fitsData)

# grab uploaded base64-encoded FITS
def get_fits(b64str):
    fitsData = cStringIO.StringIO()
    fitsData.write(base64.b64decode(b64str))
    fitsData.seek(0)
    return pyfits.open(fitsData)

# draw a line centered at ra,dec of a given length at a given angle
def draw_line(plot, theta, length, ra, dec, color='b', linewidth=1, alpha=0.7):
    theta = theta*np.pi/180.0
    length = length/2.0
    dx = np.sin(theta)*length/(np.cos(dec*np.pi/180.0)*60.0)
    dy = np.cos(theta)*length/60.0
    coords = np.array([[ra+dx, ra-dx], [dec+dy, dec-dy]])
    plot.show_lines([coords], color=color, linewidth=linewidth, alpha=alpha)
    return plot

def read_mos_xml(xmlfile):
    mask_xml = open(xmlfile).read()
    doc = xml.dom.minidom.parseString(mask_xml)
    barcode=''

    pars = doc.getElementsByTagName("parameter")
    slits = doc.getElementsByTagName("slit")
    refs = doc.getElementsByTagName("refstar")

    parameters = {}

    for par in pars:
        name = par.getAttribute("name")
        val = par.getAttribute("value")
        parameters[name] = val

    ra = float(parameters["CENTERRA"])
    dec = float(parameters["CENTERDEC"])
    pa = float(parameters["ROTANGLE"])
    propcode = parameters["PROPOSALCODE"]
    pi_email = parameters["PI"]
    barcode = "Mask #%s" % parameters["MASKNUM"]
    title = "%s (%s; %s)" % (barcode, propcode, pi_email)
    return ra, dec, pa, propcode, pi_email, barcode, title,slits, refs

# draw slits and reference boxes for MOS
def mos_plot(plot, slits, refs):
    # draw the slits
    slit_ra = []
    slit_dec = []
    slit_width = []
    slit_length = []
    for slit in slits:
        slit_ra.append(float(slit.attributes['xce'].value))
        slit_dec.append(float(slit.attributes['yce'].value))
        slit_width.append(float(slit.attributes['width'].value)/3600.0)
        slit_length.append(float(slit.attributes['length'].value)/3600.0)
      
    plot.show_rectangles(slit_ra, slit_dec, slit_width, slit_length, edgecolor='r', alpha=0.7)
    # make bigger boxes around the reference objects
    ref_ra = []
    ref_dec = []
    ref_width = []
    ref_height = []

    for ref in refs:
        ref_ra.append(float(ref.attributes['xce'].value))
        ref_dec.append(float(ref.attributes['yce'].value))
        ref_width.append(5.0/3600.0)
        ref_height.append(5.0/3600.0)

    if len(ref_ra) > 0:
        plot.show_rectangles(ref_ra, ref_dec, ref_width, ref_height, edgecolor=(1,1,0), linewidth=2, alpha=0.7)
    return plot

# set up basic plot
def init_plot(hdu, imserver, title, ra, dec):
    servname = {}
    servname['none']=''
    servname['poss2ukstu_red'] = "POSS2/UKSTU Red"
    servname['poss2ukstu_blue'] = "POSS2/UKSTU Blue"
    servname['poss2ukstu_ir'] = "POSS2/UKSTU IR"
    servname['poss1_blue'] = "POSS1 Blue"
    servname['poss1_red'] = "POSS1 Red"

    out = sys.stdout
    sys.stdout = open("/dev/null", 'w')
    plot = aplpy.FITSFigure(hdu)
    plot.show_grayscale()
    plot.set_theme('publication')
    sys.stdout = out
    plot.add_label(0.5, 1.03,
                  title,
                  relative=True, style='italic', weight='bold', size='large')
    plot.add_label(-0.05, -0.05, "%s" % servname[imserver], relative=True, style='italic', weight='bold')

    plot.add_grid()
    plot.grid.set_alpha(0.2)
    plot.grid.set_color('b')

    plot.show_circles([ra, ra], [dec, dec], [4.0/60.0, 5.0/60.0], edgecolor='g')
    plot.add_label(0.79,
                    0.79,
                    "RSS",
                    relative=True,
                    style='italic',
                    weight='bold',
                    size='large',
                    horizontalalignment='left',
                    color=(0,0,1))
    plot.add_label(0.86,
                  0.86,
                  "SCAM",
                  relative=True,
                  style='italic',
                  weight='bold',
                  size='large',
                  horizontalalignment='left',
                  color=(0,0,1))
    plot.add_label(ra,
                  dec+4.8/60.0,
                  "N",
                  style='italic',
                  weight='bold',
                  size='large',
                  color=(0,0.5,1))
    plot.add_label(ra+4.8/(np.abs(np.cos(dec*np.pi/180.0))*60),
                  dec,
                  "E",
                  style='italic',
                  weight='bold',
                  size='large',
                  horizontalalignment='right',
                  color=(0,0.5,1))
    plot = draw_line(plot, 0, 8, ra, dec, color='g', linewidth=0.5, alpha=1.0)
    plot = draw_line(plot, 90, 8, ra, dec, color='g', linewidth=0.5, alpha=1.0)
    return plot

def finderchart(xmlfile, image=None, outfile=None):
    """Given an image and an xml file, create the finder chart"""
    
    #read in the xml
    ra, dec, pa, propcode, pi_email, barcode, title, slits, refs=read_mos_xml(xmlfile)

    #read in the image
    if image:
        hdu = pyfits.open(image)
        imserver='none'
    else:
        imserver='poss2ukstu_red'
        hdu = get_dss(imserver, ra, dec)
 
    #create the plot
    plot = init_plot(hdu, imserver, title, ra, dec)

    #add the slits
    plot = mos_plot(plot, slits, refs)
   
    #save the plot
    if outfile:
        plot.save(str(outfile))
   
 

if __name__=='__main__':  
   finderchart(sys.argv[1], sys.argv[2], sys.argv[3])
