#!/usr/bin/env python
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.                                                     #

"""
MOSTOOLS is a package written that handles the edge detection,

Author                 Version      Date
-----------------------------------------------
J. P. Kotze (SAAO)     0.2       19 Nov 2010

TODO
----
- Add error checking
- Add a check to see if the amount of edges are the same.
- Finish the function comments
- take out all the hard coding of numbers.



LIMITATIONS
-----------
1. Need to add checking that the edges correspond

"""

from astropy.io import fits 
from astropy.io.fits import Column

import pywcs

import math
import numpy
import matplotlib.pyplot as pl
from matplotlib.ticker import NullFormatter
import scipy.ndimage.filters
import scipy.interpolate as si
import scipy.signal
from spectools import SALTSpecError
import saltimagetools
from xml.dom import minidom

from slitmask import SlitMask

from astropy import modeling as mod

nullfmt = NullFormatter()

MAX_REJECT = 0.5
MIN_NPIXELS = 5
GOOD_PIXEL = 0
BAD_PIXEL = 1
KREJ = 2.5
MAX_ITERATIONS = 5

# given an image flatten it for a certain x range for all the y values
# for the images from pyfits, this implies all the rows for a certain
# range of columns


def flatten(im, start, end):
    ''' flattens the image sections determined by locate_slits'''
    return im[:, start:end].sum(axis=1)


def detect_edges(im, width, sig, thres, full_report):
    '''use the sobel edge detection algorithm to detect the edges of MOS
    spectra
    * take the log of the edges to make it easier to detect the faint edges.
    '''

    # ------TEST-------
    # taking the log10 of the flattened image so that the sigma clipping is
    # improved.
    # ---- Extra debug info

    # make sure there are no 0 values before taking the log
    im[im <= 0.0] = 0.1
    I = numpy.log10(im)

    # check the images before performing the edge detetion
    # smoothedI = smooth_edges(I)
    smoothedI = I

    # perform the edge detection with the sobel image filter
    edges = scipy.ndimage.filters.sobel(smoothedI, mode='wrap')

    # NOT IMPLEMENTED
    # smooth the edge detected array using a median filter with 5 pixel size
    median_filtered_array = scipy.ndimage.filters.median_filter(
        edges, size=(5,))

    # get the median of the edge detected array using a pixel size specified
    med_edges = scipy.ndimage.filters.median_filter(edges, size=(50,))

    # create a Gaussian kernel with a width of 25 pixels and sigma of 2.2
    g = scipy.signal.gaussian(width, sig, sym=1)

    # convolve the image with the Gaussian kernel
    gaussian = scipy.ndimage.filters.convolve(edges, g)
    #  p = numpy.arange(0,len(gaussian),1)

    # clip the Gaussian array to determine where the edges begin
    s_deviation = numpy.std(med_edges)
    clip = (gaussian <= thres * s_deviation) & \
           (gaussian >= - thres * s_deviation)
    SD = numpy.std(gaussian * clip)

    # populate an array with all the left hand edges
    start_edge = (gaussian > 10 * SD) * 1

    # populate an array with all the left hand edges
    end_edge = (gaussian < -10 * SD) * 1

    # determine the coordinates of the left and right edges
    ri = 1
    S_edge = numpy.zeros(len(start_edge))
    for h in range(0, len(start_edge) - 1):
        if (start_edge[h + 1] - start_edge[h]) == 1:
            S_edge[h] = ri
            ri += 1

    li = 1
    E_edge = numpy.zeros(len(end_edge))
    for l in range(0, len(end_edge) - 1):
        if (end_edge[l + 1] - end_edge[l]) == -1:
            E_edge[l] = li
            li += 1

    # S is the array that contains the position of the left edge
    S = []

    # E is the array that contains the position of the right edge
    E = []

    for k in range(1, ri):
        for m in range(0, len(S_edge)):
            if S_edge[m] == k:
                S.append(m)

    for k in range(1, ri):
        for m in range(0, len(E_edge)):
            if E_edge[m] == k:
                E.append(m)

    if full_report:
        return S, E, edges, median_filtered_array, med_edges, S_edge, E_edge, \
            gaussian, start_edge, end_edge, I, smoothedI
    else:
        return S, E


def smooth_edges(I):
    '''
    this function smoothes the edges of the array, which in some cases cause
    problems for edge detection.
    '''

    length = numpy.shape(I)[0]
    start = I[0:40]
    end = I[length - 40:length]

    smoothed_start = scipy.ndimage.filters.median_filter(start, size=(3,))
    smoothed_end = scipy.ndimage.filters.median_filter(end, size=(3,))

    I[0:40] = smoothed_start
    I[length - 40:length] = smoothed_end

    return I


def extract_spectrum(im, left_edge, right_edge):
    ''' extract  a rectangular image section determined by the min and max from
    the spline values of each individual spectrum'''

    return im[left_edge:right_edge, :]


def image_sections(im, debug, inc=1):
    '''
    divide the image into the amount of sections determined from divide_image
    '''

    if debug:
        print 'dividing the image into sections'

    # determine the section size
    inc_size = numpy.shape(im)[1] / inc
    sec = []
    x_coor = []

    if inc == 1:
        flattened_image = flatten_image(im)
        sec.append(flattened_image)
        x_coor.append(1)

    elif inc > 1:

        for i in range(0, inc):
            flattened_image = flatten_image(
                im[:, i * inc_size:(i + 1) * inc_size])
            sec.append(flattened_image)
            x_coor.append(i * inc_size)

    return sec, x_coor


def fit_spline(x, y, new_x, order):
    '''
    this function determines the new y values for the slit positions by fitting
    a spline of order (order)
    '''

    # determine the spline parameters from the slit positions
    tck = si.splrep(x, y, k=order)
    # evaluate the spline using the new x-values
    new_y = numpy.round_(si.splev(new_x, tck, der=0))

    return new_y


def divide_image(im, sections):
    '''define the ranges to be summed when flattening the image sections
    * im - flat image data
    * sections - amount of sections that the image needs to be cut up in. These
    sections are evenly spaced.
    '''
    xw = im.shape[1]
    # add 1 so that it returns an inclusive list
    return range(0, xw + 1, int(xw / sections))


def locate_slits(im, inc, width, sig, thres, full_report=True):
    """
    this function is responsable for running the edge detection algorithm on
    flattened image sections. The slits detected from the flattened image image
    sections are appended to the slits array. All other info from the edge_detection
    function is written in the other_info array. The slits array is returned to the
    main program.

    * im - slit image data array
    * inc - array that contains the x-positions of the image sections to use
    for the flattening of the images.
    """

    slits = []
    other_info = []

    for i in range(0, len(inc) - 1):
        # send the image section to be flattened
        flattened_image = flatten(im, inc[i], inc[i + 1])
        # return all the relevant arrays from edge_detect
        S, E, edges, median_filtered_array, med_edges, S_edge, E_edge, gaussian, start_edge, end_edge, I, smoothedI = detect_edges(
            flattened_image, width, sig, thres, True)

        # append the slit position of the current section to the slits array
        slits.append([S, E])

        # append all the other info of the detected edges to the other_info
        # array
        other_info.append([edges,
                           median_filtered_array,
                           med_edges,
                           S_edge,
                           E_edge,
                           gaussian,
                           flattened_image,
                           I])

    if full_report:
        return slits, other_info
    return slits


def get_slits(slits):
    '''
    this function checks whether the detected slits for the left and right
    edges are the same length and also that the number of slits detected for
    each image section is the same. Additionally it checks whether there is
    an unexpecdedly large deviation in the slit curvature.
    '''

    num_sections = len(slits)

    # set up the calculation for the number of slits calculated
    # for each section
    section = [slits[i] for i in range(0, num_sections)]
    lengths = [[len(section[j][0]), len(section[j][1])]
               for j in range(0, len(section))]

    # first check and see if all the left and right edges are of equal length
    for k in range(0, len(lengths)):
        if lengths[k][0] != lengths[k][1]:
            msg = 'Number of edges for section %i does not match. Try changing your threshold value' % k
            raise SALTSpecError(msg)

    for m in range(0, len(lengths)):
        if lengths[m][0] == 0 or lengths[m][1] == 0:
            msg = 'No edges detected for section %i. Try changing your threshold value' % k
            raise SALTSpecError(msg)

    equal = []
    for i in range(0, len(lengths) - 1):
        if lengths[i][0] != lengths[
                i + 1][0] or lengths[i][1] != lengths[i + 1][1]:
            equal.append(False)
        else:
            equal.append(True)

    if (not all(equal)):
        msg = 'An unequal amount of slits were detected for one of the sections'
        raise SALTSpecError(msg)

    # write the slit positions into the defined slit positions format
    # such that each slit is given by:
    # [slitnumber, [bottom slit positions], [top slit positions]]
    allslits = []
    for j in range(0, len(slits[0][0][:])):
        spline_left = [slits[i][0][j] for i in range(0, len(slits))]
        spline_right = [slits[i][1][j] for i in range(0, len(slits))]
        allslits.append([j, spline_left, spline_right])

    return allslits


def extract_slits(slits, spline_x, im, order=2, padding=0):
    """
    this function extracts the individual spectra from the MOS image using the
    slit positions determined by get_slits. the slit positions are fitted using
    a spline, padded with a zero mask and then extracted.
    """
    y_dim, x_dim = numpy.shape(im)

    # create a list that will contain each of the extracted spectra
    allspec = []
    spline_positions = []

    # create the new x array, sampled at every pixel position
    new_x = numpy.arange(0, x_dim, 1)

    for j in range(0, len(slits)):
        #
        # create a mask for each of the slits
        mask = numpy.ma.make_mask_none((y_dim, x_dim))

        # get the detected left and right edges for each of the sections
        edge_left = slits[j][1]
        edge_right = slits[j][2]

        # fit a spline to the edges
        if len(spline_x) > 1:
            spline_left_edge = fit_spline(
                spline_x,
                edge_left,
                new_x,
                int(order)) - padding
            spline_right_edge = fit_spline(
                spline_x,
                edge_right,
                new_x,
                int(order)) + padding
        else:
            spline_left_edge = new_x * 0.0 + edge_left[0] - padding
            spline_right_edge = new_x * 0.0 + edge_right[0] + padding

        spline_left_edge = spline_left_edge.astype(int)
        spline_right_edge = spline_right_edge.astype(int)

        # create the mask section for the spectra
        for k in range(0, x_dim - 1):
            mask[spline_left_edge[k]:spline_right_edge[k], k] = True

        # mask everything except the region that lie within the spline values
        c = mask * im

        # determine the size of the rectangle to cut out
        miny = numpy.min(spline_left_edge)
        maxy = numpy.max(spline_right_edge)

        # extract the rectangular section from the MOS image
        ext = extract_spectrum(c, miny, maxy)

        # append the extracted data section to the allspec list
        allspec.append(ext)

        # append the spline fitted edges to the slit_position list
        spline_positions.append([spline_left_edge, spline_right_edge])

    return allspec, spline_positions


def slits_HDUtable(slit_pos, order):
    '''
    create the BinaryHDU table for the output image and the slit image
    The fits binary table format for 64-bit floats is K
    '''

    columns = []
    columns.append(Column(name='spline_order', format='K', array=[order]))
    columns.append(Column(name='slitnum', format='K', array=[len(slit_pos)]))
    for i in range(0, len(slit_pos)):
        columns.append(Column(name='slit_%i_left_edge' % i, format='K',
                              array=slit_pos[i][1]))
        columns.append(Column(name='slit_%i_right_edge' % i, format='K',
                              array=slit_pos[i][2]))

    tbhdu = fits.BinTableHDU.from_columns(columns)

    return tbhdu


def read_slits_HDUtable(tableHDU):
    '''
    read the binary fits HDU table and return the extracted slits.
    '''

    d = tableHDU.data
    slitnum = d.field('slitnum')[0]
    order = d.field('spline_order')[0]

    allslits = []
    for i in range(0, slitnum):
        left_edge = [
            d.field(
                'slit_%i_left_edge' %
                i)[j] for j in range(
                0, len(
                    d.field(
                        'slit_%i_left_edge' %
                        i)))]
        right_edge = [
            d.field(
                'slit_%i_right_edge' %
                i)[j] for j in range(
                0, len(
                    d.field(
                        'slit_%i_right_edge' %
                        i)))]
        allslits.append([i, left_edge, right_edge])

    return order, allslits


def write_outputslitfile(slit_pos, outputslitfile, order):
    '''
    write the slits to an ASCII file
    slit_pos: array of splines fitted to the detected edges. The first element
    in the array is the slit number, second element is the left edge array, the
    third element is the right edge.

    outputslitfile: name of the output file
    '''

    # open the output ASCII file for writing
    output_ascii = open(outputslitfile, 'w')

    S = slit_pos

# verify_line = '#$ specslit slit position file\n'
#    output_ascii.write(verify_line)
    # write the spline order to the output file
    order_line = '# Spline fit order : %i\n' % order
    output_ascii.write(order_line)

    # write each line in the slit_pos  array to the output file
    # each element is space seperated and the edges arrays are comma
    # comma speratedS
    for i in range(0, len(S)):
        line = '  '.join([str(S[i][0]),
                          ','.join([str(S[i][1][j])
                                    for j in range(0, len(S[i][1]))]),
                          ','.join([str(S[i][2][j]) for j in range(0, len(S[i][2]))])]) + '\n'

        output_ascii.write(line)

    # close the output file
    output_ascii.close()


def read_slits_from_ascii(simg):
    '''
    this function reads the input from an outputslitfile if it specified as a
    slit image. It will read the slit positions return an array that resembles
    the split_pos array that is passed to the extract_slits function and the
    order of the spline to be fitted.
    '''

    slitfile = open(simg, 'r')
    slits = slitfile.readlines()
    slitfile.close()

    allslits = []
    for line in slits:
        if line[0] == '#':
            order = int(line.split(':')[1])
        else:
            col = line.split()
            slitnum = int(col[0])
            left_edge = [
                float(
                    col[1].split(',')[i]) for i in range(
                    0, len(
                        col[1].split(',')))]
            right_edge = [
                float(
                    col[2].split(',')[i]) for i in range(
                    0, len(
                        col[2].split(',')))]
            allslits.append([slitnum, left_edge, right_edge])

    return order, allslits


def read_slits_from_ds9(simg, order=1):
    '''
    This function reads the input from an outputslitfile if it specified as a
    ds9 region file. It will read the slit positions return an array that resembles
    the split_pos array that is passed to the extract_slits function and the
    order of the spline to be fitted.
    '''

    # read in the slit file
    slitfile = open(simg, 'r')
    slits = slitfile.read()
    slitfile.close()

    # throw an error if the positions aren't in image format
    if slits.count('image') == 0:
        msg = "Please use 'image' for format when saving ds9 region file"
        raise SALTSpecError(msg)
    slits = slits.split('\n')

    # loop through and create all slits
    allslits = []
    alltext = []
    i = 0
    for line in slits:
        if line.startswith('#') or line.startswith('global'):
            pass
        elif line.count('box'):
            line = line.split('#')
            if len(line) == 2:
                rtext = line[1].replace(' text={', '')
                rtext = rtext.replace('}', '')
                alltext.append(rtext.split(','))
            line = line[0].replace('box(', '')
            line = line.replace(')', '')
            col = line.split(',')
            slitnum = i
            left_edge = [int(float(col[1]) - 0.5 * float(col[3]))]
            right_edge = [int(float(col[1]) + 0.5 * float(col[3]))]
            allslits.append([slitnum, left_edge, right_edge])
            i += 1

    return order, allslits, alltext


def which_ext(hdu, ext_name, debug):
    '''
    this function determines which extension in the fits file is the primary, sci
    and binary table. an integer is retunred with the position of the ext in
    the fits file.
    * hdu is an opened fits file. opened with pyfits.open
    * ext_name is the name of the wanted extension
    '''
    if debug:
        print 'determining the fits extension for %s' % ext_name

    if ext_name == 'PRIMARY':
        for i in range(0, len(hdu)):
            if hdu[i].name == 'PRIMARY':
                pos = i
                return pos

    elif ext_name == 'SCI':
        for i in range(0, len(hdu)):
            if hdu[i].name == 'SCI':
                pos = i
                return pos

    elif ext_name == 'BINTABLE':
        for i in range(0, len(hdu)):
            if hdu[i].name == 'BINTABLE':
                pos = i
                return pos

    else:
        msg = 'the hdu was not found in the fits header'
        raise SALTSpecError(msg)


def show_image_with_edges(im, splines, img, debug):

    if debug:
        print 'plotting the extracted spectra on the slit image'

    y_dim, x_dim = numpy.shape(im)
    x = numpy.arange(0, x_dim, 1)
    # collopse the image on the y-axis

    # start with a rectangular Figure
    pl.figure(num=None, figsize=(14, 10))

    z1, z2 = saltimagetools.zscale(im, 0.25)
    colmap = pl.cm.Greys

    pl.imshow(im, cmap=colmap, vmin=z1, vmax=z2, origin='lower')
    for i in range(0, len(splines)):
        left_edge = splines[i][0]
        right_edge = splines[i][1]
        pl.plot(x, left_edge, 'b-')
        pl.plot(x, right_edge, 'r-')

    pl.xlim((-10, x_dim + 10))
    pl.ylim((-10, y_dim + 10))
    pl.title('Slits detected on %s' % img)
    pl.xlabel('x - pixels')
    pl.ylabel('y - pixels')

    pl.show()


def show_image(im, fignum=None, figsize=(19, 10)):
    z1, z2 = zscale(im)
    pl.figure(num=fignum, figsize=figsize)
    colmap = pl.cm.Greys
    xmax, ymax = numpy.shape(im)
    ext = [0, xmax, 0, ymax]
    plot = pl.imshow(
        im,
        cmap=colmap,
        vmin=z1,
        vmax=z2,
        origin='lower',
        aspect='auto')

    return plot


def convert_fromsky(ra, dec, cra, cdec, position_angle=0, equinox=2000, xpixscale=0.125,
                    ypixscale=0.125, ccd_cx=3121, ccd_cy=2051):
    """Convert from ra,dec coordindates to x,y coordindate on CCD.  Currently hardwired
       but should use the CCD information or an RSS model for more accurate results
    """

    wcs = pywcs.WCS(naxis=2)
    wcs.wcs.crpix = [ccd_cx, ccd_cy]
    wcs.wcs.cdelt = numpy.array(
        [-xpixscale, ypixscale]) / 3600.  # set in degrees
    wcs.wcs.crval = [cra, cdec]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crota = [position_angle, position_angle]  # rotate SKY this amount?
    wcs.wcs.equinox = equinox

    # convert onsky coords to pix
    xp, yp = wcs.wcs_sky2pix(ra, dec, 1)

    return xp, yp


def distortion_correct(
        x, y, x_d=3235.8, y_d=2096.1, A=0.0143, B=0.0078, r_d=1888.8):
    """Distortion correct pixel position according to Ken's model
       where
       x'=x_d+(x-x_d)*(1+A*(r/r_d)**2+B*(r/r_d)**4)
       y'=y_d+(y-y_d)*(1+A*(r/r_d)**2+B*(r/r_d)**4)
    """
    # calculate the radius of the source
    r = ((x - x_d) ** 2 + (y - y_d) ** 2) ** 0.5

    # calucluate the x and y corrected positions
    xi = x_d + (x - x_d) * (1 + A * (r / r_d) ** 2 + B * (r / r_d) ** 4)
    yi = y_d + (y - y_d) * (1 + A * (r / r_d) ** 2 + B * (r / r_d) ** 4)

    return xi, yi


def convert_slits_from_mask(
        slitmask, order=1, xbin=2, ybin=2, pix_scale=0.1267, cx=0, cy=0):
    """Given a slitmask, return the x and y position
       for each of the slits on the CCD.

       The format for slit positions is [id,
    """
    slit_positions = []

    # determine the slit center position
    cra = slitmask.center_ra
    cdec = slitmask.center_dec

    # set up the x- and y-pixscale
    xpixscale = pix_scale * xbin
    ypixscale = pix_scale * ybin

    # convert the slitlets in the slit mask to slits for extraction
    for i in range(slitmask.slitlets.nobjects):
        sid = slitmask.slitlets.data[i]['name']
        sra = slitmask.slitlets.data[i]['targ_ra']
        sdec = slitmask.slitlets.data[i]['targ_dec']
        slen1 = slitmask.slitlets.data[i]['len1']
        slen2 = slitmask.slitlets.data[i]['len2']
        spriority = slitmask.slitlets.data[i]['priority']
        if spriority == -1:
            slen1 = 2.5
            slen2 = 2.5

        sx, sy = convert_fromsky(sra, sdec, cra, cdec, -slitmask.position_angle,
                                 xpixscale=xpixscale, ypixscale=ypixscale, ccd_cx=cx, ccd_cy=cy)
        sx, sy = distortion_correct(sx, sy, x_d=3235.8 /
                                    xbin, y_d=2096.1 /
                                    ybin, A=0.0143, B=0.0078, r_d=1888.8 /
                                    xbin)
        sy1 = sy + slen1 / ypixscale
        sy2 = sy - slen2 / ypixscale
        slit_positions.append([sid, sy2, sy1])

    # return the values
    return order, slit_positions


def read_slitmask_from_xml(sxml):
    """Read the slit information in from an xml file
       that was used to create the slitmask
    """

    # read in the xml
    dom = minidom.parse(sxml)

    # create the slitmask
    slitmask = SlitMask()

    # read it in
    slitmask.readmaskxml(dom)

    # return
    return slitmask


def parsexml(xmlfile):
    dom = minidom.parse(xmlfile)
    slits = dom.getElementsByTagName('slit')
    refstars = dom.getElementsByTagName('refstar')
    s1 = slits[0]  # need to loop over slits
    s1.getAttribute('xce')
    float(s1.getAttribute('yce'))

    parameters = dom.getElementsByTagName('parameter')
    p1 = parameters[2]  # need to loop over the parameters
    p1.getAttribute('name')

    # to convert the unicode to strings just do:
    # str(p1.getAttribute('name'))


def extract_fit_flux(data, mask=None, err=None, padding=1, bounds=True):
    """For each column, fit the flux in that column to a model and return the flux.  This fits a 
    flat baseline plus a guassian to the summed line profile and then fixes the mean and stddev
    of the Guassian profile.  The amplitude and baseline are then fit to each wavelength in the 
    data.

    Prior to fitting, the function examines the slit to only include good data.
    
    Parameters
    ----------
    data: np.ndarray
        2-D array of fluxes

    mask: np.ndarray
        Array of masked pixels (optional)
 
    err: np.ndarray
        Array of error values (optional)
 
    padding: int
        Padding around identify edges

    bounds: boolean
        Determine the lower and upper bounds.  If False, uses the entire array.
        
    
    Returns 
    -------
    flux
    """
    
    #determine the upper and lower bounds of where the data are useful
    if bounds:
       s = data.sum(axis=1)
       t = numpy.gradient(s)/numpy.median(s)
       v = numpy.gradient(numpy.gradient(s)/numpy.median(s))
       i_min = t.argmin()
       i_max = t.argmax()
       i_max = i_max + numpy.where(t[i_max:]<0.05)[0][0] + padding
       i_min = numpy.where(t[:i_min]>-0.05)[0][-1] - padding
       c = data[i_max:i_min, :]
       if err is not None:
           e = err[i_max:i_min, :]
       if mask is not None:
           k = err[i_max:i_min, :]
    else:
       c = 1.0 * data
       if err is None: 
          e = 0.0 * c + 1.0 
       else: 
          e = 1.0 * err
       if mask is None: 
          k = 0.0 * c
       else:
          k = 1.0 * mask
    
    #build the best fit model 
    sarr = c.sum(axis=1)
    xarr = numpy.arange(len(sarr))
    fitter=mod.fitting.LevMarLSQFitter()
    m_init=mod.models.Gaussian1D(amplitude=sarr.max(), mean=xarr[sarr.argmax()], stddev=5) + mod.models.Polynomial1D(0)
    m = fitter(m_init, xarr, sarr)
    m.stddev_0.fixed=True
    m.mean_0.fixed=True
    
    
    # fit each row to the model 
    ys, xs = c.shape
    flux = numpy.zeros(xs)
    sky = numpy.zeros(xs)
    for i in range(xs):
        if mask is not None:
           l = (k[:,i]==0)
           if l.sum() > 0.5*len(l): # at least half the pixels must
              n = fitter(m, xarr[l], c[:,i][l], weights=1.0/(e[:,i][l]))
              flux[i] = n.amplitude_0.value
              sky[i] = n.c0_1.value
           else:
              flux[i] = numpy.nan
              sky[i] = 0.0
        else:
           n = fitter(m, xarr, c[:,i])
           flux[i] = n.amplitude_0.value
           sky[i] = n.c0_1.value
        
    return flux, sky
