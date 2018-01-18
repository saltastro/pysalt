#!/usr/bin/env python
import sys
import argparse

from astropy.io import fits


header_dict={'PROPID':'50A', 'PROPOSER':'20A', 'OBJECT':'100A', 'RA':'12A', 'DEC':'12A', 'EPOCH':'E', 'EQUINOX':'E', 'DATE-OBS':'10A', 'UTC-OBS':'12A', 'TIME-OBS':'12A', 'EXPTIME':'D', 'OBSMODE':'20A', 'DETMODE':'20A', 'CCDTYPE':'8A', 'NCCDS':'I', 'CCDSUM':'5A', 'GAINSET':'6A', 'ROSPEED':'4A', 'INSTRUME':'8A', 'TELHA':'11A', 'TELRA':'11A', 'TELDEC':'12A', 'TELPA':'E', 'TELAZ':'E', 'TELALT':'E', 'TRKX':'E', 'TRKY':'E', 'TRKZ':'E', 'TRKPHI':'E', 'TRKTHETA':'E', 'TRKRHO':'E', 'COLPHI':'E', 'COLTHETA':'E', 'TELTEM':'E', 'PAYLTEM':'E', 'AMPTEM':'E', 'DETSWV':'16A', 'BLOCKID':'E', 'BVISITID':'E'}

rss_dict = {'FILTER':'8A', 'LAMPID':'8A', 'CALFILT':'8A', 'CALND':'E', 'PELLICLE':'8A', 'INSTPORT':'8A', 'CF-STATE':'20A', 'SM-STATE':'20A', 'SM-STA':'8A', 'SM-STEPS':'J', 'SM-VOLTS':'E', 'SM-STA-S':'E', 'SM-STA-V':'E', 'MASKID':'16A', 'MASKTYP':'16A', 'WP-STATE':'20A', 'HWP-CMD':'16A', 'HW-STEPS':'J', 'HWP-STA':'E', 'QWP-CMD':'16A', 'QW-STEPS':'J', 'QWP-STA':'E', 'QWP-ANG':'E', 'HWP-ANG':'E', 'SH-STATE':'20A', 'FO-STATE':'20A', 'FO-POS':'E', 'FO-VOLTS':'E', 'FO-POS-S':'E', 'FO-POS-V':'E', 'GR-STATE':'20A', 'GR-STA':'10A', 'GR-ANGLE':'E', 'GM-STEPS':'J', 'GM-VOLTS':'E', 'GR-STA-S':'E', 'GR-STA-V':'E', 'GR-STEPS':'J', 'GRATING':'8A', 'GRTILT':'E', 'BS-STATE':'24A', 'FI-STATE':'20A', 'FI-STA':'7A', 'FM-STEPS':'J', 'FM-VOLTS':'E', 'FM-STA-S':'E', 'FM-STA-V':'E', 'AR-STATE':'24A', 'AR-STA':'16A', 'CAMANG':'E', 'AR-STA-S':'E', 'AR-ANGLE':'E', 'PROC':'20A', 'PCS-VER':'4A', 'WPPATERN':'20A','CCDTEM':'D', 'DEWTEM':'E', 'CENTEM':'E'} 

hrs_dict = {'DETSIZE':'17A', 'DETNAM':'10A', 'DETSER':'10A', 'I2STAGE':'10A', 'EXP-TOT':'E', 'EXP-MEAN':'E', 'EXP-MID':'E', 'NODSHUFF':'E', 'NODPER':'E' , 'NODCOUNT':'E', 'PRE-DEW':'E' , 'PRE-VAC':'E' , 'FOC-BMIR':'E', 'FOC-RMIR':'E', 'TEM-AIR':'E' , 'TEM-VAC':'E' , 'TEM-RMIR':'E', 'TEM-COLL':'E', 'TEM-RCAM':'E', 'TEM-BCAM':'E', 'TEM-ECH':'E' , 'CCDTEMP':'E', 'TEM-OB':'E' , 'TEM-IOD':'E'}

scam_dict = {'FILPOS':'I', 'FILTER':'8A', 'CCDTEM':'D', 'DEWTEM':'E', 'CENTEM':'E'}

def create_header_dict_from_list(instrument):
    """From a list of header keywords and formats, create the dictionary of keywords

    Parameters
    ----------
    instrument: str
        Name of the instrument
 
    Returns
    -------
    fits_header_dict: dict
        Dictionary of keywords
    
    """ 
    fits_header_dict=header_dict
    if instrument == 'RSS':
       fits_header_dict.update(rss_dict)
    if instrument == 'HRS':
       fits_header_dict.update(hrs_dict)
    if instrument == 'SCAM':
       fits_header_dict.update(scam_dict)
    return fits_header_dict
   

   

def create_header_dict_from_sdb(instr, sdb):
    """
    """


def fits_header_check(image, fits_header_dict=None, missing=False):
    """Check the header values in the image

    This task will check the fits header values in the image
    and confirm that all header entries are present and that 
    they have an appropriate value

    Parameters
    ----------
    image: str
       Name of an input image

    fits_header_dict: None, dict,  or sdb_mysql
       If None, fits_header_dict will be created from the header list. If an sdb_mysql instance,
       it will be created from the sdb for the instrument.  If a dictionary, it will use that 
       dictonary

    missing: boolean
       If True, it will report on keywords in the FITS file but not in the list.

    Returns
    -------
    fits_header_check: dict
       A dictionary of all FITS keywords that are missing or incorrect.
       It will return an empty dictionary if everything is correct
    """
  
    # open the file 
    hdu = fits.open(image) 

    # determine the instrument
    instrument = hdu[0].header['INSTRUME']

    if instrument not in ['RSS', 'SALTICAM', 'HRS']:
       raise TypeError('{} is not for a SALT instrument or does not have an appropriate instrument keyword'.format(instrument))
 
    # create the fits_header_dict to compare with
    if fits_header_dict is None:
       fits_header_dict = create_header_dict_from_list(instrument)
    #elif isinstance(fits_header_dict, dict):
       #pass
    #elif isinstance(fits_header_dict, object):
    #   pass
    else:
       raise TypeError('{} is not None, dict, or sdb instance'.format(fits_header_dict))


    # check for any missing or incomplete headers
    missing_list=[]
    empty_list=[]
    for key in fits_header_dict:
        if key in hdu[0].header:
           value = hdu[0].header[key]
           if not value:
              if fits_header_dict[key] not in ["J", "E", "D", "I"]: empty_list.append(key)
        else:
           missing_list.append(key)
           
    # check instrument specific keywords that have to be specified
    wrong_list=[]
    #if hdu[0].header['CCDTYPE'] Idu[0].header['RA']

    if instrument == 'RSS':
       # check LAMPID
       if hdu[0].header['CCDTYPE'] == 'ARC' and hdu[0].header['LAMPID'].strip()=='NONE':
          wrong_list.append('LAMPID')
       if hdu[0].header['CCDTYPE'] == 'FLAT' and hdu[0].header['LAMPID'].strip()=='NONE':
          wrong_list.append('LAMPID')

    if instrument == 'HRS':
       if hdu[0].header['OBJECT'] == 'Bias' and hdu[0].header['OBSTYPE']!='Bias':
          wrong_list.append('OBSTYPE')

    # check for extension keywords

    # if missing, check for any headers not in the list
    absent_list=[]
    if missing:
       for key in hdu[0].header:
           if key not in fits_header_dict:
               absent_list.append(key)

    hdu.close()

    if missing_list or empty_list or absent_list or wrong_list:
       return missing_list, empty_list, wrong_list, absent_list

    return 

if __name__=='__main__':
   parser = argparse.ArgumentParser(description='Check SALT FITS Header')
   parser.add_argument('infile', help='SALT HRS image')
   parser.add_argument('-m', dest='missing', default=False, action='store_true', help='Warn about missing keywords')
   parser.add_argument('-d', dest='database', default=False, action='store_true', help='Use database keywords')
   args = parser.parse_args()

   infile = args.infile

   fits_header_dict = None
   if args.database:
      fits_header_dict='sdb'

   results = fits_header_check(infile, fits_header_dict=fits_header_dict, missing=args.missing)
   if results==None:
      exit()

   # print out results
   missing, empty, wrong, absent = results

   hdu = fits.open(infile)
   print('{} {}'.format(infile, hdu[0].header['OBJECT']))
   if missing: print("Keywords that are missing: {}\n".format(missing))
   if empty: print("Keywords that are empty: {}\n".format(empty))
   if wrong: print("Keywords that are wrong: {}\n".format(wrong))
   hdu.close()
