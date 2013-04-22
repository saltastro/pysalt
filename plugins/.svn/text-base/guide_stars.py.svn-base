#!/usr/bin/env python

# Author                   Version      Date
# -----------------------------------------------
# Keith Smith (Nottingham)   0.8        23 Nov 2007

# Module for finding SALT guide star
# - full changelog is in /doc/changelog.txt
# - manual is (will be?) in /doc/manual.txt


__version__ = "0.8"
__author__ = "Keith Smith"

__doc__="\nSALT guide star finder, version "+__version__ +"""

Finds guide stars suitable for use with the Southern African Large Telescope
Uses the HST Guide Star Catalogue 2.3.2 and the VizieR catalogue service

Usage: python guide_stars.py [OPTIONS] [TARGET]

TARGET should be the coordinates of the target or a SIMBAD-resolvable name
Acceptable formats are:
    colon-seperated sexagesimal eg. 15:43:17.2 -18:56:12.9
    space-seperated sexagesimal eg. '15 43 17.2' '-18 56 12.9'
    decimal eg. 15.7837 -18.9543
    target name eg. m31 or 'RV Cen'
All coordinates should be J2000 right ascension and declination

OPTIONS are as follows, arguments are compulsory for both long and short forms:
        --help              Prints this help
    -v    --verbose           Prints useful information during execution
    -d    --debug             Prints debugging information, implies -v
    -f    --filter=FILTER     Uses the FILTER filter,
                                defaults to V
    -i    --instrument=INS    Uses the INS instrument,
                                defaults to RSS
    -r    --radius=RADIUS     Excludes RADIUS arcseconds around target,
                                defaults to 2 arcseconds"""


# import required modules

import urllib2              # reading URLs
from xml.dom import minidom # XML parsing
import csv, StringIO        # CSV parsing
from numpy import *         # array
import sys, getopt                  # command line switches
import warnings

from salttime import dec2sex

# define exceptions

class GuideStarError(Exception): pass    # misc errors in this module
class VizError(Exception): pass # errors thrown by VizieR
#class URLError(Exception): pass # errors accessing the URL
class BadInput(Exception): pass # incorrect input


# define global variables, later read from the command line



# define helper functions

def usage():
    print __doc__
    raise SystemExit(2)        # 2 is the UNIX code for bad command line input apparently

def isDecimal(string):
    "returns true if the input string converts to a valid decimal"
    try:
        float(string)
        return True
    except ValueError:
        return False




def checkdms(dms):
    """Verify a sexagesimal string; returns True if valid, False if not
    """
    assert isinstance(dms, str)

    try:
        d=dms.split(':')
        for i in range(3): float(d[i])
        return True
    except:
        return False


def checkInput(targetRA, targetDec, imfilter, instrument, targetRadius, maxRadius=300):
    "validates input to the main function"

    if isDecimal(targetDec):
        # dec is decimal
        if isDecimal(targetRA)==False:  # check RA in same format
            raise BadInput, 'RA and dec appear to be in different formats'
        if targetDec>0: 
           targetDec='+'+dec2sex(float(targetDec))
        else:
           targetDec=dec2sex(float(targetDec))
        targetRA=dec2sex(float(targetRA)/15.0)
    elif isDecimal(targetRA):
        # RA is in decimal, but dec wasn't
        # pretty sure there are no valid target names which are decimals...
        # need to check this as decimals pass as valid sexa!
        raise BadInput, 'RA and dec appear to be in different formats'
    elif checkdms(targetDec) == True:
        # target already colon-seperated value
        # couldn't check this first as this function returns True
        # for decimals (for some reason), but False if space seperated sexa
        if checkdms(targetRA) == False:   # check RA in same format
            raise BadInput, 'RA and dec appear to be in different formats'
        pass    # both already in correct format
    elif checkdms(targetDec.replace(' ',':')) == True:
        # is a valid space seperated sexa, convert to colon seperated
        if (checkdmsStr(targetRA.replace(' ',':')) == False or
            (targetRA.replace(' ',':')==targetRA)):
            # check RA works colon seperated and wasn't to start with
            raise BadInput, 'RA and dec appear to be in different formats'
        
        targetDec=targetDec.replace(' ',':')
        targetRA=targetRA.replace(' ',':')
    else:
        # dec isn't decimal, or space/colon seperated sexa, or blank
        raise BadInput, 'Format of declination not recognised, was given: '+str(targetDec)

    # convert filters to those in the GSC
    if imfilter == "":
        warnings.warn('No filter specified, defaulting to V')
        imfilter = 'V'
    elif imfilter in ('j', 'F', 'V', 'N'):
        pass    # these are the native filters in the GSC
    elif imfilter in ('U', "u'", 'u'):
        warnings.warn('Specified '+ imfilter + ' filter, but GSC has very little data this far blue. Falling back to photographic Bj')
        imfilter = 'j'    # there IS a U field in GSC, but hardly any entries
    elif imfilter in ('B', "g'", 'b', 'v'):
        warnings.warn('Specified ' + imfilter + ' filter, closest GSC band is photographic Bj')
        imfilter = 'j'   # there is also a Bmag field in GSC, but less entries
    elif imfilter in 'y':
        warnings.warn('Specified ' + imfilter + ' filter, closest GSC band is photographic V')
        imfilter = 'V'
    elif imfilter in ('R', "r'"):
        warnings.warn('Specified ' + imfilter + ' filter, closest GSC band is photographic F')
        imfilter = 'F'
    elif imfilter in ('I', "i'", "z'"):
        warnings.warn('Specified ' + imfilter + ' filter, closest GSC band is photographic N')
        imfilter = 'N'
    else:
        raise BadInput, "Filter '%s' not recognised" % str(imfilter)

    # check instrument input
    instrument = instrument.lower()
    if instrument == "":
        warnings.warn('No instrument specified, defaulting to RSS')
        instrument = 'rss'
    elif instrument in ['rss','pfis']:
        instrument = 'rss'
    elif instrument in ["salticam","scam"]:
        warnings.warn('Selected SALTICAM; SALTICAM values not yet validated')
        instrument = 'scam'
    elif instrument == "hrs":
        warnings.warn('Selected HRS; HRS values not yet validated')
        instrument = 'hrs'
    else:
        raise BadInput, 'Instrument "' + str(instrument) + '" not recognised'

    # check radius
    if targetRadius == '':
        warnings.warn('No target radius specified, defaulting to 2 arcsec')
        targetRadius = 2.
    elif 0 < targetRadius < maxRadius:
        pass
    elif targetRadius > maxRadius:
        raise BadInput, 'Target radius '+str(targetRadius)+' arcsec is larger than the science FoV'
    else:
        raise BadInput, 'Target radius of ' + str(targetRadius) + 'arcsec is invalid'

    return (targetRA, targetDec, imfilter, instrument, targetRadius)



def queryUrl(url):
    "accesses the input url and returns the response"

    request = urllib2.Request(url)
    opener = urllib2.build_opener()
    request.add_header('User-Agent', 'SALT guide star finder/'+ __version__ +' www.salt.ac.za')
    try:
        data = opener.open(request).read()
    except:
        raise GuideStarError, "Could not connect to VizieR. Please check your internet connection"
    return data

def constructVizUrl(targetRA, targetDec, min_r, max_r, imfilter, min_mag, max_mag):
    "constructs the url for the VizieR guide star query"

    url = "http://vizier.u-strasbg.fr/cgi-bin/asu-xml?" # base URL for VizieR XML queries
    url += "-source=I/305/out"                  # use HST GSC 2.3.2
    url += "&-c=" + targetRA
    url += targetDec + "&-c.eq=J2000"   # target section
    url += "&-c.rm=%s,%s" % (min_r, max_r)      # anulus range (arcmin)
    url += "&-out=%smag&%smag=%s..%s" % (imfilter, imfilter, min_mag, max_mag)
    url += "&-out.max=1000"                        # max entries
    url += "&-out=_r,_RA*-c.eq,_DE*-c.eq"       # calculate 2000 RA, dec and distance
    url += ",Class&Class=0"                     # only objects flagged as stars
    url += "&-oc.form=sexa"                     # output coords in sexagesimal
    url += "&-sort=-_r"                         # sort by decreasing r
    url += "&-mime=CSV"                         # data in CSV format

    return url



def parseVizResponse(response):
    "parses the response from VizieR into into a list of data for each star"

# extract CSV table from XML

##    if "****" in response:
##        # actually puts this in
##        #<INFO ID="Errors" value="(**** indicates an Error,++++ indicates a Warning)">
##        #blah blah blah
##        #</INFO>
##        # could handle this better with a little parsing
##        errordata=""
##        for line in response:
##            if "****" in response:
##                errordata += line
##        raise VizError, "Vizier generated an error. The data returned was:\n" + errordata


    xml = minidom.parseString(response)

    table = xml.getElementsByTagName('CSV')
    if len(table)==0:
        info = xml.getElementsByTagName('INFO')
        for element in info:
            id = element.getAttribute('ID')
            if id=='Errors':
                #found an error
                raise VizError, 'Vizier generated an error. The error was:' + element.firstChild.data
        # no CSV table was found ie. no data
        # needs some error checking - might be a bad url or target
        return [],0,[]
    colsep = table[0].getAttribute('colsep')        # CSV columns seperator
    headlines = table[0].getAttribute('headlines')  # CSV header lines
    table = table[0].firstChild.data                # the CSV table itself

    # convert to from unicode
    colsep = colsep.encode('ASCII')
    headlines = int(headlines)


# extract data from CSV table

    csvreader=csv.reader(StringIO.StringIO(table), delimiter=colsep)


# this next bit is very kludgy, weird csvreader object and array manipulations
    for row in csvreader:   # can't index csvreader!
        if len(row)>0:  # blank rows in the reader!
            if 'csvtable' in locals():      # checks to see if csvtable exists
                # python has no exist() function!
                csvtable=vstack((csvtable,array([row])))
                # fragile, dimension mismatch possible
            else:
                csvtable=array([row])
                # can't append to blank arrays!

# old way of using dictionaries:
        #if len(row)>=4:
            #rowdict = {0:row[0],1:row[1],2:row[2],3:row[3],4:row[4]}
            #csvtable.append(rowdict)   # last column Class = 0 for all

    #csvtable = csvtable[int(headlines):] # strip headers


#    header = csvtable[:headlines]
#    data = csvtable[headlines:]

    n_stars=len(csvtable)-headlines

# split off the Class column, as an added bonus won't break if there isn't one
    index = (csvtable[0,:]!='Class')  # boolean index of first row != Class
    csvtable = csvtable[:,index]        # just retain columns in the index

    return csvtable, n_stars, headlines



def sortResults(table, n_stars, headlines, min_r, max_r, min_mag, max_mag):
    "Sorts the results to select the best guide stars"

    # split the header off the table
    header=table[:headlines]
    data=table[headlines:]

    top=header[0,:]

    r_index = top=='_r'
    if not True in r_index:
        raise GuideStarError, ('Could not find a radius column in search results\n' +
              'Returned columns are: ' + str(top))
    r=data[:,r_index]

    mag_index = (top=='Vmag') + (top=='Nmag') + (top=='Bmag') +\
                (top=='Fmag') + (top=='Umag') + (top=='jmag')
    # hard coded brute force
    if not True in mag_index:
        raise GuideStarError, ('Could not find a magnitude column in search results\n' +
              'Returned columns are: ' + str(top))
    mag=data[:,mag_index]

    # convert to floats
    r=r.astype(float64)
    mag=mag.astype(float64)

    r_pref=select([r<min_r, (min_r<=r) & (r<=max_r), r>max_r],
                  [min_r-r, 0., r-max_r])

    mag_pref=select([mag<min_mag, (min_mag<=mag) & (mag<=max_mag), mag>max_r],
                  [min_mag-mag, 0., mag-max_mag])

    pref = 2*r_pref + mag_pref    # weighted such that 1 mag = 0.5 arcmin

    sort_index = lexsort(keys=(mag.T,pref.T),axis=1)    # decide the order from the pref list
                                    # in order, so LOWEST first
    # note the index is crazy - lists the indexes in order, not the order of indices!
    sorted = data[sort_index]
    pref_sort = pref[sort_index]

    n_stars = min(n_stars,6)    # return at most 6 stars


    # replace some of the column headings
    top[top.tolist().index('_r')]='Offset'
    top[top.tolist().index('_RAJ2000')]='RA (J2000)'
    top[top.tolist().index('_DEJ2000')]='Dec (J2000)'

    header[0,:]=top

    # put the header back on
    data = vstack((header,sorted[:,:n_stars][0]))

    return data, n_stars

def QueryCatalog(targetRA, targetDec, abs_min_r=0, abs_max_r=10, imfilter='V', 
    abs_min_mag=0, abs_max_mag=30, header=True):
    '''Finds guide stars based upon the HST GSC 2.3.2.
    Pass target RA and Dec (J2000), and optionally
    filter, instrument, and target radius in arcsec.
 
    Parameters
    ----------
    targetRA: string 
       string specifying the target RA.  The format should either 
       be the decimal position or a colon or space separated values

    targetDec: string 
       string specifying the target RA.  The format should either 
       be the decimal position or a colon or space separated values
   
    abs_min_r: float
       Minimum radius in arcminutes for search annulus

    abs_max_r: float
       Maximum radius in arcminutes for search annulus

    imfilter: string
       Input filter for the observations

    abs_min_mag: float
       Minimum magnitude for stars (bright limit)

    abs_max_r: float
       Maximum mangitude for stars (faint limit)

    header: boolean
       Include header in array

    Returns
    -------

    data: array
       Array containing the sorted results from the query

    n_stars: int
       Number of stars returned

    headlines: int [optional]
       Number of lines in the header.  Only returned if header=True

    '''
 
    # set up query
    url = constructVizUrl(targetRA, targetDec, abs_min_r, abs_max_r, imfilter, abs_min_mag, abs_max_mag)

    # retrieve data
    response = queryUrl(url)

    #sort the data
    data, n_stars, headlines = parseVizResponse(response)

    if not header:
       return data[headlines:], n_stars
    return data, n_stars, headlines

# main function

def findGuideStars(targetRA, targetDec="", imfilter='', instrument='', targetRadius=2.):
    '''Finds SALT-suitable guide stars based upon the HST GSC 2.3.2.
    Pass target RA and Dec (J2000), and optionally
    filter, instrument, and target radius in arcsec.
 
    Parameters
    ----------
    targetRA: string 
       string specifying the target RA.  The format should either 
       be the decimal position or a colon or space separated values

    targetDec: string 
       string specifying the target RA.  The format should either 
       be the decimal position or a colon or space separated values

    imfilter: string
       Input filter for the observations

    instrument: string
       Instrument to be used for guiding: either rss or scam
 
    targetRadius:  float
       Radius in arcsec to search for targets.   Maximum value should be 300"
    
    '''

    (targetRA, targetDec, imfilter, instrument, targetRadius)= \
        checkInput(targetRA, targetDec, imfilter, instrument, targetRadius)

    # set up instrument-specific settings
    # this should go into the check input function, but too much passing
    # for the moment - will change later
    # alternatively use classes?
    instrument = instrument.lower()
    if instrument == 'rss':
        pref_max_r = 5.     # prefered outer edge of anulus in arcmin
        pref_min_r = 4.     # prefered inner edge of anulus in arcmin
        abs_max_r = 5.      # absolute limit on outer edge
        abs_min_r = 1.      # absolute limit on inner edge
        pref_max_mag = 15.  # prefered magnitude limits
        pref_min_mag = 12.
        abs_max_mag = 19.   # absolute magnitude limits
        abs_min_mag = 10
    elif instrument == 'scam':
        # salticam guider has not yet been installed - use values from rss
        # if you fill this in, please change warning message in input check above
        pref_max_r = 5.
        pref_min_r = 4.
        abs_max_r = 5.
        abs_min_r = 1.
        pref_max_mag = 15.
        pref_min_mag = 12.
        abs_max_mag = 19.
        abs_min_mag = 10.
    elif instrument == "hrs":
        # hrs has not yet been built - use values from rss
        # if you fill this in, please change warning message in input check above
        pref_max_r = 5.
        pref_min_r = 4.
        abs_max_r = 5.
        abs_min_r = 1.
        pref_max_mag = 15.
        pref_min_mag = 12.
        abs_max_mag = 19.
        abs_min_mag = 10.
    else:
        raise GuideStarError, "Instrument settings not found, but passed input checking. This shouldn't happen"


    # modify radii if the source is very large
    targetRadius = (targetRadius / 60.) # convert from arcsec to arcmin
    if targetRadius > pref_max_r:
        if targetRadius > abs_max_r:
            return "Target is larger than the guide star field of view"
        else:
            pref_max_r = abs_max_r
    if targetRadius > abs_min_r:
        abs_min_r = targetRadius    # + a bit to avoid vignetting?
        if abs_min_r > pref_min_r:
            pref_min_r = abs_min_r

    gsdata = []

    #query the catalog and return the results as an array
    data,n_stars, headlines=QueryCatalog(targetRA, targetDec, abs_min_r, abs_max_r, imfilter, abs_min_mag, abs_max_mag)

    #sort the stars
    if n_stars>0:
        data, n_stars = sortResults(data, n_stars, headlines, \
                                pref_min_r, pref_max_r, pref_min_mag, pref_max_mag)

    if n_stars>1:
        print('Finished, selected ' + str(n_stars) + ' guide stars')
        status = 'Found ' + str(n_stars) + ' stars'
    elif n_stars==1:
        print('Finished, found one guide star')
        status = 'Found 1 star'
    else:
        print('Failed, found no guide stars')
        if imfilter != 'V':
            print('Falling back to V filter...')
            status, data = findGuideStars(targetRA, targetDec, 'V', instrument, targetRadius*60.)
            status = status + ' after falling back to V'
        else:
            print('Nothing else to fall back to. Try modifying your query')
            status = 'Failed, no suitable stars found'


    return data


if __name__ == "__main__":
   # executes if module is run from the command line

   # read command line options
   try:
       opts,args = getopt.getopt(sys.argv[1:],"vdf:i:r:",
           ["verbose","debug","filter=","instrument=","radius=","help"])
   except getopt.GetoptError, inst:
       print inst
       print 'Use --help to get a list of options'
       sys.exit(2)

   ra, dec, imfilter, ins, radius = "","","","",""

   # parse them to the relevant variables
   for opt, arg in opts:
       if opt in ('--help'):
            usage()
       elif opt in ('-v','--verbose'):
            verbose=True
       elif opt in ('-d','--debug'):
            verbose=True    # implied
            debug=True
       elif opt in ('-f','--filter'):
            imfilter = arg
       elif opt in ('-i','--instrument'):
            ins = arg
       elif opt in ('-r','--radius'):
            radius = float(arg)
       else:
            print 'Unknown option: ' + opt
            usage()

   for argument in args:
        if ra=="":
            ra = argument
        elif dec=="":
            dec = argument
        else:
            #too many arguments
            raise BadInput, 'Too many arguments, takes one or two input arguments'

   if ra=="":    # no target was specified
        raise BadInput, 'No target specified'

   n_stars,data = findGuideStars(ra,dec,imfilter,ins,radius)
   print n_stars, data
   sys.exit(0)



### Testing stuff
#    print "Debugging mode: using default values"
##    verbose=True
##    debug=True
##    import simbad
##    name="hudf"
###    (ra, dec) = simbad.simbad(name)
##    (ra, dec) = (name, "")
###    (ra, dec) = ('80.9','80.9')
##    radius=2.
##    filter="B"
##    ins = 'rss'
##    status, data = findGuideStars(ra, dec, filter=filter, instrument=ins, targetRadius=radius)
##    print data
###    print status
###    for line in data:
#        print line
