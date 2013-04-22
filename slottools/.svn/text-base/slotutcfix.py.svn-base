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

"""SLOTUTC is a tool to fix the UTC time in SALTICAM slot mode imaging.  Because 
of a timing errror in OS being used, the clock loses time between seconds and
the UTC time is inaccurate between updates of the system clock (which is 
currently updated at the beginning of each second.

As such this program reads in a set of slotmode data and calculates what the 
true exposure plus deadtime is.  Once calculated, it will go through the images 
and correct the UTC times to the correct value.  This process assumes that the 
UTC time for an image right after the second turns is a fiducial time.  In 
addition, the program will fix the timing for a series of frames where the 
timing of the frames is based on the reading out the wrong number of frames.

20080119
    * Fixed bug that wrote negative times
20080217
    * Added the addition fixed to the timing

.. note::
    This information belongs in the SVN logfile.
"""

# Ensure python 2.5 compatibility
from __future__ import with_statement

import os
import sys
import time
import re
from pyraf import iraf
import saltsafekey
import saltsafeio
import salttime
import slottool
import pyfits
import numpy as np
from saltsafelog import logging 
from salterror import SaltError, SaltIOError

# Make sure the plotting functions work with an older version of matplotlib
try:
    import matplotlib.pyplot as plt
except ImportError:
    import matplotlib.pylab as plt

def slotutcfix(images,update,outfile,ampperccd,ignorexp,droplimit,inter,plotdata,logfile,verbose,debug):
    
    with logging(logfile,debug) as log:
        # set up the variables
        utc_list = []


        # is the input file specified?
        saltsafeio.filedefined('Input',images)

        # if the input file is a list, does it exist?
        if images[0] == '@':
            saltsafeio.listexists('Input',images)

        # parse list of input files and place them in order
        infiles=saltsafeio.listparse('Raw image',images,'','','')
        infiles.sort()

        # check input files exist
        saltsafeio.filesexist(infiles,'','r')

        # check to see if the output file exists and if so, clobber it
        if os.path.isfile(outfile):
            try:
                os.remove(outfile)
            except:
                raise SaltIOError('File ' + outfile + ' can not be removed')

        # open the outfile
        if outfile:
            try:
                fout=open(outfile,'w')
            except:
                raise SaltIOError('File ' + outfile + ' can not be opened')

        # get time of first exposure and basic information about the observations
        infile=infiles[0]
        struct=saltsafeio.openfits(infile)

        # check to make sure slotmode data
        detmode=saltsafekey.get('DETMODE',struct[0], infile)
        if detmode != 'Slot Mode':
            raise SaltIOError('Data are not Slot Mode Observations')

        # Check to see if SLOTUTCFIX has already been run
        # and print a warning if they have
        if saltsafekey.found('SLOTUTC', struct[0]):
            message='Data have already been processed by SLOTUTCFIX'
            log.warning(message)

        # check to make sure that it is the right version of the software
        scamver=saltsafekey.get('DETSWV', struct[0], infile)
        try:
            scamver=float(scamver.split('-')[-1])
            if 4.42 <= scamver <= 5.00:
                pass
            else:
                raise SaltError('cannot currently correct this version of the SCAM software.')

        except:
            raise SaltError('Not able to read software version')

        # requested exposure time
        req_texp=saltsafekey.get('EXPTIME',struct[0],infile)

        # how many extensions?
        nextend=saltsafekey.get('NEXTEND',struct[0],infile)

        # how many amplifiers
        amplifiers=saltsafekey.get('NCCDS',struct[0],infile)
        amplifiers = int(ampperccd*float(amplifiers))
        if ampperccd>0:
            nframes = nextend/amplifiers
            nstep=amplifiers
        else:
            nframes = nextend
            nstep=1

        # how many total frame and unique times
        ntotal=nextend*len(infiles)
        nunique=len(infiles)*nframes-ignorexp+1

        # Create arrays necessary for analysis
        id_arr=np.arange(nunique)
        utc_arr=np.zeros(nunique,dtype=float)

        # Read in each file and make a list of the UTC values
        if verbose:
            log.message('Reading in files to create list of UTC values.')

        j=0
        for n,infile in enumerate(infiles):
            # Show progress
            if verbose:
                percent=100.*float(n)/float(len(infiles))
                ctext='Percentage Complete: %.2f\r' % percent
                sys.stdout.write(ctext)
                sys.stdout.flush()

            struct=saltsafeio.openfits(infile)
            if not len(struct)-1==nextend:
                raise SaltIOError(infile,' has a different number of extensions from the first file')

            # Skip through the frames and read in the utc
            istart=1
            if infile==infiles[0]:
                istart=ignorexp*amplifiers+1
            for i in range(istart,len(struct), amplifiers):
                try:
                    utc_list.append(saltsafekey.get('UTC-OBS', struct[i], infile))
                    utc_arr[j]=slottool.getobstime(struct[i], infile)
                    j += 1
                except Exception, e:
                    raise SaltIOError('Unable to create array of UTC times.  Please check the number of extensions in the files')

            # close FITS file
            saltsafeio.closefits(struct)

        # set up the other important arrays
        try:
            diff_arr=utc_arr.copy()
            diff_arr[1:]=diff_arr[1:]-utc_arr[:-1]
            diff_arr[0]=-1
            dsec_arr=utc_arr-utc_arr.astype(int)
        except:
            raise SaltIOError('Unable to create timing arrays')

        # calculate the real exposure time
        if verbose:
            log.message('Calculating real exposure time.')

        real_expt, med_expt, t_start, t_arr, ysum_arr=calculate_realexptime(id_arr, utc_arr, dsec_arr, diff_arr, req_texp, utc_list)

        # plot the results
        if plotdata:
            if verbose:
                log.message('Plotting data.')
            plt.ion()
            plt.plot(t_arr,ysum_arr,linewidth=0.5,linestyle='-',marker='',color='b')
            plt.xlabel('Time (s)')
            plt.ylabel('Fit')

        # Calculate the corrrect values
        if verbose:
            log.message('Calculating correct values')

        i_start = abs(utc_arr-t_start).argmin()
        t_diff=utc_arr*0.0+real_expt
        nd=utc_arr*0.0
        ndrop=0
        for i in range(len(utc_arr)):
            if utc_arr[i] >= t_start:
                t_new=t_start+real_expt*(i-i_start+ndrop)
                t_diff[i]=utc_arr[i]-t_new
                while (t_diff[i]>real_expt and nd[i] < droplimit):
                    nd[i]+= 1
                    t_new=t_start+real_expt*(i-i_start+ndrop+nd[i])
                    t_diff[i]=utc_arr[i]-t_new
                if (nd[i]<droplimit):
                    ndrop += nd[i]
            else:
                t_new=t_start+real_expt*(i-i_start)
                t_diff[i]=utc_arr[i]-t_new
                while (t_diff[i]>real_expt and nd[i] < droplimit):
                    nd[i]+= 1
                    t_new=t_start+real_expt*(i-i_start-nd[i])
                    t_diff[i]=utc_arr[i]-t_new

        # calculate the corrected timestamp by counting 6 record files forward and
        # 8 recored + unrecorded files back--or just 8*t_exp forward.
        # if the object is near the end of the run, then just replace it with
        # the correct value assuming no dropped exposures.

        # first make the array of new times
        new_arr=utc_arr-t_diff

        # Next loop through them to find the corrected time
        corr_arr=utc_arr*0.0

        for i in range(len(new_arr)):
            if i+6 < len(new_arr)-1:
                corr_arr[i]=new_arr[i+6]-8*real_expt
            else:
                corr_arr[i]=new_arr[i]-2*real_expt

        t_diff=utc_arr-corr_arr

        # write out the first results
        msg="Dwell Time=%5.3f Requested Exposure Time=%5.3f Nobs = %i Dropped = %i" % (real_expt, req_texp,  nunique, ndrop)
        if verbose:
            log.message(msg)

        if outfile:
            fout.write('#'+msg+'\n')
            fout.write('#%23s %2s %12s %12s %10s %8s %4s \n' % ('File', 'N', 'UTC_old', 'UTC_new', 'UTC_new(s)', 'Diff', 'drop' ))

        # Give the user a chance to update the value
        if inter:
            message='Update headers with a dwell time of %5.3f s [y/n]? ' % real_expt
            update=saltsafeio.yn_ask(message)
            if not update:
                message='Set Dwell Time manually [y/n]? '
                update=saltsafeio.yn_ask(message)
                if update:
                    message='New Dwell Time: '
                    real_expt=saltsafeio.ask(message)
                    try:
                        real_expt=float(real_expt)
                    except Exception, e:
                        msg='Could not set user dwell time because %s' % e
                        raise SaltError(msg)

        # If requested, update the UTC times
        if update or outfile:
            if verbose:
                log.message('Updating UTC times')

            j=0
            for n,infile in enumerate(infiles):
                # Show progress
                if verbose:
                    percent=100.*float(n)/float(len(infiles))
                    ctext='Percentage Complete: %.2f\r' % percent
                    sys.stdout.write(ctext)
                    sys.stdout.flush()

                struct=saltsafeio.openupdatefits(infile)

                # Skip through the frames and read in the utc
                istart=1
                if infile==infiles[0]:
                    istart=ignorexp*amplifiers+1
                for i in range(istart,len(struct), amplifiers):
                    for k in range(0,amplifiers):
                        if update:
                            struct[i+k]=updateheaders(struct[i+k],i+k, t_diff[j], real_expt, utc_list[j], infile)

                    if outfile:
                        utc_new=saltsafekey.get('UTC-OBS', struct[i+k], infile)
                        utc_new_sec=slottool.getobstime(struct[i], infile)
                        fout.write('%25s %2i %12s %12s %7.3f  %5.4f %4i \n' % (infile, i, utc_list[j], utc_new, utc_new_sec, t_diff[j], nd[j] ))
                    j += 1

                # Housekeeping key words
                if update:
                    history = 'SALTUTCFIX -- '
                    history += 'images='+infile+' '
                    saltsafekey.housekeeping(struct[0],'SLOTUTC','UTC has been corrected',history,infile)

                # update fits file
                if update:
                    saltsafeio.updatefits(struct)
                    saltsafeio.closefits(struct)

        # close outfile
        if outfile:
            fout.close()

        # Keep plot window open
        if plotdata:
            plt.show()

def updateheaders(struct, ext, tdiff, real_expt, utc, infile):
    # exit if tdiff wasn't updated
    if tdiff == real_expt:
        msg='No adequate correction found for frame %i in file %s' % (ext, infile)
        raise SaltError(msg)

        return struct

    # calculate the new utc value
    try:
        ntime=salttime.sex2dec(utc)
        ntime=ntime-tdiff/3600.0
        newutc=salttime.dec2sex(ntime)
    except Exception,e:
        msg='Could not update UTC in %i header of image %s because %s' % (ext, infile, e)
        raise SaltError(msg)

        return struct

    # update the headers
    if utc==saltsafekey.get('UTC-OBS', struct):
        expt_string='%5.4f' % real_expt
        td_string='%5.4f' % tdiff
        if not saltsafekey.found('DUTC', struct):
            try:
                saltsafekey.put('UTC-OBS', newutc, struct, infile)
                saltsafekey.put('TIME-OBS', newutc, struct, infile)
                saltsafekey.new('DWETIME', expt_string, 'Dwell Time', struct, infile)
                saltsafekey.new('DUTC', td_string, 'Change in UTC time', struct, infile)
            except Exception, e:
                msg='Could not update %i header of image %s because %s' % (ext, infile, e)
                raise SaltIOError(msg)
        else:
            try:
                saltsafekey.put('UTC-OBS', newutc, struct, infile)
                saltsafekey.put('TIME-OBS', newutc, struct, infile)
                saltsafekey.put('DWETIME', real_expt, struct, infile)
                saltsafekey.put('DUTC', tdiff, struct, infile)
            except Exception, e:
                msg='Could not update %i header of image %s because %s' % (ext, infile, e)
                raise SaltError(msg)
    else:
        raise SaltIOError('Frame missing from list of times')

    return struct

def calculate_realexptime(id_arr, utc_arr, dsec_arr, diff_arr, req_texp, utc_list):
    """Calculates the real exposure time.
    This makes the following assumptions:
    #. That the measurement after the turn of the second is a fiducial
    #. That there is an integer number of frames between each fiducial exposure
    #. We then set up a metric which is Y=np.sum(i-int(i)) where i=dt/t_exp
    #. Then the minimum of Y is found between the requested exposure time and the median time difference
    #. And the best exposure time is the time at that minimum

    returns median exposure time and real exposure time
    """
    t_exp=0

    # calculate the median time
    try:
        t_wrong=np.median(diff_arr)
    except:
        raise SaltError('Unable to calculate median time difference')

    # Compress the arrays to find those closest to the second mark
    mask=(dsec_arr<t_wrong)*(diff_arr>0)
    t=np.compress(mask,utc_arr)
    s=np.compress(mask,dsec_arr)
    id=np.compress(mask,id_arr)

    # Now set up the components in the equation
    try:
        t_start=t[0]
        dt=t[1:]-t[0]
    except Exception, e:
        msg='Unable to set up necessary arrays because %s' % e
        raise SaltError(msg)

    # Now find the one that minimizes the equation y=np.sum(i-int(i))
    t_max=t_wrong*1.2
    if t_max < req_texp+0.05:
        t_max=req_texp+0.15
    try:
        t_arr, ysum_arr=find_real_time(dt, req_texp, t_max)
        t_real=t_arr[ysum_arr.argmin()]
    except Exception, e:
        msg='Unable to calculate real dwell time because %s' % e
        raise SaltError(msg)

    return t_real, t_wrong, t_start, t_arr, ysum_arr

def find_real_time(dt, t_min, t_max):
    """Find the real exposure+dead time

    returns float
    """
    t_e=np.arange(t_min,t_max,0.0001)
    ysum=t_e*0.0
    for i in range(len(t_e)):
        ysum[i]=ntime_func(dt, t_e[i])
    return t_e, ysum

def ntime_func(dt, t_e):
    """Merit function to determine best time
    Weighted for the number of objects in each step

    return float
    """
    y=0
    for j in range(len(dt)):
        i=dt[j]/t_e
        y += abs(i-round(i))
    return y

# -----------------------------------------------------------
# main code

parfile = iraf.osfn("slottools$slotutcfix.par")
t = iraf.IrafTaskFactory(taskname="slotutcfix",value=parfile,function=slotutcfix,pkgname='slottools')
