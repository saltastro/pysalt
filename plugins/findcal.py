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
# DAMAGES (INCte: 2007/05/26
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)    #
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      #
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN #
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE          #
# POSSIBILITY OF SUCH DAMAGE.                                              #
############################################################################


#!/usr/bin/env python

"""
FINDCAL--From the list of observations for a night,
find the different settups that will be required

Author                 Version      Date
-----------------------------------------------
S M Crawford (SAAO)    0.1          19 Jun 2011

"""

from pyraf import iraf
from pyraf.iraf import pysalt
import saltsafemysql as saltmysql

import datetime


def findcal(obsdate, sdbhost, sdbname, sdbuser, password):
    """Find all of the unique configuration settings used on a night 
       and insert it into the database
  
       obsdate--Observation date in YYYYMMDD

       sdbhost--host name for sdb

       sdbname--name of the sdb

       sdbuser--user for the sdb

       password--sdb password
    """

    #connect to the db
    sdb=saltmysql.connectdb(sdbhost, sdbname, sdbuser, password)

    #get the nightinfo id
    logic="Date='%4s-%2s-%2s'" % (obsdate[0:4], obsdate[4:6], obsdate[6:8])
    results=saltmysql.select(sdb, 'NightInfo_Id', 'NightInfo', logic)
    night_id=results[0][0]

    #select all the scam data from this obsdate
    cmd_select='FileName, FileData_Id, CCDTYPE, DETMODE, CCDSUM, GAINSET, ROSPEED, FILTER'
    cmd_table=''' FileData 
  left join FitsHeaderImage using (FileData_Id) 
  left join FitsHeaderSalticam using (FileData_Id)
'''
    cmd_logic="FileName like 'S" + obsdate+"%' and CCDTYPE='OBJECT'"

    results=saltmysql.select(sdb, cmd_select, cmd_table, cmd_logic)

    #loop through all the results and return only the Set of identical results
    caldict=create_caldict(results)
               
    #insert the results into the database
    for k in caldict:
        #first check to see if it has already been entered
        record=saltmysql.select(sdb, 'FileData_Id', 'SalticamNightlyCalibration', 'FileData_Id=%i' % k)
        if len(record)<1:
           #check for block_id
           blockid=saltmysql.select(sdb, 'Block_Id', 'FileData', 'FileData_Id=%i' % k)[0][0]

           #get the calibration types requested
           if blockid:
               request=saltmysql.select(sdb, 'SalticamCalibrationType_Id', 'SalitcamCalibration', 'Block_Id=%i' % blockid)
               for cid in request:
                   cid=cid[0]
                   cmd_insert='NightInfo_Id=%i, FileData_Id=%i, SatlicamCalibrationType_Id=%i' % (night_id, k, cid)
                   saltmysql.insert(sdb, cmd_insert, 'SalitcamNightlyCalibration')
           print k, " ".join([str(k) for k in caldict[k]])


    #select all the RSS data from this obsdate
    rssheaderlist='CCDTYPE, DETMODE, OBSMODE, CCDSUM, GAINSET, ROSPEED, FILTER, GRATING, GR_STA, AR_STA, MASKID'
    cmd_select='FileName,FileData_Id, %s' % rssheaderlist
    # CCDTYPE, DETMODE, OBSMODE, CCDSUM, GAINSET, ROSPEED, FILTER, GRATING, GR_STA, AR_STA, MASKID'
    cmd_table=''' FileData 
  left join FitsHeaderImage using (FileData_Id) 
  left join FitsHeaderRss using (FileData_Id)
  join ProposalCode using (ProposalCode_Id)
'''
    cmd_logic="FileName like 'P" + obsdate+"%' and CCDTYPE='OBJECT' and Proposal_Code not like 'CAL_SPST'"

    results=saltmysql.select(sdb, cmd_select, cmd_table, cmd_logic)

    #loop through all the results and return only the Set of identical results
    caldict=create_caldict(results)
    #insert the scam results into the database
    for k in caldict:
       #first check to see if it has already been entered
       record=saltmysql.select(sdb, 'FileData_Id', 'RssNightlyCalibration', 'FileData_Id=%i' % k) 
       if len(record)<1:
           #period for checking for SPST.  If the uses requests this, it gets set to 
           # 7 days, but the default is taken within the last month for non-requests
           period=30

           #check for block_id
           blockid=saltmysql.select(sdb, 'Block_Id', 'FileData', 'FileData_Id=%i' % k)[0][0]

           #get the calibration types requested
           if blockid:
               request=saltmysql.select(sdb, 'RssCalibrationType_Id', 'RssCalibration', 'Block_Id=%i' % blockid)
               for cid in request:
                   cid=cid[0]
                   if cid==7: period=7
                   #Arcs and lamp flats are taken with the block and not separately
                   #spectroscopic standards are handled separately so that we can check the period for 
                   #them. 
                   if cid not in [1,2,3,7]:
                       cmd_insert='NightInfo_Id=%i, FileData_Id=%i, RssCalibrationType_Id=%i' % (night_id, k, cid)
                       saltmysql.insert(sdb, cmd_insert, 'RssNightlyCalibration')
                   #for biases, check to see if there is already an entry with the same settings
                   if cid==1:
                       if not checkforbias(sdb, k):
                          cmd_insert='NightInfo_Id=%i, FileData_Id=%i, RssCalibrationType_Id=%i' % (night_id, k, cid)
                          #saltmysql.insert(sdb, cmd_insert, 'RssNightlyCalibration')
 
                     
           if caldict[k][2]=='SPECTROSCOPY':

               #Add a check to make sure that something similar to this one isn't already in there
               #insert spec photometric standard
               if not checkforspst(sdb, k, caldict[k], rssheaderlist, period=period):
                   cmd_insert='NightInfo_Id=%i, FileData_Id=%i, RssCalibrationType_Id=7' % (night_id, k)
                   saltmysql.insert(sdb, cmd_insert, 'RssNightlyCalibration')

                #write proceedure to check for different requested calibration types

       print k," ".join([str(k) for k in caldict[k]])

    return

def checkforbias(sdb, k, instr='rss'):
    """Check to see if a bias of the same type already exists in the calibration table
    """ 
    if instr=='rss':
       caltable='RssNightlyCalibration' 
       logic='RssCalibrationType_Id=1'
       fitstable='FitsHeaderRss'
    elif instr=='scam':
       caltable='SalticamNightlyCalibration' 
       logic='SalticamCalibrationType_Id=1'
       fitstable='FitsHeaderSalticam'
 
 
    #get all the bias requests
    logic='CalibrationTaken is Null and Ignored=0 and %s' %  logic
    record=saltmysql.select(sdb, 'FileData_Id', caltable, logic)

    # if no results return false
    try: 
       record=record[0]
    except:
       return False
    
    select='CCDSUM,GAINSET, ROSPEED'
    table='FileData join %s using (FileData_Id)' % fitstable
    logic='FileData_Id=%i' % k
    kdata=saltmysql.select(sdb, select, table, logic)[0]
    for i in record:
        idata=saltmysql.select(sdb, select, table, 'FileData_Id=%i' % i)[0]
        if compare_configs(kdata, idata): return True
        
    return False


def checkforspst(sdb, fid, keylist, plist, period=7):
    """Check the sdb for SPST observations.  Returns true if the data have been taken in the 
       last period.

       fid: int
           FileData_Id

       dlist:  list
          header information for a given file

       keylist: list
          list of RSS header keywords

       period: int
          number of days in the past to check if data were taken
 
       Note
       ----
       The [:-1] selection in the comparison of the lists is to exclude maskid

    """
    try:
       utstart=saltmysql.select(sdb, 'UTStart', 'FileData', 'FileData_Id=%i' % fid)[0][0]
    except Exception, e:
       print e
       return False

    #set the period for to check for the data
    utstart=utstart-datetime.timedelta(days=period)

    #select all the RSS data from this obsdate
    cmd_select='FileName,FileData_Id, %s' % plist
    cmd_table=''' FileData 
  left join FitsHeaderImage using (FileData_Id) 
  left join FitsHeaderRss using (FileData_Id)
  join ProposalCode using (ProposalCode_Id)
'''
    cmd_logic="Proposal_Code like 'CAL_SPST' and UTSTART>'%s'" % (utstart)
    results=saltmysql.select(sdb, cmd_select, cmd_table, cmd_logic)

    #compare results
    for r in results: 
        if compare_configs(keylist[-1], r[2:-1]): 
           return True

    #now check the nightly calibrations table
    logic='CalibrationTaken is Null and Ignored=0 and RssCalibrationType_Id=7' 
    record=saltmysql.select(sdb, 'FileData_Id', caltable, logic)

    # if no results return false
    try: 
       record=record[0]
    except:
       return False

    #selection the data from the list and compare results
    cmd_select='FileName,FileData_Id, %s' % plist
    cmd_table=''' FileData 
  left join FitsHeaderImage using (FileData_Id) 
  left join FitsHeaderRss using (FileData_Id)
  join ProposalCode using (ProposalCode_Id)
'''
    for i in record:
        cmd_logic='FileData_Id=%i' % i
        results=saltmysql.select(sdb, cmd_select, cmd_table, cmd_logic)[0]
        if compare_configs(keylist[-1], r[2:-1]): 
           return True
        
    return False


    return False

        

def create_caldict(results):
    """Sort through the results from the query and create the calibration dictionary

    """
    #loop through all the results and return only the Set of identical results
    caldict={}
    for r in results:
        clist=r[2:]
        for k in caldict:
           if compare_configs(clist, caldict[k]):
                 clist=[]
        if clist: caldict[r[1]]=clist
    return caldict

def compare_configs(c1, c2):
    """Compare two configurations list to see if they are the same.  
       returns true if they are the same
    """
    for i in range(len(c1)):
        if c1[i]!=c2[i]: return False
    return True


if __name__=='__main__':
   import sys
   findcal(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
