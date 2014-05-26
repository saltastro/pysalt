################################# LICENSE ##################################
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.  See LICENSE file for more details                  #
#                                                                          #
############################################################################


#!/usr/bin/env python

"""
SDBLOADOBSLOG--Given an observing log, load it into 
teh SDB

Author                 Version      Date
-----------------------------------------------
S M Crawford (SAAO)    0.1          19 Jun 2011

"""
import saltsafemysql as saltmysql

def sdbloadobslog(logstr, obsdate, sdbhost, sdbname, sdbuser, password):
    """Load logstr into the SDB"""
    print logstr
    #connect to the db
    sdb=saltmysql.connectdb(sdbhost, sdbname, sdbuser, password)

    #get the nightinfoid
    #date='YYYY-MM-DD' <- note speechmarks
    #data= your table
    #Select NightInfo_Id from NightInfo where Date=$date;
    #(get nightinfo_Id)
    #Insert into ObsLogTable values ($nightinfo_Id, $data);
    logic="Date='%4s-%2s-%2s'" % (obsdate[0:4], obsdate[4:6], obsdate[6:8])
    results=saltmysql.select(sdb, 'NightInfo_Id', 'NightInfo', logic)
    night_id=results[0][0]

    #check to see if it exists
    logic='NightInfo_Id=%i' % night_id
    results=saltmysql.select(sdb, 'NightInfo_Id', 'ObsLogTable', logic)
    if results: 
        update_cmd="ObsLogTable='%s'" % (logstr)
        saltmysql.update(sdb, update_cmd, 'ObsLogTable', logic)
    else:
        insert_cmd="NightInfo_Id=%i, ObsLogTable='%s'" % (night_id, logstr)
        saltmysql.insert(sdb, insert_cmd, 'ObsLogTable')
    
    print obsdate, night_id


