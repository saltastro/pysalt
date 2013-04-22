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


