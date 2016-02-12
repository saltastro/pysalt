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

"""saltmysql are salt tasks associated with connecting to and accessessing a mysql server database."""

import saltprint
import os, sys, string, types, re
import MySQLdb
import saltsafeio as saltio
from salttime import datatimeobs2DateTime,sex2dec
from salterror import SaltError
findNum = re.compile(r'(?:[^\d]*(\d+)[^\d]*)+')

class SALTMySQLError(SaltError):
    """Errors involving MYSQL should cause this exception to be raised."""
    pass


def connectdb(host,dbname,user,passwd):
    """Connect to  a database"""
    db=''
    try:
        db = MySQLdb.connect(host=host,db=dbname,user=user,passwd=passwd)
    except:
        message = 'Cannot connect to ' + host
        raise SALTMySQLError(message)

    return db

def connectelsview(host, dbname,user,passwd):
    """Connect to  a ELS viewer with converted to fix timestamp serialization issue
    """
    db = ''
    try:
       orig_conv = MySQLdb.converters.conversions
       conv_iter = iter(orig_conv)
       convert = dict(zip(conv_iter, [str,] * len(orig_conv.keys())))
       db = MySQLdb.connect(host=host,db=dbname,user=user,passwd=passwd, conv=convert)
    except Exception, e:
       message = 'Cannot connect to %s due to %s' % (host, str(e))
       raise SALTMySQLError(message)

    return db


def select(db, selection, table, logic):
    """Select a record from a table"""
    record=''
    exec_command    =""
    exec_command   +="SELECT "+selection
    exec_command   +=" FROM  "+table
    if len(logic)>0:
        exec_command   +=" WHERE  "+logic

    try:
        cursor = db.cursor()
        cursor.execute(exec_command)
        record = cursor.fetchall()
        cursor.close()

    except Exception, e:
        message = 'Cannot read value because %s in call:\n %s' % (e, exec_command)
        raise SALTMySQLError(message)

    return record

def insert(db, insertion, table):
    """Insert a record into the table"""
    exec_command    =""
    exec_command   +="INSERT INTO "+table
    exec_command   +=" SET  "+insertion

    try:
        cursor = db.cursor()
        cursor.execute(exec_command)
        cursor.execute("COMMIT")

    except Exception, inst:
        message = 'Cannot insert into %s because %s for call: \n %s' % (table, inst, exec_command)
        raise SALTMySQLError(message)

    return 

def update(db, insertion, table, logic):
    """update a record in a table"""
    exec_command    =""
    exec_command   +="UPDATE "+table
    exec_command   +=" SET  "+insertion
    if len(logic)>0:
        exec_command   +=" WHERE  "+logic


    try:
        cursor = db.cursor()
        cursor.execute(exec_command)
        cursor.execute("COMMIT")

    except Exception, inst:
        message = 'Cannot insert into ' + table + '\n'
        message = message + inst[1]
        message = 'Cannot update into %s because %s for call: \n %s' % (table, inst, exec_command)
        raise SALTMySQLError(message)

    return 

def delete(db, table,logic):
    """delete a record in a table"""
    try:
        exec_command    =""
        exec_command   +="DELETE FROM "+table
        if len(logic)>0:
            exec_command   +=" WHERE "+logic


        cursor = db.cursor()
        cursor.execute(exec_command)
        cursor.execute("COMMIT")

    except Exception, inst:
        message = 'Cannot delete  entry in ' + table + '\n'
        message = message + inst
        raise SALTMySQLError(message)

    return 

def createstring(hdr, key, informat, default=''):
   """Create a string from a header key given a format"""
   try:
       retstring=informat % hdr[key]
       return retstring.strip()
   except:
       return default

def getnightinfoid(db, date):
   """Get the night info id for a given date.  Create it if it does not exist"""
   logic="Date='%s'" % date
   record=select(db, 'NightInfo_Id', 'NightInfo', logic)
   if not record:
       #create new value and extract that value
       insert_command="Date='%s'" % date
       insert(db,insert_command,'NightInfo')
       record=select(db, 'NightInfo_Id', 'NightInfo', logic)

   if len(record)==1:
       return record[0][0]
   else:
       message="Too many entries in Nightinfo"
       raise SALTMySQLError(message)


def getobsstatsid(db, date):
   #get the nightinfo_id for this date
   nightinfoid=getnightinfoid(db, date)

   logic="NightInfo_Id=%i" % nightinfoid
   record=select(db, 'ObsStats_Id', 'ObsStats', logic)
   if not record:
       #create new value and extract that value
       insert_command="NightInfo_Id=%i" % nightinfoid
       insert(db,insert_command,'ObsStats')
       record=select(db, 'ObsStats_Id', 'ObsStats', logic)

   if len(record)==1:
       return record[0][0]
   else:
       message="Too many entries in ObsStats"
       raise SALTMySQLError(message)
      
   return 1

def getpropcodeid(db, propcode):
   """From a Proposal_Code, get the ProposalCode_Id"""
   
   logic="Proposal_Code='%s'" % propcode
   record=select(db, 'ProposalCode_Id', 'ProposalCode', logic)

   if not record:
       logic="Proposal_Code='None'" 
       record=select(db, 'ProposalCode_Id', 'ProposalCode', logic)

   if len(record)==1:
       return record[0][0]
   else:
       message="Too many entries in ProposalCode"
       raise SALTMySQLError(message)
      
   return 1

def getproposalcodes(obsdate, sdbhost,sdbname,sdbuser, sdbpassword):
    """Retrieve all the proposal observed on a given date"""
    db=connectdb(sdbhost,sdbname,sdbuser,sdbpassword)
    return getpropcodes(db, obsdate, clean=True)


def create_insert(db, ImageHeader, FileNameString, PipelineFileNameString):

    #set up all the image header variables
    ExposureTimeString=createstring(ImageHeader, 'EXPTIME', '%5.3f') 
    TelEpochString=createstring(ImageHeader, 'EQUINOX', '%6.2f', '2000.0')
    TargetNameString=createstring(ImageHeader, 'OBJECT', '%s') 
    InstrueNameString=createstring(ImageHeader, 'INSTRUME', '%s')
    ObsModeString=createstring(ImageHeader, 'OBSMODE', '%s')
    DetModeString=createstring(ImageHeader, 'DETMODE', '%s')
    InstrueNameString=createstring(ImageHeader, 'INSTRUME', '%s')
    ProposalCodeString=createstring(ImageHeader, 'PROPID', '%s', default='UNKNOWN')
    ObservatString=createstring(ImageHeader, 'OBSERVAT', '%s')
    ProposalString=createstring(ImageHeader, 'PROPOSAL', '%s')

    #get the ProposalCode_Id
    if not ProposalCodeString.strip():
       ProposalCodeString='NONE'
    ProposalCodeId=getpropcodeid(db, ProposalCodeString)

    #set the pipeline strings that are a little more difficult to set
    try:
        UTStartString=ImageHeader['DATE-OBS']+' '+ImageHeader['TIME-OBS']
    except:
        UTStartString='1000-01-01 00:00:00'
    if not UTStartString.strip():
       UTStartString='1000-01-01 00:00:00'

    try:
        TelRA=15.0*sex2dec(ImageHeader['TELRA'])
        TelDec=sex2dec(ImageHeader['TELDEC'])
        TelRAString='%10.7f' % TelRA
        TelDecString='%10.7f' % TelDec
    except Exception, e:  
        UTStartString='1000-01-01 00:00:00'
        TelRAString='0.00000'
        TelDecString='100.0000'
    try:
        PipelineStartString=datatimeobs2DateTime(ImageHeader['SAL-TLM'])
    except Exception, inst:
        PipelineStartString='1000-01-01 00:00:00'

    #set the number of exposures
    nexposures=1
    try:
        if ImageHeader['DETMODE'].count('SLOT'):
           nexposures=ImageHeader['NEXTEND']
           #this is hardwired but should be changed to use NAMPS
           if ImageHeader['INSTRUME'] is 'RSS':
              nexposures=nexposures/6
           else:
              nexposures=nexposures/4
    except:
        pass
     

    if not PipelineFileNameString:
           PipelineFileNameString=FileNameString
 
    try:
        name=os.path.basename(FileNameString).split('.')[0]
        date = '%s-%s-%s' % (name[1:5], name[5:7], name[7:9])
    except Exception, e:
        message='Could not determine the date for %s because %s' \
            % (FileNameString, e)
        raise SALTMySQLError(message)

    #get the obsstats
    #ObsStatsId=getobsstatsid(db, date)

    try:
       filesize=os.path.getsize(PipelineFileNameString)/1024.0
    except OSError:
       filesize=os.path.getsize(FileNameString)/1024.0
    except:
       filesize=0

    #get the block id
    try:
       BlockString=ImageHeader['BLOCKID'].strip()
       BlockString=saltio.checkfornone(BlockString)
       try:
          if int(BlockString)==0: BlockString=None
       except:
          pass
    except KeyError:
       BlockString=''
       
    #get the block id
    try:
       BlockVisitString=ImageHeader['BVISITID'].strip()
       BlockVisitString=saltio.checkfornone(BlockVisitString)
       try:
          if int(BlockVisitString)==0: BlockVisitString=None
       except:
          pass
    except KeyError:
       BlockVisitString=''



    #create the insertion command for the data
    try:
        insert_command='StepStats_Id=1,' 
        insert_command += "UTStart='"+UTStartString+"',"
        insert_command += "ProposalCode_Id=%i," % ProposalCodeId
        insert_command += "Target_Name='"+TargetNameString+"',"
        insert_command += "ExposureTime='"+ExposureTimeString+"',"
        insert_command += "TelRA='"+TelRAString+"',"
        insert_command += "TelDec='"+TelDecString+"',"
        insert_command += "TelEpoch='"+TelEpochString+"',"
        insert_command += "FileName='"+FileNameString+"',"
        insert_command += "PipelineDate='"+PipelineStartString+"',"
        insert_command += "PipelineFileName='"+PipelineFileNameString+"',"
        insert_command += "INSTRUME='"+InstrueNameString+"',"
        insert_command += "OBSMODE='"+ObsModeString+"',"
        insert_command += "DETMODE='"+DetModeString+"',"
        insert_command += "NExposures='%i',"%nexposures
        insert_command += "FileSize='%i'"%filesize
        if BlockString: insert_command += ",Block_Id='%s'"%BlockString
        if BlockVisitString: insert_command += ",BlockVisit_Id='%s'"%BlockVisitString
    except Exception, e:
        message='Could not create insert command because %s' % e
        raise SALTMySQLError(message)

    return insert_command

def createnewFileData(db, ImageHeader, FileNameString, PipelineFileNameString):
    """Create a new entry into the FileData table.

    If necessary, it will create k and put all data under there.
    """
    message=''
    except_message=''
    FileData_Id=-1


    insert_command=create_insert(db, ImageHeader, FileNameString, PipelineFileNameString)

    insert(db,insert_command,'FileData')

    #check to see if the FileData was created
    logic="FileName='%s'" % FileNameString
    records=select(db,'FileData_Id','FileData',logic)
    if len(records[0])==1 :
       FileData_Id=records[0][0]
       message='Created entry in FileData for ' + FileNameString
    else:
       message='New FileData entry cannot be accessed for %s ' % (inname)
       raise SALTMySQLError(message)

    return FileData_Id

def updateFileData(db, ImageHeader, FileData_Id, FileNameString, PipelineFileNameString):
    """Update an entry in the FileData table.

    If necessary, it will create k and put all data under there.
    """
    from salttime import datatimeobs2DateTime,sex2dec
    message=''
    except_message=''

    logic='FileData_Id=%i'%FileData_Id

    insert_command=create_insert(db, ImageHeader, FileNameString, PipelineFileNameString)

    update(db,insert_command,'FileData',logic)

    return 


def updateFitsHeaders(db, imageStruct, FileData_Id, log=None):
    """Update the fits header tables

    """
    message=''
    except_message=''
    logic='FileData_Id=%i'%FileData_Id

    try:
        ImageHeader=imageStruct[0].header
    except Exception, e:
        message='Cannot extract image header information because %s' %e
        raise SALTMySQLError(message)

    #get some of the basic information
    InstrueNameString=createstring(ImageHeader, 'INSTRUME', '%s')
    ProposalCodeString=createstring(ImageHeader, 'PROPID', '%s')
    ObservatString=createstring(ImageHeader, 'OBSERVAT', '%s')
    ProposalString=createstring(ImageHeader, 'PROPOSAL', '%s')


    #Check to see if FitsHeader exists, if not, then create
    records=select(db,'FileData_Id','FitsHeaderImage',logic)
    if not records:
       insert_command='FileData_Id=%i,' % FileData_Id
       insert_command += "PROPID='%s',"%ProposalCodeString
       insert_command += "PROPOSAL='%s',"%ProposalString
       insert_command += "OBSERVAT='%s'" % ObservatString
       insert(db,insert_command,'FitsHeaderImage')

    #update FistHeaderImage
    updatefitstable(db, 'FitsHeaderImage', ImageHeader, logic)

    #create the associate FitsHeaderPipeline table
    records=select(db,'FileData_Id','FitsHeaderPipeline',logic)
    if not records:
       insert_command='FileData_Id=%i' % FileData_Id
       insert(db,insert_command,'FitsHeaderPipeline')

    #update FitsHeaderPipeline Table
    updatefitstable(db, 'FitsHeaderPipeline', ImageHeader, logic)


    #update the fits headers, if it doesn't exist
    #create either the assoicate FitsHeaderRSS or SALTICAM tables
    insert_command='FileData_Id=%i' % FileData_Id
    if InstrueNameString=='SALTICAM':
       records=select(db,'FileData_Id','FitsHeaderSalticam',logic)
       if not records:
           insert(db,insert_command,'FitsHeaderSalticam')
       updatefitstable(db, 'FitsHeaderSalticam', ImageHeader, logic)
    elif InstrueNameString=='RSS':
       records=select(db,'FileData_Id','FitsHeaderRss',logic)
       if not records :
           insert(db,insert_command,'FitsHeaderRss')
       updatefitstable(db, 'FitsHeaderRss', ImageHeader, logic)
    elif InstrueNameString=='HRS':
       records=select(db,'FileData_Id','FitsHeaderHrs',logic)
       if not records :
           insert(db,insert_command,'FitsHeaderHrs')
       updatefitstable(db, 'FitsHeaderHrs', ImageHeader, logic)
     
    else:
       message='Did not update the Instrument tables'
       if log is not None:  
          log.warning(message)


    return 

def updatefitstable(db, Table, ImageHeader, logic):
    """Upate a fits header table"""
    try:
        exec_command='show columns from %s' % Table
        cursor = db.cursor()
        cursor.execute(exec_command)
        image_columns = cursor.fetchall()
        for column in image_columns:
            ColumnName=column[0]
            if not ColumnName=='FileData_Id':
                HeaderName=ColumnName.replace('_','-')
                if HeaderName=='DECL': HeaderName='DEC'
                if HeaderName=='PM-DEC': HeaderName='PM_DEC'
                try:
                   updateheaderinfo(db,Table,ImageHeader,ColumnName, HeaderName,column,logic)
                except KeyError:
                   pass 
                except SALTMySQLError, e:
                   pass    
    except Exception, e:
        message='Failed to update %s because %s' % (Table, e)
        raise SALTMySQLError(message)


    return 

def updateheaderinfo(db,Table,ImageHeader, ColumnName, HeaderName,column,logic):
    """Given a column, update it with header information"""

    HeaderString=ImageHeader[HeaderName]

    if type(HeaderString)==types.StringType:
        HeaderString=HeaderString.strip()
        try:
            maxChar=int(findNum.search(column[1]).groups()[0])
            if len(HeaderString) > maxChar:
                HeaderString=HeaderString[:maxChar]
            insertString="%s='%s'"% (ColumnName,HeaderString)
        except:
            insertString="%s='%s'"% (ColumnName,HeaderString)
    elif type(HeaderString)==types.FloatType:
        insertString='%s=%f'% (ColumnName,HeaderString)
    elif type(HeaderString)==types.IntType:
        insertString='%s=%i'% (ColumnName,HeaderString)
    else:
        #print HeaderString, ColumnName
        #print type(HeaderString)
        return

    update(db,insertString,Table,logic)

    return 

def getpiptusername(sdb, pid):
   #log into database and determine username from proposer name
   state_select='distinct p.Username'
   state_from='Proposal as pr join ProposalCode as c using  (ProposalCode_Id) join ProposalContact as pc using (Proposal_Id) join Investigator as i on (pc.Contact_Id=i.Investigator_Id) join PiptUser as p using (PiptUser_Id)'
   state_logic="pr.current=1 and  c.Proposal_Code='%s'" % pid
   record=select(sdb,state_select,state_from,state_logic)

   #check to see if it was successful and raise an error if not
   if len(record)<1:
       message='SALTFTP--Unable to find username for %s' % pid
       raise SaltError(message)

   #only returns the first one
   return record[0][0]

def getpropcodes(sdb, obsdate, clean=True):
    """Return all proposal codes for an observing date"""
    state_select='Distinct Proposal_Code'
    state_from='FileData join ProposalCode using (ProposalCode_Id)'
    state_logic="FileName like '%"+obsdate+"%'"
    record=select(sdb,state_select,state_from,state_logic)
    if len(record)==0:
       return []
    else:
       pids=[x[0] for x in record]
    if clean:
       pids=saltio.removebadpids(pids)
       pids=saltio.removeengineeringpids(pids)
    return pids


def getpiptemail(sdb, username):
   #log into database and determine username from proposer name
   state_select='Email'
   state_from='Investigator join PiptUser using (PiptUser_Id)'
   state_logic="Username='%s'" % username
   record=select(sdb,state_select,state_from,state_logic)

   #check to see if it was successful and raise an error if not
   if len(record)<1:
       message='SALTFTP--Unable to find username for %s' % pid
       raise SaltError(message)

   #only returns the first one
   return record[0][0]


