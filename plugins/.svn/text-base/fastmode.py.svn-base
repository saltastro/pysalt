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
FASTMODE -- Provide real time response for each program.  First, this checks to see 
if the proposal asks for fasttime response.  If so, the task then copies the data
to the fast directory for that proposal on saltpipe.   If it is the first object
data for that proposal, then it will also send an email to the contact PI

Author                 Version      Date
-----------------------------------------------
S M Crawford (SAAO)    0.1          18 Jan 2012

"""
import os
import saltsafemysql as saltmysql
import saltsafestring as saltstring
from saltsafeio import email as sendemail


def runfast(filename, propcode, obsdate, server, readmefile, sdbhost, sdbname, sdbuser, password):
   """Handle fast data delivery for the proposal.  For a given filename
   """
   
   if propcode is None or propcode=='None': return

   #first check in the sdb if fast data delivery is needed
   sdb=saltmysql.connectdb(sdbhost,sdbname, sdbuser, password)
   select_term='Distinct Surname, email, username, ProposalCode_Id'
   from_term='''
Block join Pointing using (Block_Id) 
      join PipelineConfig using (Pointing_Id)
      join Proposal using (Proposal_Id)
      join ProposalCode using (ProposalCode_Id)
      join PipelineDataAccessMethod using (PipelineDataAccessMethod_Id)
      join ProposalContact using (Proposal_Id)
      join Investigator on (Investigator_Id=Contact_Id)
      join PiptUser using (PiptUser_Id)
   '''
   where_term="Proposal_Code like '%s' and current=1 and DataAccessMethod='Fast'" \
              % (propcode)
   print 'Select %s from %s where %s' % (select_term, from_term, where_term)
   try:
       record=saltmysql.select(sdb, select_term, from_term, where_term)
   except Exception, e:
       print e
       return None
   print "Checking for fast data"
   print record
 
   if record:
      surname, email, username, propid= record[0]
      print surname, email, username, propid
   else:
      return

   #second if so, then copy the data to the contact PI directory
   #on saltpipe under the fast directory. 
   #rawfilename=getrawfilename(filename)
   y=os.system('scp %s sa@saltpipe:/salt/ftparea/%s/fast%s/' % (filename, username, obsdate))
   if y==256:
      y=os.system('ssh sa@saltpipe mkdir /salt/ftparea/%s/fast%s' % (username, obsdate))
      y=os.system('scp %s sa@saltpipe:/salt/ftparea/%s/fast%s/' % (filename, username, obsdate))
  
   if y!=0:
      print "Problem with copying file %s to /salt/ftparea/%s/fast%s/"  % (filename,  username, obsdate)

   #copy the reduced data
   y=os.system('scp mbxp%s sa@saltpipe:/salt/ftparea/%s/fast%s/' % (os.path.basename(filename), username, obsdate))
   #check the type of data it is and copy over an ancillery data as well

   #if it is the first object file, check to see if an email has been
   #sent, and if not, send email
   #try to copy the spectroscopic data
   print filename, filename.startswith('P')
   if os.path.basename(filename).startswith('P'):
       sfilename='smbxp%s.txt' % (os.path.basename(filename).split('.fits')[0])
       print sfilename
       try:
           y=os.system('scp %s sa@saltpipe:/salt/ftparea/%s/fast%s/' % (sfilename, username, obsdate))
       except Exception, e:
           print e
   if os.path.basename(filename).startswith('S'):
       try:
           sfilename='mbxp%s.cat' % (os.path.basename(filename).split('.fits')[0])
           print sfilename
           y=os.system('scp %s sa@saltpipe:/salt/ftparea/%s/fast%s/' % (sfilename, username, obsdate))
       except Exception, e:
           print e
 

   #check to see if an email has been sent
   select_term='PipelineStatus'
   from_term='''
PipelineProposalStatistics 
join PipelineStatus using (PipelineStatus_Id)
join NightInfo using (NightInfo_Id)
join ProposalCode using (ProposalCode_Id)
'''
   where_term="Proposal_Code like '%s' and Date='%s-%s-%s'" % (propcode, obsdate[0:4], obsdate[4:6], obsdate[6:8])
   print select_term, from_term, where_term
   try:
     record=saltmysql.select(sdb, select_term, from_term, where_term)[0][0]
   except:
     record=None
   print record

   if record=='FastEmail':
      return
   else:
      #insert information into the database
      nightinfoid=saltmysql.getnightinfoid(sdb, obsdate)
      insert_term="NightInfo_Id=%i, ProposalCode_Id=%i, PipelineStatus_Id=8" % (nightinfoid, propid)
      table_term="PipelineProposalStatistics"
      saltmysql.insert(sdb, insert_term, "PipelineProposalStatistics")

      #send email
      sender='sa@salt.ac.za'
      recipient=email
      bcc='crawfordsm@gmail.com'
      subject='SALT data available for %s' % propcode
      message=open(readmefile).read()
      message=message.replace('OBSDATE', obsdate)
      sendemail(server,'sa',password,sender,recipient,bcc, subject,message)
      
   

   sdb.close()
   return

def getrawfilepath(filename):
    """Given a raw file name, returns the path on the SALT server of the raw file"""
    if filename.count('S'):
       ddir='salt/scam/'
       i=filename.index('S')
       filedate=saltstring.filedate(filename[i:])
    elif filename.count('P'):
       ddir='salt/rss/'
       i=filename.index('P')
       filedate=saltstring.filedate(filename[i:])
    print ddir, filedate
    return '%s%s/%s/raw/%s' % (ddir, filedate[0:4], filedate[4:8], filename)
