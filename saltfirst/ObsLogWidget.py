"""
ObsLogwidget is a Qt4 Widget for displaying the night log for the night
"""
from PyQt4 import QtGui, QtCore

#headerList= ["Name", "UTC-OBS", "OBJECT", "PROPOSAL", "EXPTIME", "CCDSUM", "GAINSET", "ROSPEED", "OBSMODE", "DETMODE", "FILTER", "GRATING", "GR-ANGLE", "AR-ANGLE", "FOCUS", "MASKTYPE"]
headerList= ["Name", "TIME-OBS", "OBJECT", "PROPID", "EXPTIME", "CCDSUM", "GAINSET", "ROSPEED", "OBSMODE", "DETMODE", "CCDTYPE", "FILTER", "GRATING", "GR-ANGLE", "AR-ANGLE", "MASKID", "LAMPID", "FOCUS", "TELAZ", "TELALT", "SEEING", "NSOURCES", "BMEAN", "BMIDPT", "BSTD"]
printList=["Name", "OBJECT", "PROPID", "EXPTIME", "CCDSUM", "GAINSET", "ROSPEED", "OBSMODE", "DETMODE", "CCDTYPE", "FILTER"]


class ObsLogWidget(QtGui.QWidget):


   def __init__(self, obsdict=None, minrow=5, obsdate=None, parent=None):
        super(ObsLogWidget, self).__init__(parent)
        self.minrow=minrow
        self.obsdate=obsdate

        self.set_obsdict(obsdict)

        #Add a widget table
        self.obstable=QtGui.QTableWidget()
        self.obstable.setRowCount(self.nrow)
        self.obstable.setColumnCount(len(headerList))
        self.obstable.setHorizontalHeaderLabels(headerList)

        #set the rows

        for i in range(self.nrow):
           self.setrow(i, resize=False)  


        for j in range(0, len(headerList)):
           self.obstable.resizeColumnToContents(j)


        #set the selection and set the action
        if len(self.obsdict)>1:
           self.obstable.item(len(self.obsdict)-1, 0).setSelected(True)
        QtCore.QObject.connect(self.obstable, QtCore.SIGNAL('cellClicked (int,int)'), self.cellselected)
   

        #add a print button
        self.printButton = QtGui.QPushButton("Print")
        self.printButton.clicked.connect(self.printfornightlog)
	self.printButton.clicked.connect(self.printobslog)

        #add a calibration button
        self.calsButton = QtGui.QPushButton("Calibrations")
        self.calsButton.clicked.connect(self.findcals)
 

        # Set up the layout
        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addWidget(self.obstable)
        mainLayout.addWidget(self.calsButton)
        mainLayout.addWidget(self.printButton)
        self.setLayout(mainLayout)

   def cellselected(self, i, j):
       newfile=self.obstable.item(i,j).text()
       print i,j,newfile, 
       self.emit(QtCore.SIGNAL("cellclicked(QString)"), newfile)

   def set_obsdict(self, obsdict):
       if obsdict is None:
           self.obsdict={}
       else:
           self.obsdict=obsdict

       self.nrow=max(len(obsdict), self.minrow)

   def printobslog(self):
       """Print the observation log"""
       dpi=6
       print "Obslog", self.obstable.rowCount()
       hstr='#'
       #set up the output file
       if self.obsdate:
          fileout=self.obsdate+'.obslog'
          fout=open(fileout, 'w')
       for j in range(self.obstable.columnCount()):
           sform='%'+str(int(self.obstable.columnWidth(j)/dpi))+'s '
           hstr += sform % self.obstable.horizontalHeaderItem(j).text()
       #print hstr
       fout.write(hstr+'\n')
       for i in range(len(self.obsdict)):
           ostr=''
           for j in range(self.obstable.columnCount()):
               sform='%'+str(int(self.obstable.columnWidth(j)/dpi))+'s '
               ostr += sform % self.obstable.item(i, j).text()
           #print ostr
           fout.write(ostr+'\n')
       fout.close()

   def printfornightlog(self):
       """Print the obslog in the format for the night log"""
       hdrstr='%5s %8.8s %-20.20s %-15s %5s %3s %3s %3s %7s %7s %6s %6s %7s\n' \
              % ('#FILE', 'UT', 'OBJECT', 'PROPID', 'EXPT', 'BIN', 'DET', 'OBS', \
              'FILTER', 'GRATING', 'GR-ANG', 'AR-ANG', 'SLIT')
       outstr=''  
       for i in range(len(self.obsdict)):
           k=self.obsdict.order()[i]
           name=k[0]+k[9:13]
           timeobs=self.obsdict[k][headerList.index('TIME-OBS')].strip()
           objname=self.obsdict[k][headerList.index('OBJECT')].strip()
           propid=self.obsdict[k][headerList.index('PROPID')].strip()
           exptime=self.obsdict[k][headerList.index('EXPTIME')]
           ccdsum=self.obsdict[k][headerList.index('CCDSUM')]
           gainset=self.obsdict[k][headerList.index('GAINSET')].strip().upper()
           rospeed=self.obsdict[k][headerList.index('ROSPEED')].strip().upper()
           obsmode=self.obsdict[k][headerList.index('OBSMODE')].strip().upper()
           detmode=self.obsdict[k][headerList.index('DETMODE')].strip().upper()
           filtername=self.obsdict[k][headerList.index('FILTER')].strip()
           grating=self.obsdict[k][headerList.index('GRATING')].strip()
           graang=self.obsdict[k][headerList.index('GR-ANGLE')]
           arang=self.obsdict[k][headerList.index('AR-ANGLE')]
           slitname=self.obsdict[k][headerList.index('MASKID')]
           xbin,ybin=ccdsum.split()
           sbin='%sx%s' % (xbin, ybin)
           dmode=self.detset(gainset, rospeed)
           omode=self.obsset(obsmode, detmode)
           try:
              graang='%5.2f' % float(graang)
           except:
              graang=str(graang)
    
           try:
              arang='%5.2f' % float(arang)
           except:
              arang=str(arang)
       
           #36 char
           outstr +='%5s %8.8s %-20.20s %15s %5.2f %3s %3s %3s %7s %7s %6s %6s %7s\n' \
            % (name, timeobs, objname, propid, exptime, sbin, dmode, omode, \
               filtername, grating, graang, arang, slitname )
       logstr=hdrstr+outstr
       self.emit(QtCore.SIGNAL("updateobslogdb(QString)"), logstr)

   def findcals(self):
       self.emit(QtCore.SIGNAL("updatecals(QString)"), '')

   def obsset(self,obsmode, detmode):
       obs=''
       if obsmode=='IMAGING': obs='IM'
       if obsmode=='SPECTROSCOPY': obs='SP'
       if obsmode=='FABRY-PEROT': obs='FP'
       det=''
       if detmode=='NORMAL': det='N'
       if detmode=='SLOT-MODE': det='S'
       if detmode=='FRAME TRANSFER': det='F'
       if detmode=='DRIFT': det='D'
       return obs+det

   def detset(self,gainset, rospeed):
       g=''
       if gainset=='BRIGHT': g='B'
       if gainset=='FAINT': g='F'
       r=''
       if rospeed=='SLOW': r='S'
       if rospeed=='FAST': r='F'
       return g+r

   def setrow(self, i, resize=True):
       """Set all the values in a row from the obsdictionary"""
       if i >= len(self.obsdict): return
       k=self.obsdict.order()[i]
       nameItem=QtGui.QTableWidgetItem(k)
       self.obstable.setItem(i, 0, nameItem)
       self.obstable.resizeColumnToContents(0)
       for j in range(1, len(self.obsdict[k])):
           item=self.parseItem(self.obsdict[k][j])
           self.obstable.setItem(i, j, item)
           if resize: self.obstable.resizeColumnToContents(j)

   def addobsdict(self, name, obslist):
       """Add an item ot the obsdict"""
       self.obsdict[name]=obslist
       i=len(self.obsdict)
       if i>self.minrow:
          self.obstable.insertRow(i-1)
       self.setrow(i-1)


   def parseItem(self, x):
       """Parse an object so it can be entered into the table"""
       if isinstance(x, str):
           return QtGui.QTableWidgetItem(x)
       elif isinstance(x, float):
           return QtGui.QTableWidgetItem('%f' % x)
       elif isinstance(x, int):
           return QtGui.QTableWidgetItem('%i' % x)
       return QtGui.QTableWidgetItem('')
    



