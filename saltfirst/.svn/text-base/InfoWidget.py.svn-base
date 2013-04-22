"""
InfoWidget is a Qt4 Widget for displaying information about an image
"""
from PyQt4 import QtGui, QtCore
from ObsLogWidget import headerList, printList

class InfoWidget(QtGui.QWidget):
   def __init__(self, name, imlist, parent=None):
        super(InfoWidget, self).__init__(parent)
        self.imlist=imlist
        #set up the information panel
        self.infopanel=QtGui.QWidget()

        #add the name of the file
        self.NameLabel = QtGui.QLabel("Filename:")
        self.NameLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised )
        self.NameValueLabel = QtGui.QLabel("%s" % name)
        self.NameValueLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken )

        #add target and proposal information

        #set up the info panel layout
        infoLayout=QtGui.QGridLayout(self.infopanel)
        infoLayout.addWidget(self.NameLabel, 0, 0, 1, 1)
        infoLayout.addWidget(self.NameValueLabel, 0, 1, 1, 1)

        #add all teh other fields
        self.ValueList=[]
        for i, k in enumerate(headerList[1:]):
            Label = QtGui.QLabel("%s:" % k)
            Label.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised )
            ValueLabel = QtGui.QLabel("%s" % self.getitem(k))
            ValueLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken )
            infoLayout.addWidget(Label, i+1, 0, 1, 1)
            infoLayout.addWidget(ValueLabel, i+1, 1, 1, 1)
            self.ValueList.append(ValueLabel)


        # Set up the layout
        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addWidget(self.infopanel)
        self.setLayout(mainLayout)


   def update(self, name, imlist):
       self.imlist=imlist
       self.NameValueLabel.setText(name)
       for i, k in enumerate(headerList[1:]):
           self.ValueList[i].setText("%s" % self.getitem(k))

   def getitem(self, key):
       i=headerList.index(key)
       try:
           value=str(self.imlist[i])
       except IndexError:
           value=''
       return value

