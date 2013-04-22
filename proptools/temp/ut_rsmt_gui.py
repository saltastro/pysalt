from PyQt4 import QtGui, QtCore

from rsmt_gui import Ui_MainWindow
import sys


class ApplicationWindow(QtGui.QMainWindow):
  def __init__(self):
    # Setup widget
    QtGui.QMainWindow.__init__(self)

    # Set main widget
    self.main = QtGui.QWidget(self)

    self.ui=Ui_MainWindow()
    self.ui.setupUi(self)

    #vl = QtGui.QVBoxLayout()
    #vl.addWidget(self.ui)

    #self.main.setLayout(vl)

if __name__ == "__main__":
    #App = QtGui.QApplication([])
    f = Ui_MainWindow()
    print dir(f)
    #exit=App.exec_()

    App = QtGui.QApplication(sys.argv)
    aw = ApplicationWindow()
    aw.show()

    # Start application event loop
    exit=App.exec_()

