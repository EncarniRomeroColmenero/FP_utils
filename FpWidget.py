import numpy as np
import os, errno
from PyQt4 import QtGui,QtCore
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4 import NavigationToolbar2QT
from FpRingWidget import FpRingWidget
from FpParallWidget import FpParallWidget
from pyraf import iraf
import pyfits
from iraf import pysalt
from saltgui import MplCanvas
from erc_ring import fit_rings

class FpWidget (QtGui.QWidget):

    def __init__(self,filenumber,flatnumber=None,parent=None):
        super(FpWidget,self).__init__(parent)
        
        self.filenumber=filenumber

        if flatnumber:
            self.flatnumber=flatnumber        
        else: 
            self.flatnumber=0

        
        # add the tabs for Steve

        self.ringTab=FpRingWidget(self.filenumber,self.flatnumber)
        self.parallTab=FpParallWidget()

        self.fpTabWidget=QtGui.QTabWidget()
        self.fpTabWidget.addTab(self.ringTab,'FP Ring')
        self.fpTabWidget.addTab(self.parallTab,'FP parallelism')
        
       # Set up the main layout
        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addWidget(self.fpTabWidget)
        self.setLayout(mainLayout)

