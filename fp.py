#GUI stuff
import sys
from PyQt4 import QtGui
from FpWidget import FpWidget


#------- MAIN ------

def main():

    #create GUI
        app=QtGui.QApplication([])        
#        aw=FpWidget(12,1)
        aw=FpWidget(13)
#        aw.setMinimumHeight(400)
#        aw.setMinimumWidth(1850)
        aw.show()



	sys.exit(app.exec_())

#-----------------------------------------------------------



if __name__ == "__main__": 

    main()
