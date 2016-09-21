#GUI stuff
import sys, argparse
from PyQt4 import QtGui

from pyraf import iraf
from iraf import pysalt
sys.path.insert(0, '/home/ccd/erc/FP_utils/')

from FpWidget import FpWidget


#------- MAIN ------

def main(ring,flat):

        ring=ring
        flat=flat


    #create GUI
        app=QtGui.QApplication([])        
#        aw=FpWidget(12,1)
        aw=FpWidget(ring,flat)
#        aw.setMinimumHeight(400)
#        aw.setMinimumWidth(1850)
        aw.show()



        sys.exit(app.exec_())

#-----------------------------------------------------------



if __name__ == "__main__": 


    parser = argparse.ArgumentParser(description='Argument list for fp.py')
    parser.add_argument('-flat','--flat', help='Flat file',required=False)
    parser.add_argument('-ring','--ring',help='Ring file', required=True)

    args = parser.parse_args()

    ring=int(args.ring)

    try:
            flat=int(args.flat)
    except:
            flat=0

    main(ring,flat)
