#!/usr/bin/env python

import sys
import numpy as np
import pylab as pl

from scipy import optimize as opt

from ring import *

if len(sys.argv) != 2:
    print "Usage: parallel <X|Y>"
    print ""
    exit()

axis = sys.argv[1]

if axis.upper() == 'X':
    print "Doing X parallelism fit"
    col = 2
elif axis.upper() == 'Y':
    print "Doing Y parallelism fit"
    col = 3
else:
    print "Specify either X or Y."
    exit()

outfile = "%s_rings.dat" % axis.upper()

pl.figure()

x, rms = np.loadtxt(outfile, usecols=(col, 8), unpack=True)
init = [rms.max(), 0.2, x.mean()]
fit = opt.fmin_powell(focus_func, init, args=(rms, x),
                      ftol=0.00001, full_output=False, disp=False)
print ""
print "Min RMS = %.3f" % fit[0]
print "Best %s = %.2f" % (axis.upper(), fit[2])
print "Scale = %.3f" % fit[1]
pl.subplot(111)
pl.scatter(x, rms)
xp = range(int(x.min()), int(x.max()))
f = fit[0] - fit[1] * (xp - fit[2]) ** 2
pl.plot(xp, f)
pl.xlabel(axis.upper())
pl.show()
