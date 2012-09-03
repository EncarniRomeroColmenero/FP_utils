#!/usr/bin/env python

import sys
import os
import numpy as np
import pylab as pl

from scipy import optimize as opt

from ring import *

if len(sys.argv) != 3:
    print "Usage: parallel <X|Y> <range>"
    print ""
    print "E.g.:  parallel X 48-58"
    print "       parallel Y 48,59-69,75-83"
    exit()

axis = sys.argv[1]
flist = sys.argv[2]

if axis.upper() == 'X':
    print "Doing X parallelism fit"
    col = 2
elif axis.upper() == 'Y':
    print "Doing Y parallelism fit"
    col = 3
else:
    print "Specify either X or Y."
    exit()

try:
    list = []
    files = flist.split(',')
    for file in files:
        nums = file.split('-')
        if len(nums) == 1:
            list.append(int(nums[0]))
        else:
            for i in range(int(nums[0]), int(nums[1])+1):
                list.append(i)

except:
    print "Specify file range as <low>-<high>, e.g. 25-35."
    exit()
    
date = os.getcwd().split('/')[-1]

outfile = "%s_rings.dat" % axis.upper()
out = open(outfile, "w")
pl.figure()

for i in list:
    fits = "mbxpP%s%04d.fits" % (date, i)
    print fits
    hdu = pyfits.open(fits)
    (data, header) = (hdu[1].data, hdu[0].header)
    etalon = int(header['ET-STATE'].split()[3])
    etwave_key = "ET%dWAVE0" % etalon
    name_key = "ET%dMODE" % etalon
    x_key = "ET%dX" % etalon
    y_key = "ET%dY" % etalon
    z_key = "ET%dZ" % etalon
    x = header[x_key]
    y = header[y_key]
    z = header[z_key]
    etname = header[name_key]
    cenwave = float(header[etwave_key])
    
    good, wave, prof, fit, pars = fit_rings(fits, flatfile="/home/ccd/FP_utils/sky_flat.dat")

    if good:
        resid = prof - fit
        rms = resid.std()
        max = prof.max()
        out.write("%s %s %5d %5d %5d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n" % (fits,
                                                                              etname,
                                                                              x,
                                                                              y,
                                                                              z,
                                                                              max,
                                                                              pars['R'][0],
                                                                              pars['Amplitude'][0],
                                                                              rms,
                                                                              pars['Gamma'][0],
                                                                              pars['FWHM'][0]))
    
        pfile = "%s_prof.dat" % fits
        np.savetxt(pfile, np.transpose((wave, prof, fit, resid)))

out.close()

x, peaks = np.loadtxt(outfile, usecols=(col, 7), unpack=True)
init = [peaks.max(), 0.2, x.mean()]
fit = opt.fmin_powell(focus_func, init, args=(peaks, x), ftol=0.00001, full_output=False, disp=False)
print ""
print "Peak flux = %.3f" % fit[0]
print "Best %s = %.2f" % (axis.upper(), fit[2])
print "Scale = %.3f" % fit[1]
pl.subplot(111)
pl.scatter(x, peaks)
xp = range(int(x.min()), int(x.max()))
f = fit[0] - fit[1]*(xp-fit[2])**2
pl.plot(xp, f)
pl.xlabel(axis.upper())
pl.show()

