#!/usr/bin/env python

import sys
import os
import numpy as np

#import pylab as pl
#from scipy import optimize as opt

from ring import *

if len(sys.argv) != 3:
    print "Usage: calib <linewave><range>"
    print ""
    print "E.g.:  calib 6678.28 48-58"
    print "       calib 6678.28 48,59-69,75-83"

    exit()

line = sys.argv[1]
flist = sys.argv[2]

trim_rad = 470

try:
    list = []
    files = flist.split(',')
    for file in files:
        nums = file.split('-')
        if len(nums) == 1:
            list.append(int(nums[0]))
        else:
            for i in range(int(nums[0]), int(nums[1]) + 1):
                list.append(i)

except:
    print "Specify file range as <low>-<high>, e.g. 25-35."
    exit()

date = os.getcwd().split('/')[-1]

outfile = "calib_rings.dat"
out = open(outfile, "a")

for i in list:
    fits = "mbxpP%s%04d.fits" % (date, i)
    hdu = pyfits.open(fits)
    (data, header) = (hdu[1].data, hdu[0].header)
    dateobs = header['DATE-OBS']
    time = header['TIME-OBS']
    lamp = header['LAMPID']
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
    binning = int(header['CCDSUM'].split()[0])

    ysize, xsize = data.shape

    # cut FP image down to square
    fp_im = data[:, (xsize-ysize)/2:(xsize+ysize)/2]

    # mask those gaps
    fp_im = ma.masked_less_equal(data[:, (xsize-ysize)/2:(xsize+ysize)/2], 0.0)

    # define center based on FP ghost imaging with special mask
    xc = 2054 / binning
    yc = 2008 / binning
    # we use the 4x4 version of trim_rad since the vast majority of FP
    # data is taken with 4x4 binning
    trim_rad *= 4/binning

    mask_val = 0.0

    # get the radial profile and flatten it with a default QTH flat profile
    prof, r = FP_profile(fp_im, xc, yc, trim_rad=trim_rad, mask=mask_val)
    prof = flatprof(prof, binning)

    # find the peaks and bail out if none found
    npeaks, peak_list = find_peaks(prof, width=40)
    if npeaks > 0:
        cen_peak = centroid(prof, peak_list[0], 20)
        if np.isnan(cen_peak):
            cen_peak = peak_list[0]
        print "%s ring at R = %f" % (fits, cen_peak)
        out.write("%s %s %s %s %s %f %s\n" % (fits, dateobs, time, lamp, line, cen_peak, z))

out.close()
