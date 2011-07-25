#!/usr/bin/env python

import sys
import os
import pyfits
import math
import numpy as np
import numpy.ma as ma
import scipy.ndimage as nd
import pylab as pl
from radialProfile import azimuthalAverage
from qvoigt import qvoigt
import ds9
from scipy import optimize as opt

# wrap calculating the radial profile to include options for using r or r**2 for the profile.
# want to display vs. r, but fit voigt profile vs. r**2.  also provide options to trim radius
# to exclude outside the FOV and to set a mask value.
def FP_profile(im, xcen, ycen, Quadratic=False, trim_rad=None, mask=0):
    prof = azimuthalAverage(im, center=[xcen,ycen], maskval=mask)
    r = np.arange(len(prof))
    if trim_rad:
        prof = prof[0:trim_rad]
        r = r[0:trim_rad]
    if Quadratic:
        return prof, r**2
    else:
        return prof, r

# wrapper to fit a voigt profile to given data.
def fit_func(p, data, rsq):
    pos = p[0]
    amp = p[1]
    fwhm = p[2]
    gam = p[3]
    back = p[4]
    return np.sum( (data - back - qvoigt(rsq,amp,pos,fwhm,gam))**2 )

# find peaks in a 1D array using the maximum filter technique given a specified filter width.
def find_peaks(arr, width=50):
    arr_fil = nd.filters.maximum_filter1d(arr, width)
    # original and filtered array are only equal at the peaks.
    # peaks is an array of indices of the peak locations
    peaks = np.nonzero( abs(arr-arr_fil) > 20.0 )[0]
    # look for transitions to define peaks
    pkdiff = np.diff(peaks)
    peak_ind = np.where(pkdiff > 1)[0]
    peak_list = []
    flux_list = []
    pk_start = 0
    npeaks = 0
    print peaks
    print pkdiff
    print peak_ind
    if len(peak_ind) < 3:
        for p in peak_ind:
            peak_list.append(peaks[p])
            flux_list.append(arr[peaks[p]])
            npeaks = npeaks + 1
    else:
        for p in peak_ind:
            p = p+1
            ave = ma.sum(arr[pk_start:p]*peaks[pk_start:p])/ma.sum(peaks[pk_start:p])
            peak_list.append(ave)
            flux_list.append(arr[int(ave)])
            npeaks = npeaks + 1
            pk_start = p
    # sort peaks in descending order of flux
    ind = np.argsort(np.array(flux_list))[::-1]
    peak_list = np.array(peak_list)[ind]
    return npeaks, peak_list

print "Opening ds9...."
disp = ds9.ds9()

f = pyfits.open(sys.argv[1])
hdu = f
(data, header) = (hdu[0].data, hdu[0].header)
f.close()

ysize, xsize = data.shape

# cut FP image down to square
fp_im = data[:,(xsize-ysize)/2:(xsize+ysize)/2]
disp.set_np2arr(fp_im, dtype=np.int32)
# mask those gaps
fp_im = ma.masked_less_equal(data[:,(xsize-ysize)/2:(xsize+ysize)/2], 0.0)

# first guess is the center of the aperture (assume 4x4 binning here)
xc=513
yc=502
ap_rad = 470.0
mask_val = 0.0

# find brightest ring and refine center
minsize = 10
pk_width = 30

#    xup = nd.filters.median_filter(fp_im[yc,xc:], size=minsize)
xup = fp_im[yc-3:yc+3,xc:].sum(axis=0)
n_xup, xup_list = find_peaks(xup)
print "XUP Nrings = %d" % n_xup
for peak in xup_list:
    print "\t R = %f" % peak

#    xdown = nd.filters.median_filter(fp_im[yc,:xc], size=minsize)
xdown = fp_im[yc-3:yc+3,:xc].sum(axis=0)
n_xdown, xdown_list = find_peaks(xdown[::-1])
print "XDOWN Nrings = %d" % n_xdown
for peak in xdown_list:
    print "\t R = %f" % peak
    
#    yup = nd.filters.median_filter(fp_im[yc:,xc], size=minsize)
yup = fp_im[yc:,xc-3:xc+3].sum(axis=1)
n_yup, yup_list = find_peaks(yup)
print "YUP Nrings = %d" % n_yup
for peak in yup_list:
    print "\t R = %f" % peak
    
#    ydown = nd.filters.median_filter(fp_im[:yc,xc], size=minsize)
ydown = fp_im[:yc,xc-3:xc+3].sum(axis=1)
n_ydown, ydown_list = find_peaks(ydown[::-1])
print "YDOWN Nrings = %d" % n_ydown
for peak in ydown_list:
    print "\t R = %f" % peak

xc_new = 0.5*(xc+xup_list[0] + xc-xdown_list[0])
yc_new = 0.5*(yc+yup_list[0] + yc-ydown_list[0])
if np.sqrt( (xc-xc_new)**2 + (yc-yc_new)**2 ) > 15:
    print "Center moved too far."
else:
    xc = xc_new
    yc = yc_new
    
print xc, yc

prof, r = FP_profile(fp_im, xc, yc, Quadratic=False, trim_rad=ap_rad, mask=mask_val)
rsq = r**2

# first use max location as center, then centroid from there
max_r_i = prof.argmax()
low_i = max_r_i-5
high_i = max_r_i+6
max_r = np.sum(r[low_i:high_i]*prof[low_i:high_i])/prof[low_i:high_i].sum()
pmax = prof[max_r_i]

print "Max is %f at R %f (%d)" % (pmax, max_r, max_r_i)
disp.set("regions command {circle %f %f %f # color=red}" % (xc, yc, max_r))

init = [max_r**2, pmax, 150.0, 10.0, 500.0]

fit = opt.fmin_powell(fit_func, init, args=(prof, rsq), xtol=0.0001, ftol=0.00001)

print "R = %.3f, Amp = %.3f, FWHM = %.3f, Gam = %.3f, Background = %.3f" % (np.sqrt(fit[0]), fit[1], np.sqrt(fit[2]), fit[3], fit[4])

fit_v = fit[4] + qvoigt(rsq, fit[1], fit[0], fit[2], fit[3])

pl.figure()
pl.plot(r, prof)
pl.plot(r, fit_v)
pl.plot([max_r, max_r], [0, pmax])
pl.show()

