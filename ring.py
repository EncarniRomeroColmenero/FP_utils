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
    # peaks is an array of indices of the peak locations.  
    peaks = np.nonzero( arr == arr_fil )[0]
    npeaks = len(peaks)
    flux_list = []
    for p in peaks:
        flux_list.append(arr[p])
    # sort peaks in descending order of flux
    ind = np.argsort(np.array(flux_list))[::-1]
    peak_list = peaks[ind]
    return npeaks, peak_list

# quick-n-dirty centroider
def centroid(data, pos, width, x=None):
    size = len(data)
    if x == None:
        x = np.arange(size)
    l = pos-width
    h = pos+width
    if h >= size-1: h = size-1
    dat = data[l:h]
    return np.sum(dat*x[l:h])/np.sum(dat)

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
cutsize = 3
cenwidth = 20

# here we find peaks in the 4 cardinal directions, centroid the brightest one, and use the
# results to fit for the ring center.  use image slices since median filtering doesn't appear
# to be mask-aware.  the gaps were creating spurious peaks.  
xup = fp_im[yc-cutsize:yc+cutsize,xc:].sum(axis=0)
n_xup, xup_list = find_peaks(xup)
xup_peak = centroid(xup, xup_list[0], cenwidth)
print "XUP:  R = %f" % xup_peak

xdown = fp_im[yc-cutsize:yc+cutsize,:xc].sum(axis=0)
n_xdown, xdown_list = find_peaks(xdown[::-1])
xdown_peak = centroid(xdown, xdown_list[0], cenwidth)
print "XDOWN:  R = %f" % xdown_peak
    
yup = fp_im[yc:,xc-cutsize:xc+cutsize].sum(axis=1)
n_yup, yup_list = find_peaks(yup)
yup_peak = centroid(yup, yup_list[0], cenwidth)
print "YUP:  R = %f" % yup_peak

ydown = fp_im[:yc,xc-cutsize:xc+cutsize].sum(axis=1)
n_ydown, ydown_list = find_peaks(ydown[::-1])
ydown_peak = centroid(ydown, ydown_list[0], cenwidth)
print "YDOWN:  R = %f" % ydown_peak

xc_new = 0.5*(xc+xup_peak + xc-xdown_peak)
yc_new = 0.5*(yc+yup_peak + yc-ydown_peak)

if np.sqrt( (xc-xc_new)**2 + (yc-yc_new)**2 ) > 25:
    print "Center moved too far."
else:
    xc = xc_new
    yc = yc_new
    
print xc, yc

prof, r = FP_profile(fp_im, xc, yc, Quadratic=False, trim_rad=ap_rad, mask=mask_val)
rsq = r**2

# first use max location as center, then centroid from there
max_r_i = prof.argmax()
max_r = centroid(prof, max_r_i, cenwidth, x=r)
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

