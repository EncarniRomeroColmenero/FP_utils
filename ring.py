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
from qvoigt import qvoigt, voigt_fwhm
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

# wrapper to fit a set of voigt profiles to given data.
def fit_func(p, data, rsq):
    back = p[0]
    func = 0.0
    for i in range(1,len(p),4):
        pos = p[i]
        amp = p[i+1]
        fwhm = p[i+2]
        gam = p[i+3]
        func = func +  qvoigt(rsq,amp,pos,fwhm,gam)
    return np.sum( (data - back - func)**2 )

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

# here we find peaks in the 4 cardinal directions, centroid the brightest one, and use the
# results to fit for the ring center.  use image slices since median filtering doesn't appear
# to be mask-aware.  the gaps were creating spurious peaks.  
def find_center(im, xc, yc, cutsize=5, tolerance=25):
    xup = im[yc-cutsize:yc+cutsize,xc:].sum(axis=0)
    n_xup, xup_list = find_peaks(xup)
    xup_peak = centroid(xup, xup_list[0], cenwidth)

    xdown = im[yc-cutsize:yc+cutsize,:xc].sum(axis=0)
    n_xdown, xdown_list = find_peaks(xdown[::-1])
    xdown_peak = centroid(xdown, xdown_list[0], cenwidth)
    
    yup = im[yc:,xc-cutsize:xc+cutsize].sum(axis=1)
    n_yup, yup_list = find_peaks(yup)
    yup_peak = centroid(yup, yup_list[0], cenwidth)

    ydown = im[:yc,xc-cutsize:xc+cutsize].sum(axis=1)
    n_ydown, ydown_list = find_peaks(ydown[::-1])
    ydown_peak = centroid(ydown, ydown_list[0], cenwidth)

    xc_new = 0.5*(xc+xup_peak + xc-xdown_peak)
    yc_new = 0.5*(yc+yup_peak + yc-ydown_peak)

    if np.sqrt( (xc-xc_new)**2 + (yc-yc_new)**2 ) > tolerance:
        print "Center moved too far."
        return xc, yc
    else:
        return xc_new, yc_new

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

if __name__=='__main__':
    print "Opening ds9...."
    disp = ds9.ds9()

    f = pyfits.open(sys.argv[1])
    hdu = f
    (data, header) = (hdu[0].data, hdu[0].header)
    etalon = int(header['ET-STATE'].split()[3])
    etwave_key = "ET%dWAVE0" % etalon
    cenwave = float(header[etwave_key])
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
    ap_rad = 485.0
    mask_val = 0.0

    # convert r-squared to angstroms.
    lam_rsq = (cenwave/6563.0)*24.0/ap_rad**2

    # find brightest ring and refine center
    cutsize = 3
    cenwidth = 20

    prof, r = FP_profile(fp_im, xc, yc, Quadratic=False, trim_rad=500, mask=mask_val)
    rsq = cenwave - lam_rsq*r**2

    npeaks, peak_list = find_peaks(prof, width=40)

    print "Found %d rings at:" % npeaks
    for peak in peak_list:
        cen_peak = centroid(prof, peak, cenwidth)
        if np.isnan(cen_peak):
            cen_peak = peak

        print "\t R %f" % cen_peak
        disp.set("regions command {circle %f %f %f # color=red}" % (xc, yc, cen_peak))

    pl.figure()
 
    max_r = peak_list[0]
    pmax = prof[max_r]
    back = 2500.0
    fwhm = 1.0
    gam = 1.0
    init = [back]

    # keep 6 brightest
    if len(peak_list) > 6:
        peaks = peak_list[0:6]
    else:
        peaks = peak_list
    
    for peak in peaks:
        if peak > 30:
            # position
            init.append(cenwave-lam_rsq*peak**2)
            # amplitude
            init.append(prof[peak])
            # FWHM
            init.append(fwhm)
            # gamma
            init.append(gam)

    #fit = opt.fmin_slsqp(fit_func, init, args=(prof, rsq), bounds=bounds)
    #fit = opt.fmin_tnc(fit_func, init, args=(prof, rsq), bounds=bounds, approx_grad=True)
    fit = opt.fmin_powell(fit_func, init, args=(prof, rsq), ftol=0.00001)

    fit_v = fit[0]
    print "Background = %f" % fit_v
    for i in range(1,len(fit),4):
        fwhm = voigt_fwhm(fit[i+2], fit[i+3])
        print "\tR = %.3f, Amp = %.3f, Gauss FWHM = %.3f, Gamma = %.3f, FWHM = %.3f" % (fit[i], fit[i+1], fit[i+2], fit[i+3], fwhm)
        fit_v = fit_v + qvoigt(rsq, fit[i+1], fit[i], fit[i+2], fit[i+3])

    pl.plot(rsq, prof)
    pl.plot(rsq, fit_v)
#    pl.plot([max_r, max_r], [0, pmax])
    pl.show()

