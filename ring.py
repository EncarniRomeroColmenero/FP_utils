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


# wrap calculating the radial profile to include options for using r or r**2
# for the profile.  want to display vs. r, but fit voigt profile vs. r**2.
# also provide options to trim radius to exclude outside the FOV and to
# set a mask value.
def FP_profile(im, xcen, ycen, trim_rad=None, mask=0):
    prof = azimuthalAverage(im, center=[xcen, ycen], maskval=mask)
    r = np.arange(len(prof))
    if trim_rad:
        prof = prof[0:trim_rad]
        r = r[0:trim_rad]
    return prof, r


# wrapper to fit a set of voigt profiles to given data.
def fit_func(p, data, rsq):
    back = p[0]
    func = 0.0
    for i in range(1, len(p), 4):
        pos = p[i]
        amp = p[i+1]
        fwhm = p[i+2]
        gam = p[i+3]
        func = func + qvoigt(rsq, amp, pos, fwhm, gam)
    return np.sum((data - back - func) ** 2)


# parabolic function to fit for focus and parallelism
def focus_func(p, data, x):
    peak = p[0]
    a = p[1]
    b = p[2]
    func = peak - a * (x - b) ** 2
    return np.sum((data - func) ** 2)


# find peaks in a 1D array using the maximum filter technique
# given a specified filter width.
def find_peaks(arr, width=50):
    arr_fil = nd.filters.maximum_filter1d(arr, width)
    # original and filtered array are only equal at the peaks.
    # peaks is an array of indices of the peak locations.
    peaks = np.nonzero(arr == arr_fil)[0]
    npeaks = len(peaks)
    flux_list = []
    for p in peaks:
        flux_list.append(arr[p])
    # sort peaks in descending order of flux
    ind = np.argsort(np.array(flux_list))[::-1]
    peak_list = peaks[ind]
    p = []
    for peak in peak_list:
        if peak > 30:
            p.append(peak)
    npeaks = len(p)
    return npeaks, p


# here we find peaks in the 4 cardinal directions, centroid the brightest one,
# and use the results to fit for the ring center.  use image slices since
# median filtering doesn't appear to be mask-aware.  the gaps were creating
# spurious peaks.
def find_center(im, xc, yc, cutsize=5, tolerance=25):
    xup = im[yc-cutsize:yc+cutsize, xc:].sum(axis=0)
    n_xup, xup_list = find_peaks(xup)
    xup_peak = centroid(xup, xup_list[0], cenwidth)

    xdown = im[yc-cutsize:yc+cutsize, :xc].sum(axis=0)
    n_xdown, xdown_list = find_peaks(xdown[::-1])
    xdown_peak = centroid(xdown, xdown_list[0], cenwidth)

    yup = im[yc:, xc-cutsize:xc+cutsize].sum(axis=1)
    n_yup, yup_list = find_peaks(yup)
    yup_peak = centroid(yup, yup_list[0], cenwidth)

    ydown = im[:yc, xc-cutsize:xc+cutsize].sum(axis=1)
    n_ydown, ydown_list = find_peaks(ydown[::-1])
    ydown_peak = centroid(ydown, ydown_list[0], cenwidth)

    xc_new = 0.5*(xc+xup_peak + xc-xdown_peak)
    yc_new = 0.5*(yc+yup_peak + yc-ydown_peak)

    if np.sqrt((xc-xc_new)**2 + (yc-yc_new)**2) > tolerance:
        print "Center moved too far."
        return xc, yc
    else:
        return xc_new, yc_new


# quick-n-dirty centroider
def centroid(data, pos, width, x=None):
    size = len(data)
    if x is None:
        x = np.arange(size)
    l = pos-width
    h = pos+width
    if h >= size-1:
        h = size-1
    dat = data[l:h]
    return np.sum(dat*x[l:h])/np.sum(dat)

# flat-field a profile using analytic fit
def flatprof(profile, binning):
    apix = 0.1267*binning
    r = np.arange(len(profile))*apix
    a = 1.00475
    b = -0.00027789
    c = -1.26563e-07
    d = 9.76049e-10
    e = -2.22584e-12
    flat = a + b*r + c*r**3 + d*r**4 + e*r**5
    return profile/flat

def fit_rings(file, trim_rad=470, disp=None):
    hdu = pyfits.open(file)
    (data, header) = (hdu[1].data, hdu[0].header)
    etalon = int(header['ET-STATE'].split()[3])
    etwave_key = "ET%dWAVE0" % etalon
    name_key = "ET%dMODE" % etalon
    etname = header[name_key]
    cenwave = float(header[etwave_key])
    binning = int(header['CCDSUM'].split()[0])
    if header['OBSMODE'] != 'FABRY-PEROT':
        return False, np.empty(1), np.empty(1), np.empty(1), np.empty(1)

    ysize, xsize = data.shape

    # cut FP image down to square
    fp_im = data[:, (xsize-ysize)/2:(xsize+ysize)/2]
    if disp:
        disp.set_np2arr(fp_im, dtype=np.int32)
    # mask those gaps
    fp_im = ma.masked_less_equal(data[:, (xsize-ysize)/2:(xsize+ysize)/2], 0.0)

    # define center based on FP ghost imaging with special mask
    xc = 2052 / binning
    yc = 2008 / binning
    mask_val = 0.0
    f = {}
    if cenwave < 5200:
        f['MR'] = 22149.6
        f['LR'] = 22795.92
    if cenwave > 6500 and cenwave < 6600:
        f['MR'] = 22713.0
        f['LR'] = 24191.40
    if cenwave >= 6600 and cenwave < 6700:
        f['MR'] = 22848.0
        f['LR'] = 24169.32
    if cenwave >= 6700 and cenwave < 6900:
        f['MR'] = 22828.92
        f['LR'] = 24087.68
    else:
        f['MR'] = 22828.92
        f['LR'] = 24400.32

    # find brightest ring and refine center
    cutsize = 3
    cenwidth = 20

    prof, r = FP_profile(fp_im, xc, yc, trim_rad=trim_rad, mask=mask_val)
    prof = flatprof(prof, binning)

    wave = cenwave/np.sqrt(1.0 + (r*binning/f[etname])**2)

    npeaks, peak_list = find_peaks(prof, width=40)
    if npeaks < 1:
        print "No peaks found."
        return False, np.empty(1), np.empty(1), np.empty(1), np.empty(1)

    print "Found %d rings at:" % npeaks
    for peak in peak_list:
        cen_peak = centroid(prof, peak, cenwidth)
        if np.isnan(cen_peak):
            cen_peak = peak

        print "\t R %f" % cen_peak
        if disp:
            disp.set("regions command {circle %f %f %f # color=red}" %
                     (xc, yc, cen_peak))

    max_r = peak_list[0]
    pmax = prof[max_r]
    back = 100.0
    fwhm = 6.0
    gam = 2.0
    init = [back]
    bounds = [(0.0, 2.0e8)]

    # keep 6 brightest
    if len(peak_list) > 6:
        peaks = peak_list[0:6]
    else:
        peaks = peak_list

    for peak in peaks:
        # position
        init.append(cenwave/np.sqrt(1.0 + (peak*binning/f[etname])**2))
        bounds.append((cenwave-30, cenwave+30))
        # amplitude
        init.append(prof[peak])
        bounds.append((0.0, 1.0e8))
        # FWHM
        init.append(fwhm)
        bounds.append((0.1, 30.0))
        # gamma
        init.append(gam)
        bounds.append((0.0, 5.0))

    #fit = opt.fmin_slsqp(fit_func, init, args=(prof, wave),
    #                     bounds=bounds)
    fit, nfeval, rc = opt.fmin_tnc(fit_func, init, args=(prof, wave), epsilon=0.0001,
                                   bounds=bounds, approx_grad=True, disp=0, maxfun=1000)
    #fit = opt.fmin_powell(fit_func, init, args=(prof, wave),
    #                      ftol=0.00001, full_output=False, disp=False)

    print "Return code: %s" % opt.tnc.RCSTRINGS[rc]

    pars = {}
    fit_v = fit[0]
    print "\nBackground = %f" % fit[0]
    pars['Background'] = fit[0]
    pars['R'] = []
    pars['Amplitude'] = []
    pars['Gauss FWHM'] = []
    pars['Gamma'] = []
    pars['FWHM'] = []
    for i in range(1, len(fit), 4):
        fwhm = voigt_fwhm(fit[i+2], fit[i+3])
        pars['R'].append(fit[i])
        pars['Amplitude'].append(fit[i+1])
        pars['Gauss FWHM'].append(fit[i+2])
        pars['Gamma'].append(fit[i+3])
        pars['FWHM'].append(fwhm)
        fit_v = fit_v + qvoigt(wave, fit[i+1], fit[i], fit[i+2], fit[i+3])

    return True, wave, prof, fit_v, pars

if __name__ == '__main__':
    filenum = sys.argv[1]

    try:
        filenum = int(filenum)
    except:
        print "Specify file as a number, e.g. 25"
        exit()

    date = os.getcwd().split('/')[-1]

    fits = "mbxpP%s%04d.fits" % (date, filenum)

    good, rsq, prof, fit, pars = fit_rings(fits)
    resid = prof - fit
    rms = resid.std()
    print "\tR = %.3f, Amp = %.3f, RMS = %.3f, Gamma = %.3f, FWHM = %.3f" % \
          (pars['R'][0],
           pars['Amplitude'][0],
           rms,
           pars['Gamma'][0],
           pars['FWHM'][0])
    pl.figure()
    pl.subplot(211)
    pl.plot(rsq, prof, label="Profile")
    pl.plot(rsq, fit, label="Fit")
    pl.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
              ncol=2, mode="expand", borderaxespad=0.)
    pl.subplot(212)
    pl.plot(rsq, resid, label="Profile - Fit")
    pl.legend(loc=1)
    pl.show()
