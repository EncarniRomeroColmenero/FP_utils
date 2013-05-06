#!/usr/bin/env python

import sys
import os
import ring
import pyfits

import numpy as np
import pylab as pl
import pandas as pd
import numpy.ma as ma

from scipy import optimize as opt


def calib_func(p, data, line):
    """
    likelihood function used to fit wavelength as a function of gap, Z, and
    ring radius, R.  pass data as list-like (R, Z).

    Parameters
    ----------
    P : array-like containing fit parameters
    data: tuple of numpy arrays containing radii and Z values, (R, Z)
    line: wavelengths of the calibration lines

    Returns
    -------
    likelihood as float
    """
    a = p[0]
    b = p[1]
    f = p[2]
    r = data[0]
    z = data[1]
    if len(p) == 4:
        t = p[3]
        td = data[2]
        func = (a + t * td + b * z) / np.sqrt(1.0 + (r / f) ** 2)
    else:
        func = (a + b * z) / np.sqrt(1.0 + (r / f) ** 2)
    return np.sum((line - func) ** 2)


def calibrate(l, r, z, time):
    """
    function to take lambda, r, and Z data, fit calib_func to them,
    and return/plot the results

    Parameters
    ----------
    l : numarray-like containing wavelengths of calibration lines
    r : numarray-line containing ring radii
    z : numarray-like containing Z values
    time : numarray-like containing timedelta values in hours

    Returns
    -------
    results as dict
    """
    if np.abs(time).max() > 0.25:
        init = [6680.0, 0.2, 21000.0, 0.0]
        fit = opt.fmin_powell(calib_func, init,
                              args=((r, z, time), l), ftol=0.00001)
        results = dict(zip(['A', 'B', 'F', 'T'], fit))
    else:
        init = [6680.0, 0.2, 21000.0]
        fit = opt.fmin_powell(calib_func, init,
                              args=((r, z), l), ftol=0.00001)
        results = dict(zip(['A', 'B', 'F'], fit))
        results['T'] = 0.0
    return results


# for running from the command line which we do for now in normal operation
if __name__ == '__main__':

    if len(sys.argv) != 2:
        print "Usage: calib <file range>"
        print ""
        print "E.g.:  calib 48-58"
        print "       calib 48,59-69,75-83"

        exit()

    flist = sys.argv[1]

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
        print "Specify file range as <low>-<high>, e.g. 25-35 or 23,25-35,40."
        exit()

    date = os.getcwd().split('/')[-1]

    outfile = "calib_rings_%s.dat" % flist

    if not os.path.exists(outfile):
        out = open(outfile, "w")
        out.write("# time  lamp  line  R  Z\n")

        for i in list:
            fits = "mbxgpP%s%04d.fits" % (date, i)
            if not os.path.exists(fits):
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

            lines = ring.linelist[lamp]
            line = lines[(lines > cenwave - 25) & (lines < cenwave)][0]

            ysize, xsize = data.shape

            # cut FP image down to square
            fp_im = data[:,
                         (xsize - ysize) / 2:
                         (xsize + ysize) / 2]

            # mask those gaps
            fp_im = ma.masked_less_equal(data[:,
                                              (xsize - ysize) / 2:
                                              (xsize + ysize) / 2], 0.0)

            # define center based on FP ghost imaging with special mask
            xc = 2054 / binning
            yc = 2008 / binning
            # we use the 4x4 version of trim_rad since the vast majority of FP
            # data is taken with 4x4 binning
            trim_rad *= 4 / binning

            mask_val = 0.0

            # get the radial profile
            prof, r = ring.FP_profile(fp_im, xc, yc,
                                      trim_rad=trim_rad, mask=mask_val)

            # find the peaks and bail out if none found
            npeaks, peak_list = ring.find_peaks(prof, width=40)
            if npeaks > 0:
                cen_peak = ring.centroid(prof, peak_list[0], 20)
                if np.isnan(cen_peak):
                    cen_peak = peak_list[0]
                print "%s ring at R = %f" % (fits, cen_peak)
                if cen_peak < 452.0:
                    out.write("%s %s %s %s %f %s\n" % (dateobs, time, lamp,
                                                       line, cen_peak, z))

        out.close()

    # now analyze the data
    data = pd.read_table(outfile, sep='\s*', parse_dates=[[0, 1]])

    l = data.line
    r = 4.0 * data.R
    z = data.Z
    time = data['#_time']
    td = time - time[0]
    # convert time deltas to floats and from ns to hours
    td = td.astype(np.float) / (3600 * 1.0e9)
    results = calibrate(l, r, z, td)
    for k, v in results.items():
        print "%s : %f" % (k, v)

    a = results['A']
    b = results['B']
    f = results['F']
    t = results['T']
    func = (a + t * td + b * z) / np.sqrt(1.0 + (r / f) ** 2)
    resid = func - l
    print "RMS: %.3f" % resid.std()
    pl.scatter(z, resid)
    pl.xlabel("Z")
    pl.ylabel("$\delta\lambda$ ($\AA$)", fontsize=16)
    pl.title("FP Calibration Residuals")
    pl.show()
