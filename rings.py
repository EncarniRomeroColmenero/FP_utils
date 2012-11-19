#!/usr/bin/env python

import sys
import glob
import numpy as np
import pylab as pl

from ring import *

files = glob.glob(sys.argv[1])
out = open("rings.dat", "w")
pl.figure()

for file in files:
    print file
    hdu = pyfits.open(file)
    (data, header) = (hdu[0].data, hdu[0].header)
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

    good, wave, prof, fit, pars = fit_rings(file, flatfile=None)
    if good:
        resid = prof - fit
        rms = resid.std()
        max = prof.max()
        out.write(
            "%s %s %5d %5d %5d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n" %
            (file,
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

        pfile = "%s_prof.dat" % file
        np.savetxt(pfile, np.transpose((wave, prof, fit, resid)))
