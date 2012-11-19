#!/usr/bin/env python

import sys
import pyfits
import numpy as np
from ring import *

file = sys.argv[1]

hdu = pyfits.open(file)
(data, header) = (hdu[-1].data, hdu[0].header)

binning = int(header['CCDSUM'].split()[0])
ysize, xsize = data.shape

# cut FP image down to square
fp_im = data[:, (xsize - ysize) / 2:(xsize + ysize) / 2]

# mask those gaps
fp_im = ma.masked_less_equal(
    data[:, (xsize - ysize) / 2:(xsize + ysize) / 2], 0.0)

# first guess is the center of the aperture (assume 4x4 binning here)
xc = 4 * 513 / binning
yc = 4 * 502 / binning

prof, r = FP_profile(fp_im, xc, yc, trim_rad=470, mask=0.0)

prof = prof / prof[0]

np.savetxt("flat.dat", np.transpose(prof))

pl.figure()
pl.plot(r, prof, label="Profile")
pl.legend(loc=2)
pl.show()
