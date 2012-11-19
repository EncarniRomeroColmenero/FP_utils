import numpy as np
import numpy.ma as ma


def azimuthalAverage(image, center=None, maskval=0):
    """
    calculate the azimuthally averaged radial profile.

    image - 2D image
    center - [x,y] pixel coordinates used as the center. the default is
             None which then uses the center of the image (including
             fractional pixels).
    maskval - threshold value for including data in the profile
    """

    # calculate the indices from the image
    y, x = np.indices(image.shape)

    # default to image center if no center given
    if not center:
        center = np.array([(x.max() - x.min()) / 2.0,
                           (x.max() - x.min()) / 2.0])

    r = np.hypot(x - center[0], y - center[1])

    # get sorted radii and sort image accordingly
    ind = np.argsort(r.flat)
    i_sorted = image.flat[ind]

    # for FP data we need to at least mask out data at
    # 0 or less so the gaps get ignored.
    # also want to mask out area outside of aperture
    # so use given maskval to do that.
    i_ma = ma.masked_less_equal(i_sorted, maskval)
    mask = ma.getmask(i_ma)

    # remove masked data points from further analysis
    r_sorted = ma.compressed(ma.array(r.flat[ind], mask=mask))
    i_mask = ma.compressed(i_ma)

    # get the integer part of the radii (bin size = 1)
    r_int = r_sorted.astype(int)

    # find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]  # assumes all radii represented
    rind = np.where(deltar)[0]       # location of changed radius
    nr_tot = rind[1:] - rind[:-1]    # total number of points in radius bin

    # cumulative sum to figure out sums for each radius bin
    csim = ma.cumsum(i_mask, dtype=float)
    tbin = csim[rind[1:]] - csim[rind[:-1]]

    # calculate and return profile of mean within each bin
    radial_prof = tbin / nr_tot

    return radial_prof
