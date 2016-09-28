import numpy as np

from astropy.io import fits


def create_fits_with_header(header, f):
    """
    Create a FITS file with header information and an "empty" image.

    This functions creates a FITS file whose image consists of an array of zeros. The header details of the FITS file
    are taken from the given header dictionary. String representations are used for all the values.

    The FITS file is written to the file-like object f. An example usage would be:

    with open('/path/to/file.fits', 'w') as f:
        create_fits_with_header(dict(B=38), f)

    Params:
    -------
    header: dict
       Dictionary of header values for the FITS file.
    f: file-like object
       File-like object to which the FITS file content is output.
    """

    image = np.zeros(100)
    hdu_header = fits.Header(cards=header)

    hdu = fits.PrimaryHDU(image, header=hdu_header)

    hdu.writeto(f)


