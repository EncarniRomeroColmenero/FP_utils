import math

from astropy.io import fits


def A(ring_file, R, W_line):
    """Calculate the A value for ring.

    A is calculated for the ring parameters in the given FITS file and the given ring radius R as well as emission line
    wavelength W_line.

    Params:
    -------
    ring_file: file path or file object or file-like object
        A file-like object or a file path for the FITS file from which to get ring parameters.
    R: float
        Ring radius in unbinned pixels.
    W_line: float
        Lab wavelength of the the emission line used (in Angstroms).

    Returns:
    --------
    float:
        The value of A.
    """

    # read headers from file
    with fits.open(ring_file) as hdulist:
        prihdr = hdulist[0].header

    missing_card = lambda key: 'Card missing in FITS header: {key}'.format(key=key)

    # make sure an etalon state is defined
    ETALON_STATE = 'ET-STATE'
    if ETALON_STATE not in prihdr:
        raise ValueError(missing_card(ETALON_STATE))

    # parse the etalon state
    etalon_state = prihdr[ETALON_STATE]
    if etalon_state == 'S2 - Etalon 1':
        etalon = 1
    elif etalon_state == 'S3 - Etalon 2':
        etalon = 2
    elif etalon_state == 'S4 - Etalon 1 & 2':
        etalon = 1
    else:
        raise ValueError('Unsupported value for {key} card: {value}'.format(key=ETALON_STATE, value=etalon_state))

    # make sure all ring parameters are defined
    ETALON_B = 'ET{etalon}B'.format(etalon=etalon)
    ETALON_F = 'ET{etalon}F'.format(etalon=etalon)
    ETALON_Z = 'ET{etalon}Z'.format(etalon=etalon)
    for key in [ETALON_B, ETALON_F, ETALON_Z]:
        if key not in prihdr:
            raise ValueError(missing_card(key))

    # collect the ring parameters
    B = prihdr[ETALON_B]
    F = prihdr[ETALON_F]
    Z = prihdr[ETALON_Z]

    # calculate A
    W_central = W_line * math.sqrt(1 + (float(R) / F) ** 2)
    return W_central - B * Z