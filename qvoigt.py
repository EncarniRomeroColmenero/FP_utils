import numpy as np
from scipy import special


def qvoigt(x, amp, pos, fwhm, gam):
    """
    analytic approximation for a voigt profile.

    V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))
    z = (x+i*gam)/(sig*sqrt(2))
    """

    tmp = 1 / special.wofz(
        np.zeros((len(x))) + 1j * np.sqrt(np.log(2.0)) * gam).real
    tmp = tmp * amp * special.wofz(
        2 * np.sqrt(np.log(2.0)) * (x - pos) / fwhm + 1j *
        np.sqrt(np.log(2.0)) * gam).real
    return tmp


def voigt_fwhm(gfwhm, gam):
    """
    approximation to calculate FWHM of a voigt profile given its
    gaussian FWHM, gfwhm, and lorentz parameter, gam.  good to at
    worst 0.0325% for pure lorentzian.
    """
    lfwhm = 2 * gam
    return 0.5346 * lfwhm + np.sqrt(0.2166 * lfwhm ** 2 + gfwhm ** 2)
