import numpy as N
from scipy import special

def voigt(x,amp,pos,fwhm,shape):
     """\
     voigt profile

     V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))
     z = (x+i*gam)/(sig*sqrt(2))
     """

     tmp = 1/special.wofz(N.zeros((len(x))) \
           +1j*N.sqrt(N.log(2.0))*shape).real
     tmp = tmp*amp* \
           special.wofz(2*N.sqrt(N.log(2.0))*(x-pos)/fwhm+1j* \
           N.sqrt(N.log(2.0))*shape).real
     return tmp
