#VLA primary beam response
import numpy
from astropy.io import fits

def getpixsize(hdr, cunit, cdelt):
    """Get the pixel size from a header and convert
       it to arcsec if necessary"""
    if hdr[cunit].lower()[:3] == 'deg':
        pixsize = abs(hdr[cdelt])*3600
    elif hdr[cunit].lower()[:6] == 'arcsec':
        pixsize = abs(hdr[cdelt])
    elif hdr[cunit].lower()[:6] == 'arcmin':
        pixsize = abs(hdr[cdelt])*60
    else:
        mes = ("Make sure that CUNIT1/CUNIT2 is"
              " either in degree, arcmin or arcsec")
        raise UnrecognizedUnits(mes)
    return pixsize

class Getpbresponse:
    """Calculates beam response of the VLA as 
    a function of radial distance.

    Parameters
    ----------
    inpmap: string
        name of 2D fits file from which coordinates are taken
    rad: array
        input radius in arcminute, i.e, the maximum radius  
    
    Methods:
    --------
    get_vla_pbc()
        return 2D fits file containing the primary beam response values. 

    Examples
    --------
    response = Getpbresponse("fits_file.fits", rad=30, output="vla_pbcmap.fits") #30 arcminute
    vla_beam = response.get_vla_pbc()
    """
    def __init__(self, inpmap, rad, output="vla_pbcmap.fits"):
        self.img, self.hdr = fits.getdata(inpmap, header=True)
        self.img = numpy.squeeze(self.img)
        self.nx = self.img.shape[-1]
        self.ny = self.img.shape[-2]
        Y, X = numpy.indices((self.ny, self.nx))
        self.X0 = int(self.hdr['CRPIX1']) - 1
        self.Y0 = int(self.hdr['CRPIX2']) - 1
        self.pix1 = abs(getpixsize(self.hdr, "CUNIT1", "CDELT1"))
        self.pix2 = abs(getpixsize(self.hdr, "CUNIT2", "CDELT2"))
        self.radmaxarcmin = rad 
        self.x = X - self.X0
        self.y = Y - self.Y0
        self.output=output

    def get_vla_pbc(self):
        rad = numpy.hypot(self.x, self.y) # radius in pixels
        radarcmin = rad * self.pix1 / 60 # radius in arcmin
        pbcfac = self.vla_pbc(radarcmin)
        pbcfac[radarcmin > 29.8] = self.vla_pbc(29.8)
        pbcfac[radarcmin > self.radmaxarcmin] = numpy.nan
        fits.writeto(self.output, pbcfac,
                        self.hdr, overwrite=True)
        
    def vla_pbc(self, r):
        G1 = -1.343E-3
        G2 = 6.579E-7
        G3 = -1.186E-10
        f = 1.4 # 1.4 GHz
        dist = r * f
        return (1 + G1 * dist**2 + G2 * dist**4
                + G3 * dist**6)