class CUBISM_spectra:
    """Class to extract information from CUBISM spectra"""
    def __init__(self, path):
        self.path = path

def poly_wcs_ds9_region(path, wcs):
    from regions import read_ds9
    from astropy.wcs import WCS
    regions = read_ds9(path, errors='warn')
    regions = list(map(lambda x: x.to_pixel, regions))
    return regions 

class IRAC_img:
    from astropy.io import fits
    from astropy.wcs import WCS
    def __init__(self, path):
        self.path = path
        def hdu(path = path):
            hdu_IRAC = fits.open(path, mode='update')
            hdu_IRAC[0].header['CTYPE1']  = 'RA---TAN-SIP'
            hdu_IRAC[0].header['CTYPE2']  = 'DEC--TAN-SIP'
            hdu_IRAC.flush()
            return WCS(hdu_IRAC[0]), hdu_IRAC[0].data
        self.wcs, self.img = hdu(self.path)