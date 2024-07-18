import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import ez_spectral_analysis as ezpz

def test_spectra_upload():
    wavelength, flux = ezpz.spectra_read(path = 'spectra_ez_download/data/s2d/r.HARPN.2015-07-29T09-28-46.751_S2D_A.fits')
    assert len(wavelength)!=0
    assert len(flux)!=0

test_spectra_upload()

def test_spectra_plot():
    ezpz.plot_raw_data(path = 'spectra_ez_download/data/s2d/r.HARPN.2015-07-29T09-28-46.751_S2D_A.fits',orders=[16, 17, 18])

test_spectra_plot()