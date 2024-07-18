import os
import sys
import inspect
import numpy as np

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import ez_spectral_analysis as ezpz
import pytest

def test_angstrom_units():

    # read in an example spectrum
    path_to_fits_files ='spectra_ez_download/data/s2d/' 
    fits_files = os.listdir(path_to_fits_files)

    # Read .fits data
    # ---------------------------------------------------------------------------
    wavelengths, fluxes = ezpz.spectra_read(path=path_to_fits_files+fits_files[0])

    for wavelength in wavelengths:
        min_wavel, max_wavel = np.min(wavelengths), np.max(wavelengths)
        assert min_wavel >= 1000, "Wavelengths need to be in angstrom. Please convert your wavelengths."
        assert max_wavel <= 10000, "Wavelengths need to be in angstrom. Please convert your wavelengths."


