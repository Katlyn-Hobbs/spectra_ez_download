import pytest
import ../ez_spectral_analysis as ezpz

def test_angstrom_units():

    # read in an example spectrum
    path_to_fits_files ='data/s2d/' 
    fits_files = os.listdir(path_to_fits_files)

    # Read .fits data
    # ---------------------------------------------------------------------------
    wavelengths, fluxes = ezpz.spectra_read(path=path_to_fits_files+fits_files[0])

    for wavelength in wavelengths:
        min_wavel, max_wavel = np.min(wavelengths), np.max(wavelengths)
        assert min_wavel >= 1000
        assert max_wavel <= 10000