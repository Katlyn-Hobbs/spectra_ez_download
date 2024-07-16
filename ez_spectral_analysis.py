import numpy as np
from astropy.io import fits 
import pandas as pd

def spectra_read(path, flux_col = 1, wl_col = 4):

    hdul = fits.open(path)
    flux_data = hdul[flux_col].data
    wavelength_data = hdul[wl_col].data
    
    return wavelength_data, flux_data
  


