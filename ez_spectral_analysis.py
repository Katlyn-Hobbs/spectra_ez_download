import numpy as np
from astropy.io import fits 
import pandas as pd

def spectra_read(path):
    hdul = fits.open(path)
    flux_data = hdul[1].data
    wavelength_data = hdul[4].data

  


