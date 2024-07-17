import numpy as np
from astropy.io import fits 
import pandas as pd
import matplotlib.pylab as plt

def spectra_read(path, flux_col = 1, wl_col = 4):

    hdul = fits.open(path)
    flux_data = hdul[flux_col].data
    wavelength_data = hdul[wl_col].data
    
    return wavelength_data, flux_data


def plot_raw_data(path, orders=[16, 17, 18]):
    wavelengths, fluxes =spectra_read(path=path)

    fig, axes = plt.subplots(len(orders), 1,figsize=(9, len(orders)*2))
    for index in range(0, len(orders)):

        ax = axes[index]
        ax.plot(wavelengths[orders[index]], fluxes[orders[index]], label='order = '+str(orders[index]),
                color='purple')
        ax.set_ylabel("Flux")
        ax.set_title("Order = "+str(orders[index]))
    
    ax.set_xlabel(r"Wavelength ($\AA$)")
    plt.tight_layout()
    plt.show()





  


