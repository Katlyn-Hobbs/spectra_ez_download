import numpy as np
from astropy.io import fits 
import pandas as pd
import matplotlib.pylab as plt

# new packages. needs to be added to requirements and readthedocs.yaml
from specutils.spectra import Spectrum1D, SpectralRegion
from specutils.fitting import fit_generic_continuum
from astropy import units as u
from astropy.modeling import models

def spectra_read(path, instrument):
    """
    Read spectra

    Open fits files and return the wavelength and spectra for each wavelength order.

    Args:
        path (str): string. Specified path of where a spectrum is located
        instrument (str): String that specifies the instrument. Currently, the only supported instrumented is 'HARPS-N'.

    Returns:
        arrays: wavelength array, flux array

    """

    if instrument =='HARPS-N':
      flux_col = 1, wl_col = 4
      hdul = fits.open(path)
      flux_data = hdul[flux_col].data
      wavelength_data = hdul[wl_col].data
      return wavelength_data, flux_data
    else:
       print("Error: "+instrument+" is currently not supported by spectra_ez_download. Please submit an issue on github.")




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


def continuum_norm(path, plot=False, median_window=3, orders=None):
  """
  normalize the continuum of spectra

    Opens fits files, cycles through each wavelength order (or the user-specified orders) and normalizes the continuum using specutils

    Args:
        path (str): string. Specified path of where a spectrum is located
        plot (bool): boolean. If True, plots will be generated. Default is False
        orders (array): array. array of wavelength orders to be normalized. If not specified, all wavelength orders will be used.
        median_window (int): integer. Needs to be an odd integer. The width of the median smoothing kernel used to filter the data before fitting the continuum.

    Returns:
        arrays: wavelength array, normalized flux array
  """

  wavelengths, fluxes =spectra_read(path=path)

  norm_fluxes = []
  
  if orders == None:
    for order in range(0, len(fluxes)):
      spectrum = Spectrum1D(flux=fluxes[order].astype(np.float64)*u.Jy, 
                            spectral_axis=wavelengths[order].astype(np.float64)/10*u.nm)
      continuum_g_x = fit_generic_continuum(spectrum, median_window=median_window, 
                                            model=models.Chebyshev1D(4, c0=0., c1=0., c2=0., c3=0.))
      y_continuum_fitted = continuum_g_x(wavelengths[order]/10*u.nm)

      spec_normalized = spectrum / y_continuum_fitted
      norm_fluxes.append(spec_normalized.flux)

      if plot:
        f, axes = plt.subplots(2, 1, figsize=(12, 4))
        ax = axes[0]  
        ax.plot(wavelengths[order]/10*u.nm, fluxes[order]*u.Jy, color='purple')  
        ax.plot(wavelengths[order]/10*u.nm, y_continuum_fitted)  
        ax.set_ylabel('Flux')  
        ax.set_title("Continuum Fitting")

        ax = axes[1]
        ax.plot(spec_normalized.spectral_axis, spec_normalized.flux, color='purple')    
        ax.set_title("Continuum normalized spectrum")
        ax.set_ylabel('Normalized flux')  
        ax.set_xlabel('Wavelengths (nm)')  
        f.tight_layout 
  else:
    for order in orders:
      spectrum = Spectrum1D(flux=fluxes[order].astype(np.float64)*u.Jy, 
                            spectral_axis=wavelengths[order].astype(np.float64)/10*u.nm)
      continuum_g_x = fit_generic_continuum(spectrum, median_window=median_window, 
                                            model=models.Chebyshev1D(4, c0=0., c1=0., c2=0., c3=0.))
      y_continuum_fitted = continuum_g_x(wavelengths[order]/10*u.nm)

      spec_normalized = spectrum / y_continuum_fitted
      norm_fluxes.append(spec_normalized.flux)

      if plot:
        f, axes = plt.subplots(2, 1, figsize=(12, 4))
        ax = axes[0]  
        ax.plot(wavelengths[order]/10*u.nm, fluxes[order]*u.Jy, color='purple')  
        ax.plot(wavelengths[order]/10*u.nm, y_continuum_fitted)
        ax.set_ylabel('Flux')   
        ax.set_title("Continuum Fitting")

        ax = axes[1]
        ax.plot(spec_normalized.spectral_axis, spec_normalized.flux, color='purple')
        ax.set_xlabel('Wavelengths (nm)')
        ax.set_ylabel('Normalized flux')        
        ax.set_title("Continuum normalized spectrum")
        f.tight_layout    

  return wavelengths.astype(np.float64), norm_fluxes        
         

      







  


