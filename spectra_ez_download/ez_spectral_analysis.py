import numpy as np
from astropy.io import fits 
import pandas as pd
import matplotlib.pylab as plt

# new packages. needs to be added to requirements and readthedocs.yaml
from specutils.spectra import Spectrum1D, SpectralRegion
from specutils.fitting import fit_generic_continuum
from astropy import units as u
import astropy.constants as ac
from astropy.modeling import models
from scipy.interpolate import interp1d
from scipy import interpolate

class ez_spectral_analysis(object):
    """
    Completes simple spectral analysis for user-uploaded solar spectra of HARPS-N.

    Args:
        path (str): path to local data
        instrument (str): specifies the instrument used. Currently, the only supported instrumented is 'HARPS-N'.   
        orders (list, optional): list of orders to be inspected. Default is [16, 17, 18]. Beware of noisy orders and tellurics.
        plot (bool, optional): produce plots from analysis. 
        median_window (int, optional): integer. Needs to be an odd integer. The width of the median smoothing kernel used to filter the data before fitting the continuum. 
        
    """

    def __init__(self, path, instrument, orders = [16, 17, 18], plot = False, median_window = 3):
        self.path = path
        self.instrument = instrument
        self.orders = orders
        self.plot = plot
        self.median_window = median_window
        
        
    def spectra_read(self):
        """
        Read spectra
        Open fits files and return the wavelength and spectra for each wavelength order.
    
        Returns:
            arrays: wavelength array, flux array
        """
    
        if self.instrument =='HARPS-N':
          self.flux_col = 1
          self.wl_col = 4
          hdul = fits.open(self.path)
          flux_data = hdul[self.flux_col].data
          wavelength_data = hdul[self.wl_col].data
          return wavelength_data, flux_data
        else:
           print("Error: "+self.instrument+" is currently not supported by spectra_ez_download. Please submit an issue on github.")

    def plot_raw_data(self):
        """
        plot the raw spectra

        Create a figure for the raw spectra for the specified wavelength orders. If orders are not specified, 
        the default behavior is to plot orders 16, 17, 18.
        """
        wavelengths, fluxes = ez_spectral_analysis.spectra_read(self)
    
        fig, axes = plt.subplots(len(self.orders), 1,figsize=(9, len(self.orders)*2))
        for index in range(0, len(self.orders)):
    
            ax = axes[index]
            ax.plot(wavelengths[self.orders[index]], fluxes[self.orders[index]], label='order = '+str(self.orders[index]),
                    color='purple')
            ax.set_ylabel("Flux")
            ax.set_title("Order = "+str(self.orders[index]))
        
        ax.set_xlabel(r"Wavelength ($\AA$)")
        plt.tight_layout()
        plt.show()
    
    
    def continuum_norm(self):
      """
      normalize the continuum of spectra
    
      Opens fits files, cycles through each wavelength order (or the user-specified orders) and normalizes the continuum using specutils
    
      Returns:
            arrays: wavelength array, normalized flux array
      """
    
      wavelengths, fluxes = ez_spectral_analysis.spectra_read(self)
    
      norm_fluxes = []
      
      if self.orders == None:
        for order in range(0, len(fluxes)):
          spectrum = Spectrum1D(flux=fluxes[order].astype(np.float64)*u.Jy, 
                                spectral_axis=wavelengths[order].astype(np.float64)/10*u.nm)
          continuum_g_x = fit_generic_continuum(spectrum, median_window=self.median_window, 
                                                model=models.Chebyshev1D(4, c0=0., c1=0., c2=0., c3=0.))
          y_continuum_fitted = continuum_g_x(wavelengths[order]/10*u.nm)
    
          spec_normalized = spectrum / y_continuum_fitted
          norm_fluxes.append(spec_normalized.flux)
    
          if self.plot:
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
            f.tight_layout() 
      else:
        for order in self.orders:
          spectrum = Spectrum1D(flux=fluxes[order].astype(np.float64)*u.Jy, 
                                spectral_axis=wavelengths[order].astype(np.float64)/10*u.nm)
          continuum_g_x = fit_generic_continuum(spectrum, median_window=self.median_window, 
                                                model=models.Chebyshev1D(4, c0=0., c1=0., c2=0., c3=0.))
          y_continuum_fitted = continuum_g_x(wavelengths[order]/10*u.nm)
    
          spec_normalized = spectrum / y_continuum_fitted
          norm_fluxes.append(spec_normalized.flux)
    
          if self.plot:
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
            f.tight_layout()    
    
      return wavelengths.astype(np.float64), norm_fluxes        
    
    def BC_correction(self):
        """
        Apply a barycentric correction
        
          Apply a barycentric correction to wavlelength values based on Barycentric velocity and BJD values read from .fits header. Reads in flux and wavelength info from continuum_norm().
    
          Returns:
              arrays: corrected wavelengths, flux
        """
        # load in wavelength and flux:
        wavelength, flux = ez_spectral_analysis.continuum_norm(self)
        # load in BERV and BJD from header:
        hdul = fits.open(self.path)
        hdr = hdul[0].header
        berv = hdr['*BERV'][0] # km/s
        berv_ms = berv*1000
        bjd = hdr['*BJD'][0] # Barycentric Julian date (TDB)

        rv_corr = 1 + berv_ms/ac.c.value # corrective value

        # iterate through all the orders and apply the same barycentric correction to the wavelengths    
        wavelength_corrected_list = []
        for order in range(len(wavelength)):
            wavelength_corrected = wavelength[order]*rv_corr # angstroms
            wavelength_corrected_list.append(wavelength_corrected)
    
        return wavelength_corrected_list, flux

    def interpolation(self, template):   
        """
        Apply a spline interpolation.
            Perform and apply a spline interpolation to shift each spectrum to a common wavelength grid, determined by the input template wavelength grid from a .fits file. 
    
          Args:
              template (.fits file): Path to .fits file that will be used as the template wavelength grid for interpolation. 
    
          Returns:
              arrays: wavelengths, interpolated flux
        """
        berv_wl, berv_flux = ez_spectral_analysis.BC_correction(self)
        hdul = fits.open(template)

        wavelength_template = hdul[self.wl_col].data
        flux_interpolated = []
        for order in range(len(berv_wl)):
            spline_coeff = interpolate.splrep(berv_wl[order], berv_flux[order])
            flux_interpolated.append(interpolate.splev(wavelength_template[order], spline_coeff))
    
        return wavelength_template, flux_interpolated
             
    
          
    
    
    




  


