def interpolation(template, flux_col = 1, wl_col =4, berv_flux, berv_wl):
    from scipy.interpolate import interp1d
    from scipy import interpolate
    from astropy.io import fits

    hdul = fits.open(template)


    # flux_template = hdul[flux_col].data
    wavelength_template = hdul[wl_col].data
    spline_coeff = interpolate.splrep(berv_wl[j], berv_flux)
    flux_interp_segment = interpolate.splev(wavelength_template[i], spline_coeff)

    return berv_wl, flux_interp_segment