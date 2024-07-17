# analyze spectra!

tasks:
1. write function that does the continuum normalization correction - zoe
   
   a. input (filepath, wavelength_for_all_orders, flux_for_all_orders)
   
   b. output (wavelength_for_all_orders, normalized_flux_for_all_orders)

3. write function that performs barycentric correction for all spectral orders - nina

   a. input(wavelength, flux, date/time, BERV from fits header)

   b. output(wavelength, flux)

5. Wavelength interpolation function - kate

   a. input(template, flux)

6. write function that computes weighted daily average - kate

7. write a simple plotting function - all
