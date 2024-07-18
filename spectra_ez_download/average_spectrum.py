def average_spectrum(all_flux, wavelength):

    spectrum = fits.open(file_path)
            flux_data = spectrum[1].data['flux']
            wavelength_data = spectrum[4].data['wavelength']

            snr_values = []

            for i in range(40, 50):
                order_key = f'HIERARCH TNG QC ORDER{i} SNR'
                snr_value = spectrum[0].header[order_key]
                snr_values.append(snr_value)

            snr = np.mean(snr_values)

            for order_index in range(len(flux_data)):
                #mask1 = (wavelength_data[order_index] >= spectral_range[0]) & (wavelength_data[order_index] <= spectral_range[1])
                flux_per_order = flux_data[order_index]
                wavelength_per_order = wavelength_data[order_index]

                flux_data_all.append(flux_per_order)
                wavelength_data_all.append(wavelength_per_order)

            allsnr.append(snr)

        num_orders = len(flux_data_all)
        weights = [1 / (1 / snr_value) ** 2 for snr_value in allsnr]

        weighted_average_spectrum = np.average(flux_data_all, axis=0, weights=weights) #I want to taske the avwerage for each order tho so maybe this would be up?
        wavelengths = wavelength_data_all[0]  #