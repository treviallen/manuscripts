1) run extract_iris_stress_drop_data.py
2) run get_cwb_seed.py & extract_legacy_cwb_data.py
3) run scan_record_snr.py to get pick files
3) run get_sn_fft_spectra.py - makes "fft_data.pkl"
4) calc_regional_geometric_spreading.py - makes simple attenuation model (currently just at 0.75 Hz)
5) run calc_site_kappa.py - estimates site kappa-zero
6) run calc_site_kappa_dist.py - estimates distance dependent kappa
7) run fit_brune_spectra.py

Testing:
review_event_residuals.py