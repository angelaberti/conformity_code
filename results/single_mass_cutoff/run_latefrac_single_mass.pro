PRO run_latefrac_single_mass
	zlist 		= [0.2, 0.4, 0.6, 0.8, 1.0]
	dz_coeffList 	= [1.0]

;	latefrac_single_mass_hist, 1.0, 0.2, 1.0

	FOR coeff_index=0,n_elements(dz_coeffList)-1 DO BEGIN
	  latefrac_single_mass_hist, dz_coeffList[coeff_index], 0.2, 1.0

	  FOR z_index=0,n_elements(zlist)-2 DO BEGIN
		latefrac_single_mass_hist, dz_coeffList[coeff_index], zlist(z_index), zlist(z_index+1)
	  ENDFOR
	ENDFOR
END
