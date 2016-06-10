PRO run_latefrac_default_params
	zlist 		= [0.2, 0.4, 0.6, 0.8, 1.0]
	dz_coeffList 	= [0.5, 1., 1.5, 2, 2.5, 3.]

;	latefrac_default_params_hist, 1.0, 0.2, 1.0
;	latefrac_default_params, 1.0, 0.2, 0.4

	FOR coeff_index=0,n_elements(dz_coeffList)-1 DO BEGIN

	  FOR z_index=0,n_elements(zlist)-2 DO BEGIN
		latefrac_default_params_hist, dz_coeffList[coeff_index], zlist(z_index), zlist(z_index+1)
	  ENDFOR

	  latefrac_default_params_hist, dz_coeffList[coeff_index], 0.2, 1.0

	ENDFOR
END
