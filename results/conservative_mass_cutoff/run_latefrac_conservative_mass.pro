PRO run_latefrac_conservative_mass
	zlist 		= [0.2, 0.4, 0.6, 0.8, 1.0]
;	dz_coeffList 	= [0.5, 1., 1.5, 2, 2.5, 3.]
;	dz_coeffList	= [1.0]

	latefrac_conservative_mass_hist, 0.2, 1.0

	FOR z_index=0,n_elements(zlist)-2 DO BEGIN
		latefrac_conservative_mass_hist, zlist(z_index), zlist(z_index+1)
	ENDFOR
END
