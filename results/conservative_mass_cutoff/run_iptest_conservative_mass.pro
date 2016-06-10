PRO run_iptest_conservative_mass
;	dz_coeffList = [0.5, 1.5, 2., 2.5, 3.]
	dz_coeffList = [1.0]	

	FOREACH dz_coeff, dz_coeffList DO BEGIN
		iptest_conservative_mass, dz_coeff
	ENDFOREACH
END
