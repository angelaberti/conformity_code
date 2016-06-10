PRO run_iptest
	dz_coeffList = [0.5, 1.5, 2., 2.5, 3.]

	FOREACH dz_coeff, dz_coeffList DO BEGIN
		iptest_default_params, dz_coeff
	ENDFOREACH
END
