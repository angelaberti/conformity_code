PRO run_latefrac_mass_bins
;	latefrac_mass_thirds_hist, 0.8, 1

	FOREACH dm,[0.8,1.0] DO latefrac_mass_thirds_hist, dm, 3
;	FOR t=2,3 DO FOREACH dm,[0.8,1.0] DO PRINT, 'latefrac_mass_thirds_hist, ', dm, ', ', t
END
