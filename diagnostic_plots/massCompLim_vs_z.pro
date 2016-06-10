PRO massCompLim_vs_z

	zbins = [0.20, 0.30, 0.40, 0.50, 0.65, 0.80, 1.00]

	zbinsMidpts = [0.2]
	FOR i=1,n_elements(zbins)-3 DO BEGIN
		zbinsMidpts = [zbinsMidpts, [mean([zbins[i],zbins[i+1]])]]
	ENDFOR
	zbinsMidpts = [zbinsMidpts,[1.0]]
;	PRINT, zbinsMidpts

	default_mean = [9.172,9.492,9.792,10.138,10.508,10.886]

	!p.multi=0
	erase

	PLOT, zbinsMidpts, default_mean, xrange=[0.2,1.0], yrange=[8,12]
	OPLOT, zbinsMidpts, default_mean+0.5
	OPLOT, zbinsMidpts, 0*indgen(7)+10
END
