PRO plot_completeness
	xmin=0.2
	xmax=1.0

	dR = 0.1 ; Mpc
	dz_coeff = 2.0
	targPhi = -3.7
	stringPHI = strtrim(string(-1.*targPhi, format='(f20.1)'), 2)

	dataAll = MRDFITS('~/results/default_parameters/allAboveMassCompLim-0.5.fits',1)
	dataAll = dataAll[WHERE(dataAll.targ_weight GE 1.)]

	dataIP  = MRDFITS('~/results/default_parameters/IP_data/zerodSFQ_IP_dz2.0_dm0.5.fits',1)
	dataIP  = dataIP[WHERE((dataIP.IP EQ 1) AND (dataIP.targ_weight GE 1.))]

	dataAllcons = MRDFITS('~/results/conservative_mass_cutoff/allAboveMassCompLim-0.0.fits',1)
	dataAllcons = dataAllcons[WHERE(dataAllcons.targ_weight GE 1.)]

	dataIPcons  = MRDFITS('~/results/match_IP_sample_rigorous/correlation_test/allAboveMassCompLim+0.5_IP.fits',1)
	dataIPcons  = dataIPcons[WHERE((dataIPcons.IP EQ 1) AND (dataIPcons.targ_weight GE 1.))]

	dataIP_M13  = MRDFITS('~/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits',1)
	dataIP_M13  = dataIP_M13[WHERE((dataIP_M13.IP EQ 1) AND (dataIP_M13.targ_weight GE 1.))]

	fields = dataAll[uniq(dataAll.field, sort(dataAll.field))].field

	ERASE
	!P.MULTI=[6,3,2]
;	!P.MULTI=1

;	PLOT, dataAll.zprimus, dataAll.mstar, xrange=[0.2,1.0], yrange=[8.5,12], /NODATA
;	OPLOT, dataAll.zprimus, dataAll.mstar, PSYM=3, color=cgColor('red')
;	OPLOT, dataIP_M13.zprimus, dataIP_M13.mstar, PSYM=3
;	OPLOT, [xmin,xmax],[10,10], LINESTYLE=0, color=cgColor('green'), THICK=3
;	OPLOT, [xmin,xmax],[10.5,10.5], LINESTYLE=0, color=cgColor('green'), THICK=3
;	OPLOT, [xmin,xmax],[11,11], LINESTYLE=0, color=cgColor('green'), THICK=3

;	PLOT, dataAllcons.zprimus, dataAllcons.mstar, xrange=[0.2,1.0], yrange=[8.5,12], /NODATA
;	OPLOT, dataAllcons.zprimus, dataAllcons.mstar, PSYM=3, color=cgColor('red')
;	OPLOT, dataIPcons.zprimus, dataIPcons.mstar, PSYM=3
;	OPLOT, [xmin,xmax],[10.1,10.1], LINESTYLE=0, color=cgColor('green'), THICK=3
;	OPLOT, [xmin,xmax],[10.5,10.5], LINESTYLE=0, color=cgColor('green'), THICK=3
;	OPLOT, [xmin,xmax],[11,11], LINESTYLE=0, color=cgColor('green'), THICK=3

	FOR f=0,n_elements(fields)-1 DO BEGIN	
	dataAllconsField = dataAllcons[WHERE(dataAllcons.field EQ fields[f])]
	dataIPconsField = dataIPcons[WHERE(dataIPcons.field EQ fields[f])]
	
	PLOT, dataAllconsField.zprimus, dataAllconsField.mstar, xrange=[0.2,1.0], yrange=[8.5,12], /NODATA
	XYOUTS, 0.6, 8.6, STRTRIM(fields[f],2), CHARSIZE=2, ALIGNMENT=0.5
	OPLOT, dataAllconsField.zprimus, dataAllconsField.mstar, PSYM=3, color=cgColor('red')
	OPLOT, dataIPconsField.zprimus, dataIPconsField.mstar, PSYM=3
	OPLOT, [xmin,xmax],[10.1,10.1], LINESTYLE=0, color=cgColor('green'), THICK=3
	OPLOT, [xmin,xmax],[10.4,10.4], LINESTYLE=0, color=cgColor('green'), THICK=3
	OPLOT, [xmin,xmax],[10.7,10.7], LINESTYLE=0, color=cgColor('green'), THICK=3
	OPLOT, [xmin,xmax],[11,11], LINESTYLE=0, color=cgColor('green'), THICK=3

	ENDFOR

;	PLOT, dataAllcons.zprimus, dataAllcons.mstar, xrange=[0.2,1.0], yrange=[8.5,12], /NODATA
;	OPLOT, dataAllcons.zprimus, dataAllcons.mstar, PSYM=3, color=cgColor('red')
;	OPLOT, dataIPcons.zprimus, dataIPcons.mstar, PSYM=3
;	OPLOT, [xmin,xmax],[10.1,10.1], LINESTYLE=0, color=cgColor('green'), THICK=3
;	OPLOT, [xmin,xmax],[10.4,10.4], LINESTYLE=0, color=cgColor('green'), THICK=3
;	OPLOT, [xmin,xmax],[10.7,10.7], LINESTYLE=0, color=cgColor('green'), THICK=3
;	OPLOT, [xmin,xmax],[11.,11.], LINESTYLE=0, color=cgColor('green'), THICK=3

END
