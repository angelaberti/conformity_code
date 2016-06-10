PRO z_vs_IPfrac_SFstatus, outputFormat
        IF (string(outputFormat) eq 'ps') THEN BEGIN
		SET_PLOT, 'ps'
		DEVICE, file="../figures/z_vs_IPfrac_SFstatus.ps", /LANDSCAPE
        ENDIF ELSE BEGIN
                SET_PLOT, 'X'
        ENDELSE

	!p.multi=0

;	!x.margin=[7,2.5] ; left, right
;	!y.margin=[3.5,1] ; bottom, top

	binsize = 0.1
	linethickness = 1

	data = mrdfits("../results/default_parameters/IP_data/zerodSFQ_IP_dv1500kms_dm0.5.fits", 1)
	
	dataIP  = data[where(data.IP eq 1)]
	dataNIP = data[where(data.IP eq 0)]

;	plothist, dataIP.zprimus,  bin=binsize, linestyle=0, thick=linethickness, yrange=[0,1700], xtitle="redshift"
;	plothist, dataNIP.zprimus, bin=binsize, linestyle=2, thick=linethickness, /overplot

;	LEGEND, ["IP", "non-IP"], linestyle=[0,2], box=0, charsize=1

	zmin = 0.2
	zmax = 1.0
	dz = binsize

	zArray = [zmin]
	FOR z_index=0,(zmax-zmin)/dz DO BEGIN
		zArray = [zArray, zmin+float(z_index)*dz]
	ENDFOR

	data_zbinIPSFarray 	= []
	data_zbinIPQarray 	= []
	data_zbinSFarray 	= []
	data_zbinQarray 	= []

	FOR z_index=1,(zmax-zmin)/dz DO BEGIN

		data_zbin 	= data[where( (data.zprimus ge zArray[z_index]) AND $
					      (data.zprimus lt zArray[z_index+1]) )]
		data_zbinIP  	= data_zbin[where(data_zbin.IP eq 1)]
;		data_zbinNIP 	= data_zbin[where(data_zbin.IP eq 0)]

		data_zbinIPSF 	= n_elements(data_zbinIP[where(data_zbinIP.SFQ eq 1)])
		data_zbinIPQ  	= n_elements(data_zbinIP[where(data_zbinIP.SFQ eq 0)])
		data_zbinSF 	= n_elements(data_zbin[where(data_zbin.SFQ eq 1)])
		data_zbinQ  	= n_elements(data_zbin[where(data_zbin.SFQ eq 0)])

		data_zbinIPSFarray 	= [data_zbinIPSFarray, 	float(data_zbinIPSF)]
		data_zbinIPQarray 	= [data_zbinIPQarray, 	float(data_zbinIPQ)]
		data_zbinSFarray 	= [data_zbinSFarray, 	float(data_zbinSF)]
		data_zbinQarray 	= [data_zbinQarray, 	float(data_zbinQ)]

		PRINT, zArray[z_index], zArray[z_index+1], float(data_zbinSF), float(data_zbinQ), $
							   float(data_zbinIPSF)/float(data_zbinSF), float(data_zbinIPQ)/float(data_zbinQ)
	ENDFOR
	
;	print, data_zbinIPSFarray, data_zbinQarray

	plot, zArray+dz/2., float(data_zbinIPSFarray/data_zbinSFarray), color=cgColor('White'), xtitle="redshift", ytitle="IP fraction by SF status", xrange=[zmin,zmax], yrange=[0,1], /nodata
	oplot, zArray+dz/2., float(data_zbinIPSFarray/data_zbinSFarray), color=cgColor('Blue')
	oplot, zArray+dz/2., data_zbinIPQarray/data_zbinQarray, color=cgColor('Red')
	LEGEND, ["Late-Type", "Early-Type"], linestyle=[0,0], color=[cgColor('Blue'),cgColor('Red')], box=0

        IF (string(outputFormat) eq 'ps') THEN BEGIN
                DEVICE, /CLOSE
        ENDIF
       	SET_PLOT, 'X'
END
