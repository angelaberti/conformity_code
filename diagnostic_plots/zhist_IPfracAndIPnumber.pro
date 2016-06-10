PRO zhist_IPfracAndIPnumber, outputFormat
        IF (string(outputFormat) eq 'ps') THEN BEGIN
                SET_PLOT, 'ps'
		DEVICE, file="../figures/z_vs_IPfracAndIPnumber.ps", /LANDSCAPE
        ENDIF ELSE BEGIN
                SET_PLOT, 'X'
        ENDELSE

;	!p.multi=[0,2,2]

;	!x.margin=[7,2.5] ; left, right
;	!y.margin=[3.5,1] ; bottom, top

	binsize = 0.02	
	linethickness = 1

	data = mrdfits("../results/default_parameters/IP_data/zerodSFQ_IP_dz1.0_dm0.5.fits", 1)
	
	dataIP  = data[where(data.IP eq 1)]
	dataNIP = data[where(data.IP eq 0)]

	plothist, dataIP.zprimus,  bin=binsize, linestyle=0, thick=linethickness, yrange=[0,1700], xtitle="redshift"
	plothist, dataNIP.zprimus, bin=binsize, linestyle=2, thick=linethickness, /overplot

	LEGEND, ["IP", "non-IP"], linestyle=[0,2], box=0, charsize=1

	zmin = 0.2
	zmax = 1.0
	dz = binsize

	zArray = [zmin]
	FOR z_index=0,(zmax-zmin)/dz DO BEGIN
		zArray = [zArray, zmin+float(z_index)*dz]
	ENDFOR
;	print, zArray		

	IPfrac_zbinned = []
	n_IParray = []
	FOR z_index=1,(zmax-zmin)/dz DO BEGIN
		data_zbin = data[where( (data.zprimus ge zArray[z_index]) AND $
					(data.zprimus lt zArray[z_index+1]) )]
		data_zbinIP  = data_zbin[where(data_zbin.IP eq 1)]
		data_zbinNIP = data_zbin[where(data_zbin.IP eq 0)]

		n_IP 	= n_elements(data_zbinIP)
		n_total = n_elements(data_zbin)
		
		n_IParray = [n_IParray, float(n_IP)]
		IPfrac_zbinned = [IPfrac_zbinned, float(n_IP)/float(n_total)]

		PRINT, zArray[z_index], zArray[z_index+1], n_IP, n_total, float(n_IP)/float(n_total)
	ENDFOR

	plot, zArray+dz/2., n_IParray, xtitle="redshift", ytitle=textoidl("N_{IP}"), xrange=[zmin,zmax], yrange=[0,1500]
	plot, zArray+dz/2, IPfrac_zbinned, xtitle="redshift", ytitle="IP fraction", xrange=[zmin,zmax], yrange=[0,1]

        IF (string(outputFormat) eq 'ps') THEN BEGIN
                DEVICE, /CLOSE
        ENDIF
       	SET_PLOT, 'X'
END
