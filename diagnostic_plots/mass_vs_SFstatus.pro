PRO mass_vs_SFstatus, outputFormat
        IF (string(outputFormat) eq 'ps') THEN BEGIN
                SET_PLOT, 'ps'
;		DEVICE, file="mass_vs_SFstatus.ps", /LANDSCAPE
        ENDIF ELSE BEGIN
                SET_PLOT, 'X'
        ENDELSE

	binsize = 0.05
	linethickness = 2

	zmin = 0.2
	zmax = 1.0

	data = mrdfits("../results/zerodSFQ_all_cart.fits", 1)

	dataSF = data[where(data.SFQ eq 1)]
	dataQ  = data[where(data.SFQ eq 0)]
	
	plot, dataQ.zprimus, dataQ.mstar, psym=3, xrange=[zmin,zmax], yrange=[8,12.5], /nodata, xtitle="redshift", ytitle=textoidl("log (M/M_{\odot})")
	oplot, dataQ.zprimus, dataQ.mstar, psym=3, color=cgColor('Red')
	oplot, dataSF.zprimus, dataSF.mstar, psym=3, color=cgColor('Blue')

	dataSFzbinned = []
	dataQzbinned  = []
	zbins = []

	FOR z_index=0,(zmax-zmin)/binsize-1 DO BEGIN
		SFpoint = mean(dataSF[where( (dataSF.zprimus gt zmin+z_index*binsize) and $
					     (dataSF.zprimus le zmin+(z_index+1)*binsize) )].mstar)
		Qpoint  = mean( dataQ[where( (dataQ.zprimus gt zmin+z_index*binsize) and $
					     (dataQ.zprimus le zmin+(z_index+1)*binsize) )].mstar)
		zbin = mean([zmin+z_index*binsize, zmin+(z_index+1)*binsize])
		dataSFzbinned = [dataSFzbinned, SFpoint]
		dataQzbinned  = [dataQzbinned,  Qpoint]
		zbins = [zbins, zbin]
		print, zbin, SFpoint, Qpoint
	ENDFOR

	oplot, zbins, dataQzbinned, linestyle=0, thick=linethickness, color=cgColor('Black')
	oplot, zbins, dataSFzbinned, linestyle=0, thick=linethickness, color=cgColor('Black')

;	plothist, dataSF.zprimus,  bin=binsize, linestyle=0, thick=linethickness, yrange=[0,1600], xtitle="redshift"
;	plothist, dataQ.zprimus, bin=binsize, linestyle=2, thick=linethickness, /overplot

	LEGEND, ["Blue: Star-Forming", "Red:  Quiescent", "Solid Lines: mean of each population"], linestyle=!null, box=0

        IF (string(outputFormat) eq 'ps') THEN BEGIN
                DEVICE, /CLOSE
        ENDIF
       	SET_PLOT, 'X'
END
