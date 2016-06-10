PRO massIPhist, outputFormat
        IF (string(outputFormat) eq 'ps') THEN BEGIN
                SET_PLOT, 'ps'
		DEVICE, file="../figures/massIPhist.ps", /LANDSCAPE
        ENDIF ELSE BEGIN
                SET_PLOT, 'X'
        ENDELSE

	binsize = 0.1	
	linethickness = 1
	ymin = 0
	ymax = 3000
	
	data = mrdfits("../results/default_parameters/IP_data/zerodSFQ_IP_dz1.0_dm0.5.fits", 1)
	
	dataSF = data[where(data.SFQ eq 1)]
	dataQ  = data[where(data.SFQ eq 0)]

	plothist, dataSF.mstar,  bin=binsize, linestyle=0, thick=linethickness, yrange=[ymin,ymax], xtitle="log (M/Msolar)"
	plothist, dataQ.mstar, bin=binsize, linestyle=2, thick=linethickness, /overplot

	LEGEND, ["Star-Forming", "Quiescent"], linestyle=[0,2], box=0

        IF (string(outputFormat) eq 'ps') THEN BEGIN
                DEVICE, /CLOSE
        ENDIF
       	SET_PLOT, 'X'
END
