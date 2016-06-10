PRO zhist_SFstatus, outputFormat
        IF (string(outputFormat) eq 'ps') THEN BEGIN
		SET_PLOT, 'ps'
		DEVICE, file="../figures/zhist_SFstatus.ps", /LANDSCAPE
        ENDIF ELSE BEGIN
                SET_PLOT, 'X'
        ENDELSE

	binsize = 0.02	
	linethickness = 1

	data = mrdfits("../results/zerodSFQ_all_cart.fits", 1)
	
	dataSF = data[where(data.SFQ eq 1)]
	dataQ  = data[where(data.SFQ eq 0)]

	plothist, data.zprimus,  bin=binsize, linestyle=0, thick=linethickness, yrange=[0,2500], xtitle="redshift"
	plothist, dataSF.zprimus, bin=binsize, linestyle=1, thick=linethickness, /overplot
	plothist, dataQ.zprimus, bin=binsize, linestyle=2, thick=linethickness, /overplot

	LEGEND, ["All", "Star-Forming", "Quiescent"], linestyle=[0,1,2], box=0

        IF (string(outputFormat) eq 'ps') THEN BEGIN
                DEVICE, /CLOSE
        ENDIF
  	SET_PLOT, 'X'
END
