PRO zhist_SFandIPstatus, outputFormat
        IF (string(outputFormat) eq 'ps') THEN BEGIN
		SET_PLOT, 'ps'
		DEVICE, file="../figures/zhist_SFandIPstatus.ps", /LANDSCAPE
        ENDIF ELSE BEGIN
                SET_PLOT, 'X'
        ENDELSE

	binsize = 0.02
	linethickness = 2

	A = FINDGEN(17) * (!PI*2/16.)
	USERSYM, 0.25*COS(A), 0.25*SIN(A), /FILL
	ptsymbol = 8

	zmin = 0.2
	zmax = 1.0

	ymin = 0
	ymax = 1500
	
;	path = '../results/conservative_mass_cutoff/IP_data/'
	path = '../results/default_parameters/IP_data/'
	data = mrdfits(path+'zerodSFQ_IP_dv1500kms_dm0.5.fits', 1)

	dataSF = data[where(data.SFQ eq 1)]
	dataQ  = data[where(data.SFQ eq 0)]

	dataSF_IP  = dataSF[where(dataSF.IP eq 1)].zprimus
	dataSF_NIP = dataSF[where(dataSF.IP eq 0)].zprimus
	dataQ_IP   = dataQ[where(dataQ.IP eq 1)].zprimus
	dataQ_NIP  = dataQ[where(dataQ.IP eq 0)].zprimus

	plothist, dataQ_IP, bin=binsize, color=cgColor('Red'), linestyle=0, thick=linethickness, yrange=[ymin,ymax], xtitle="redshift", xrange=[0.18,1.02], $
;	  title="Conservative Mass Cutoff"
	  title="Mass Limit - 0.5 dex"
	plothist, dataQ_NIP, bin=binsize, color=cgColor('Red'), linestyle=2, thick=linethickness, /overplot
	plothist, dataSF_IP, bin=binsize, color=cgColor('Blue'), linestyle=0, thick=linethickness, /overplot
	plothist, dataSF_NIP, bin=binsize, color=cgColor('Blue'), linestyle=2, thick=linethickness, /overplot
	
	LEGEND, ["ET IP", "ET non-IP", "LT IP", "LT non-IP"], colors=[cgColor('Red'),cgColor('Red'),cgColor('Blue'),cgColor('Blue')], $
	  linestyle=[0,2,0,2], box=0 ;, position=[0.7,7500]

        IF (string(outputFormat) eq 'ps') THEN BEGIN
                DEVICE, /CLOSE
        ENDIF
       	SET_PLOT, 'X'
END
