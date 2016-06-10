PRO z_vs_mass_SFandIPstatus, outputFormat
        IF (string(outputFormat) eq 'ps') THEN BEGIN
		SET_PLOT, 'ps'
		DEVICE, file="../figures/z_vs_mass_SFandIPstatus.ps", /LANDSCAPE
        ENDIF ELSE BEGIN
                SET_PLOT, 'X'
        ENDELSE

	binsize = 0.02
	linethickness = 5

	A = FINDGEN(17) * (!PI*2/16.)
	USERSYM, 0.25*COS(A), 0.25*SIN(A), /FILL
	ptsymbol = 8

	zmin = 0.2
	zmax = 1.0

	ymin = 8
	ymax = 12
	
	data = mrdfits("../results/default_parameters/IP_data/zerodSFQ_IP_dv1500kms_dm0.5.fits", 1)

	dataSF = data[where(data.SFQ eq 1)]
	dataQ  = data[where(data.SFQ eq 0)]

	dataSF_IP  = dataSF[where(dataSF.IP eq 1)]
	dataSF_NIP = dataSF[where(dataSF.IP eq 0)]
	dataQ_IP   = dataQ[where(dataQ.IP eq 1)]
	dataQ_NIP  = dataQ[where(dataQ.IP eq 0)]
	
	dataSF_IPzbinned = []
	dataQ_IPzbinned  = []
	dataSF_NIPzbinned = []
	dataQ_NIPzbinned  = []
	zbins = []

	FOR z_index=0,(zmax-zmin)/binsize-1 DO BEGIN
		SF_IPpoint = mean(dataSF_IP[where( (dataSF_IP.zprimus gt zmin+z_index*binsize) and $
					     (dataSF_IP.zprimus le zmin+(z_index+1)*binsize) )].mstar)
		Q_IPpoint  = mean( dataQ_IP[where( (dataQ_IP.zprimus gt zmin+z_index*binsize) and $
					     (dataQ_IP.zprimus le zmin+(z_index+1)*binsize) )].mstar)
		SF_NIPpoint = mean(dataSF_NIP[where( (dataSF_NIP.zprimus gt zmin+z_index*binsize) and $
					     (dataSF_NIP.zprimus le zmin+(z_index+1)*binsize) )].mstar)
		Q_NIPpoint  = mean( dataQ_NIP[where( (dataQ_NIP.zprimus gt zmin+z_index*binsize) and $
					     (dataQ_NIP.zprimus le zmin+(z_index+1)*binsize) )].mstar)

		zbin = mean([zmin+z_index*binsize, zmin+(z_index+1)*binsize])

		dataSF_IPzbinned  = [dataSF_IPzbinned,  SF_IPpoint]
		dataQ_IPzbinned   = [dataQ_IPzbinned,   Q_IPpoint]
		dataSF_NIPzbinned = [dataSF_NIPzbinned, SF_NIPpoint]
		dataQ_NIPzbinned  = [dataQ_NIPzbinned,  Q_NIPpoint]

		zbins = [zbins, zbin]
		print, zbin, SF_IPpoint, SF_NIPpoint, Q_IPpoint, Q_NIPpoint
	ENDFOR

	FOR i=0,4 DO BEGIN
		plot, dataQ.zprimus, dataQ.mstar, psym=ptsymbol, xrange=[zmin,zmax], yrange=[ymin,ymax], /nodata, xtitle="redshift", ytitle=textoidl("log (M/M_{\odot})")
		IF i eq 0 THEN $
			LEGEND, ["Early-Type IP", "Early-Type non-IP", "Late-Type IP", "Late-Type non-IP"], linestyle=[0,2,0,2], $
			  color=[cgColor('Red'),cgColor('Red'),cgColor('Blue'),cgColor('Blue')], box=0, position=[0.6,9], $
;			  color=cgColor('Black'), box=0, position=[0.6,9], $
			  thick=linethickness ELSE $
		IF i eq 1 THEN $
		  oplot, dataQ_IP.zprimus, dataQ_IP.mstar, psym=ptsymbol, color=cgColor('Red') ELSE $
		IF i eq 2 THEN $
		  oplot, dataQ_NIP.zprimus, dataQ_NIP.mstar, psym=ptsymbol, color=cgColor('Orange') ELSE $
		IF i eq 3 THEN $
		  oplot, dataSF_IP.zprimus, dataSF_IP.mstar, psym=ptsymbol, color=cgColor('Blue') ELSE $
		IF i eq 4 THEN BEGIN
		  oplot, dataSF_NIP.zprimus, dataSF_NIP.mstar, psym=ptsymbol, color=cgColor('Cyan')
		ENDIF

		oplot, zbins, dataQ_IPzbinned, linestyle=0, thick=linethickness, color=cgColor('Red')
		oplot, zbins, dataQ_NIPzbinned, linestyle=2, thick=linethickness, color=cgColor('Red')
		oplot, zbins, dataSF_IPzbinned, linestyle=0, thick=linethickness, color=cgColor('Blue')
		oplot, zbins, dataSF_NIPzbinned, linestyle=2, thick=linethickness, color=cgColor('Blue')
			
;		IF i ne 0 THEN BEGIN
;			LEGEND, ["Early-Type IP", "Early-Type non-IP", "Late-Type IP", "Late-Type non-IP"], linestyle=[0,2,0,2], $
;			  color=[cgColor('Red'),cgColor('Red'),cgColor('Blue'),cgColor('Blue')], box=0, position=[0.6,9], $
;			  color=cgColor('Black'), box=0, position=[0.6,9], $
;			  thick=linethickness
;		ENDIF
	ENDFOR

        IF (string(outputFormat) eq 'ps') THEN BEGIN
                DEVICE, /CLOSE
        ENDIF
       	SET_PLOT, 'X'
END
