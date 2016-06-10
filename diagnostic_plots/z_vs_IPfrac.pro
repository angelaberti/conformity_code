PRO z_vs_IPfrac, outputFormat
        IF (string(outputFormat) eq 'ps') THEN BEGIN
		PS_OPEN, '~/figures/z_vs_IPfrac', THICK=5, /ENCAP
		DEVICE, /INCH, XS=8, YS=6
        ENDIF ELSE BEGIN
                SET_PLOT, 'X'
        ENDELSE

	binsize = 0.05

	A = FINDGEN(17) * (!PI*2/16.)
	USERSYM, 0.5*COS(A), 0.5*SIN(A), /FILL
	ptsymbol = 8

	zmin = 0.2
	zmax = 1.0

	ymin = 0.0
	ymax = 1.0
	
	data = mrdfits('~/conformity/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits', 1)

	dataSF = data[where(data.SFQ eq 1)]
	dataQ  = data[where(data.SFQ eq 0)]

	dataSF_IP  = dataSF[where(dataSF.IP eq 1)]
	dataSF_NIP = dataSF[where(dataSF.IP eq 0)]

	dataQ_IP   = dataQ[where(dataQ.IP eq 1)]
	dataQ_NIP  = dataQ[where(dataQ.IP eq 0)]
	
	plot, dataQ.zprimus, dataQ.mstar, psym=8, xrange=[zmin,zmax], yrange=[ymin,ymax], /nodata, $
	  xtitle="redshift", ytitle=textoidl("IP fraction")
	LEGEND, ["Late-Type", "Early-Type"], psym=[8,8], color=[cgColor('Blue'),cgColor('Red')], box=0

	IPfracSFzbinned = []
	IpfracQzbinned  = []
	zbins = []

	FOR z_index=0,(zmax-zmin)/binsize-1 DO BEGIN
	
		SFpoint = n_elements(dataSF[where( (dataSF.zprimus gt zmin+z_index*binsize) and $
					(dataSF.zprimus le zmin+(z_index+1)*binsize) )])
		Qpoint  = n_elements(dataQ[ where( (dataQ.zprimus  gt zmin+z_index*binsize) and $
					(dataQ.zprimus  le zmin+(z_index+1)*binsize) )])
		SF_IPpoint = n_elements(dataSF_IP[where( (dataSF_IP.zprimus gt zmin+z_index*binsize) and $
					      (dataSF_IP.zprimus le zmin+(z_index+1)*binsize) )])
		Q_IPpoint  = n_elements(dataQ_IP[where( (dataQ_IP.zprimus gt zmin+z_index*binsize) and $
					     (dataQ_IP.zprimus le zmin+(z_index+1)*binsize) )])

		IPfracSFpoint = float(SF_IPpoint)/float(SFpoint)
		IPfracQpoint  = float(Q_IPpoint)/float(Qpoint)

		zbin = mean([zmin+z_index*binsize, zmin+(z_index+1)*binsize])

		IPfracSFzbinned = [IPfracSFzbinned, IPfracSFpoint]
		IPfracQzbinned  = [IPfracQzbinned,  IPfracQpoint]

		zbins = [zbins, zbin]
		print, zbin, float(SFpoint), float(SF_IPpoint), float(Qpoint), float(Q_IPpoint), IPfracSFpoint, IPfracQpoint
	ENDFOR

	oplot, zbins, IPfracQzbinned, linestyle=0, color=cgColor('Red')
	oplot, zbins, IPfracSFzbinned, linestyle=0, color=cgColor('Blue')
	oplot, zbins, IPfracQzbinned, psym=ptsymbol, color=cgColor('Red')
	oplot, zbins, IPfracSFzbinned, psym=ptsymbol, color=cgColor('Blue')

        IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
       	SET_PLOT, 'X'
END
