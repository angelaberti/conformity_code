PRO z_vs_latefrac, outputFormat
        IF (string(outputFormat) EQ 'ps') THEN BEGIN
		PS_OPEN, '~/figures/z_vs_latefrac', THICK=5, /ENCAP
		DEVICE, /INCH, XS=8, YS=6
        ENDIF ELSE BEGIN
                SET_PLOT, 'X'
        ENDELSE

	binsize = 0.1

	A = FINDGEN(17) * (!PI*2/16.)
	USERSYM, 0.5*COS(A), 0.5*SIN(A), /FILL
	ptsymbol = 8

	!p.multi=0.5
	!p.charsize=1.5

	zmin = 0.2
	zmax = 1.0

	ymin = 0.0
	ymax = 1.05
	
	dataAll = mrdfits('~/conformity/results/zerodSFQ_all_cart.fits', 1)
	dataAllcomp = mrdfits('~/conformity/results/conservative_mass_cutoff/allAboveMassCompLim-0.0.fits', 1)
	dataIP	= mrdfits('~/conformity/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits', 1)
	dataIP = dataIP[WHERE(dataIP.IP EQ 1)]

	PLOT, findgen(10), findgen(10), xrange=[0.2,1.0], yrange=[ymin,ymax], /NODATA, xtitle='Redshift', ytitle='Late-type fraction';, title='Late-type fraction vs. redshift for entire sample'
	
	colors = ['Magenta', 'Red', 'Blue', 'green', 'purple']

;	FOR f=0,4 DO BEGIN
;		data = dataAll[where(dataAll.field EQ fields[f])]
		data = dataAll
		dataComp = dataAllcomp

		fracSF	= []
		zbins	= []

		fracSF_weight	= []
		fracSFcomp_weight = []
		IPfracSF_weight	= []
		IPfracSF_low 	= []
		IPfracSF_high	= []
		IPfracSF_weight_low 	= []
		IPfracSF_weight_high	= []

		FOR z_index=0,(zmax-zmin)/binsize-1 DO BEGIN
			zbin = mean([zmin+z_index*binsize, zmin+(z_index+1)*binsize])

			dataAllzbin	= data[where( (data.zprimus gt zmin+z_index*binsize) AND (data.zprimus le zmin+(z_index+1)*binsize) )]
			dataCompzbin	= dataComp[where( (dataComp.zprimus gt zmin+z_index*binsize) AND (dataComp.zprimus le zmin+(z_index+1)*binsize) )]
			dataIPzbin	= dataIP[where( (dataIP.zprimus gt zmin+z_index*binsize) AND (dataIP.zprimus le zmin+(z_index+1)*binsize) )]
			dataIPlow_zbin	= dataIPzbin[WHERE((dataIPzbin.mstar GE 10.) AND (dataIPzbin.mstar LT 10.5), /NULL)]
			dataIPhigh_zbin	= dataIPzbin[WHERE((dataIPzbin.mstar GE 10.5) AND (dataIPzbin.mstar LT 11.), /NULL)]

			PRINT, z_index, n_elements(dataIPlow_zbin), n_elements(dataIPhigh_zbin)

			nAllzbin 	= float(n_elements(dataAllzbin))

			fracSFzbin	= n_elements(dataAllzbin[where(dataAllzbin.SFQ EQ 1)])/nAllzbin
			fracSFzbin_wght	= TOTAL(dataAllzbin[where(dataAllzbin.SFQ EQ 1)].targ_weight)/TOTAL(dataAllzbin.targ_weight)
			fracSFcompzbin_wght = TOTAL(dataCompzbin[where(dataCompzbin.SFQ EQ 1)].targ_weight)/TOTAL(dataCompzbin.targ_weight)
			IPfracSFzbin_wght = TOTAL(dataIPzbin[where(dataIPzbin.SFQ EQ 1)].targ_weight)/TOTAL(dataIPzbin.targ_weight)
			IPfracSFzbin_low  = n_elements(dataIPlow_zbin[where(dataIPlow_zbin.SFQ EQ 1, /NULL)])/FLOAT(n_elements(dataIPlow_zbin))
			IPfracSFzbin_high = n_elements(dataIPhigh_zbin[where(dataIPhigh_zbin.SFQ EQ 1, /NULL)])/FLOAT(n_elements(dataIPhigh_zbin))
			IPfracSFzbin_wght_low  = TOTAL(dataIPlow_zbin[where(dataIPlow_zbin.SFQ EQ 1, /NULL)].targ_weight)/TOTAL(dataIPlow_zbin.targ_weight)
			IPfracSFzbin_wght_high = TOTAL(dataIPhigh_zbin[where(dataIPhigh_zbin.SFQ EQ 1, /NULL)].targ_weight)/TOTAL(dataIPhigh_zbin.targ_weight)

			fracSF		= [fracSF, fracSFzbin]
			fracSF_weight	= [fracSF_weight, fracSFzbin_wght]
			fracSFcomp_weight = [fracSFcomp_weight, fracSFcompzbin_wght]
			IPfracSF_weight	= [IPfracSF_weight, IPfracSFzbin_wght]
			IPfracSF_low	= [IPfracSF_low, IPfracSFzbin_low]
			IPfracSF_high	= [IPfracSF_high, IPfracSFzbin_high]
			IPfracSF_weight_low	= [IPfracSF_weight_low, IPfracSFzbin_wght_low]
			IPfracSF_weight_high	= [IPfracSF_weight_high, IPfracSFzbin_wght_high]

			zbins = [zbins, zbin]
		ENDFOR

;		OPLOT, zbins, fracSF, COLOR=cgColor(colors[f]), LINESTYLE=0
		OPLOT, zbins, fracSF_weight, LINESTYLE=3, color=cgColor(colors[0])
		PRINT, zbins, fracSF_weight;, LINESTYLE=3, color=cgColor(colors[0])
		OPLOT, zbins, IPfracSF_weight, LINESTYLE=0, color=cgColor(colors[2])
		OPLOT, zbins, fracSFcomp_weight, LINESTYLE=2, color=cgColor(colors[1])
		OPLOT, zbins, IPfracSF_weight_low, LINESTYLE=4, color=cgColor(colors[3])
		OPLOT, zbins, IPfracSF_weight_high, LINESTYLE=5, color=cgColor(colors[4])
;		OPLOT, zbins, IPfracSF_low, LINESTYLE=1, color=cgColor(colors[3])
;		OPLOT, zbins, IPfracSF_high, LINESTYLE=1, color=cgColor(colors[4])

		LEGEND, ['All PRIMUS', 'All above M13 completeness limit', 'IP only (also above M13 limit)', textoidl('10 < log(M_{*IP}/M_{\odot}) < 10.5'), textoidl('10.5 < log(M_{*IP}/M_{\odot}) < 11')], $
		COLOR=cgColor(colors), LINESTYLE=[3,2,0,4,5], BOX=0, /BOTTOM, /LEFT
;	END

        IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
       	SET_PLOT, 'X'
END
