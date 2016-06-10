; Q1 lowest masses, Q4 highest masses

PRO massDiv, outputFormat
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '../../figures/massCompLimits', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	; default parameters entire sample and IP only
	data	= mrdfits('../default_parameters/IP_data/zerodSFQ_IP_dz1.0_dm0.5.fits',1)
	dataIP	= data[where(data.IP eq 1)]

	; all sources above default completeness limit
	dataAbove = mrdfits('../allAboveMassCompLim-0.5.fits',1)

	; initialize
	dataIP_Q1 = []
	dataIP_Q2 = []
	dataIP_Q3 = []
	dataIP_Q4 = []

	FOR i=0,79 DO BEGIN
		zmin = 0.2+0.01*float(i)
		zmax = 0.2+0.01*(float(i)+1.)
		dataIP_zcut	= dataIP[where( (dataIP.zprimus ge zmin) AND (dataIP.zprimus le zmax) )]
		dataAbove_zcut	= dataAbove[where( (dataAbove.zprimus ge zmin) AND (dataAbove.zprimus le zmax) )]

;		PRINT, zmin, zmax, " All above: ", n_elements(dataAbove_zcut), " IP: ", n_elements(dataIP_zcut)

		; median mass of all sources above completeness limit (divides Q2 and Q3)
		m23 = median(dataAbove_zcut.mstar)

		; median mass of lower half (divides Q1 and Q2)
		dataAbove_zcutLow  = dataAbove_zcut[where(dataAbove_zcut.mstar lt m23)]
		m12 = median(dataAbove_zcutLow.mstar)

		; median mass of upper half (divides Q3 and Q4)
		dataAbove_zcutHigh = dataAbove_zcut[where(dataAbove_zcut.mstar ge m23)]
		m34 = median(dataAbove_zcutHigh.mstar)

		dataIP_zcutQ1 	= dataIP_zcut[where( dataIP_zcut.mstar le m12 )]
		dataIP_zcutQ2	= dataIP_zcut[where( (dataIP_zcut.mstar ge m12) AND (dataIP_zcut.mstar le m23) )]
		dataIP_zcutQ3	= dataIP_zcut[where( (dataIP_zcut.mstar ge m23) AND (dataIP_zcut.mstar le m34) )]
		dataIP_zcutQ4 	= dataIP_zcut[where( dataIP_zcut.mstar ge m34 )]

		dataIP_Q1 = [dataIP_Q1, dataIP_zcutQ1]
		dataIP_Q2 = [dataIP_Q2, dataIP_zcutQ2]
		dataIP_Q3 = [dataIP_Q3, dataIP_zcutQ3]
		dataIP_Q4 = [dataIP_Q4, dataIP_zcutQ4]
	ENDFOR

        mass_binwidth	= 0.2
        z_binwidth 	= 0.04
        mass	= dataAbove[where((dataAbove.zprimus ge 0.2) AND (dataAbove.zprimus le 1.0))].mstar
        z	= dataAbove[where((dataAbove.zprimus ge 0.2) AND (dataAbove.zprimus le 1.0))].zprimus

        n_massbins = ceil(4./mass_binwidth)
        n_zbins	= floor(0.8/z_binwidth)

        hist2d = []
        FOR s_index=0,n_zbins-1 DO BEGIN
                z_low  = 0.2 + z_binwidth*float(s_index)
                z_high = 0.2 + z_binwidth*float(s_index+1)

                in_z_range = dataAbove[where((dataAbove.zprimus ge z_low) AND (dataAbove.zprimus le z_high), /null)]
;		PRINT, "z_low:", z_low, " z_high:", z_high, " # in range:", n_elements(in_z_range)

                z_row = []
                FOR m_index=0,n_massbins-1 DO BEGIN
                        m_low   = 8. + mass_binwidth*float(m_index)
                        m_high  = 8. + mass_binwidth*float(m_index+1)

                        IF in_z_range ne !null THEN $
                        in_unit = n_elements(in_z_range[where((in_z_range.mstar ge m_low) AND (in_z_range.mstar le m_high), /null)]) $
                        ELSE in_unit = 0

                        z_row = [z_row, in_unit]
                ENDFOR

                hist2d = [[hist2d], [z_row]]
        ENDFOR

        masspts = 8. + mass_binwidth*(findgen(n_massbins)+1.)
        zpts 	= 0.2 + z_binwidth*(findgen(n_zbins)+1.)

	print, n_massbins, n_zbins, n_elements(masspts), n_elements(zpts)

;	PLOTHIST, dataIP.mstar, bin=binsize, yrange=[0,500], xtitle=textoidl('Log (M_{IP} / M_{\odot})'), ytitle='Number', color=cgColor('White'), /NODATA

	!p.charsize=1.4
	PLOT, dataIP.zprimus, dataIP.mstar, psym=3, xrange=[0.2,1.0], yrange=[8.25,11.75], xtitle='Redshift', ytitle=textoidl('Log (M_{IP} / M_{\odot})'), title='Mass Completeness Limits for IP Candidates', /NODATA

	dataList = [dataIP_Q1, dataIP_Q2, dataIP_Q3, dataIP_Q4]
	colorList = ['Red', 'Blue', 'Green', 'Magenta']

;	FOR i=0,3 DO BEGIN
;		IF (i eq 0) THEN pdata=dataIP_Q1 ELSE $
;		IF (i eq 1) THEN pdata=dataIP_Q2 ELSE $
;		IF (i eq 2) THEN pdata=dataIP_Q3 ELSE $
;		IF (i eq 3) THEN pdata=dataIP_Q4

;		OPLOT, pdata.zprimus, pdata.mstar, psym=3, color=cgColor(colorList[i])

;		PLOTHIST, pdata[where(pdata.SFQ eq 0)].mstar, bin=binsize, linestyle=1, color=cgColor(colorList[i]), /OVERPLOT
;		PLOTHIST, pdata[where(pdata.SFQ eq 1)].mstar, bin=binsize, linestyle=0, color=cgColor(colorList[i]), /OVERPLOT
;	ENDFOR
	levels = float(max(transpose(hist2d)))*findgen(10)/9.
	PRINT, levels
	CONTOUR, transpose(hist2d), zpts-z_binwidth/2, masspts, levels=levels, color=cgColor('Gray'), /OVERPLOT
;	CONTOUR, min_curve_surf(transpose(hist2d)), levels=levels, color=cgColor('Gray'), /OVERPLOT

        zbins = [0.20, 0.30, 0.40, 0.50, 0.65, 0.80, 1.00]

        zbinsMidpts = [0.2]
        FOR i=1,n_elements(zbins)-3 DO BEGIN
                zbinsMidpts = [zbinsMidpts, [mean([zbins[i],zbins[i+1]])]]
        ENDFOR
	zbinsMidpts = [zbinsMidpts,[1.0]]

        cons_mean = [9.172,9.492,9.792,10.138,10.508,10.886]

	colors=[cgColor('Blue'), cgColor('Red'), cgColor('Black')]

        OPLOT, zbinsMidpts, cons_mean-0.5, linestyle=0, color=colors[0]
	OPLOT, zbinsMidpts, cons_mean, linestyle=0, color=colors[1]
        OPLOT, zbinsMidpts, 0*indgen(7)+10, linestyle=2, color=colors[2]

	LEGEND, ["Default", "Conservative", textoidl('Constant (M_{IP} > 10^{10} M_{\odot})')], linestyle=[0,0,2], color=colors, box=0, /BOTTOM, /RIGHT

;	mwrfits, dataIP_Q1, 'IP_data/dataIP_Q1.fits', /create
;	mwrfits, dataIP_Q2, 'IP_data/dataIP_Q2.fits', /create
;	mwrfits, dataIP_Q3, 'IP_data/dataIP_Q3.fits', /create
;	mwrfits, dataIP_Q4, 'IP_data/dataIP_Q4.fits', /create

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
