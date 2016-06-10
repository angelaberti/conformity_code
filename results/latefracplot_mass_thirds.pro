;PRO latefracplot, dataSet, binning, outputFormat
PRO latefracplot_mass_thirds, dataSet, binning, input_xmax, outputFormat

; data sets 	
;	4 mass thirds (0.8 dex)		dz=2
;	5 mass thirds (1.0 dex)		dz=2	
; binning
;	0 zall (default)
;	1 zbins
; outputFormat:	'ps' OR anything else	

 specs={paramTypeDir:  ['variable_mass+0.5dex',		$; 0
			'conservative_mass_cutoff',	$; 1
			'conservative+0.3dex',		$; 2
			'single_mass_cutoff',		$; 3
			'mass_bins',			$; 4
		 	'mass_bins',			$; 5
			'variable_mass+0.8dex', 	$; 7
			'variable_mass+1.0dex'], 	$; 8

	IPdataDir:     ['variable_mass+0.5dex',		$; 0
			'conservative_mass_cutoff',	$; 1
			'conservative+0.3dex',		$; 2
			'single_mass_cutoff',		$; 3
			'mass_bins',			$; 4
			'mass_bins',			$; 5
			'variable_mass+0.8dex', 	$; 7
			'variable_mass+1.0dex'], 	$; 8

	paramType:     ['variable mass limit + 0.5 dex',	$; 0
		 	'conservative IP mass cutoff',		$; 1
			'conservative IP mass cutoff + 0.3 dex',$; 2
			'single IP mass cutoff',		$; 3
			'IP mass thirds (0.8 dex diff)',	$; 4
			'IP mass thirds (1.0 dex diff)',	$; 5
			'variable mass limit + 0.8 dex',	$; 7
			'variable mass limit + 1.0 dex'],	$; 8

	latefracSuffix:['zerodSFQ_IP_dz2.0',		$
			'zerodSFQ_IP_dz2.0_dm0.0',	$
			'zerodSFQ_IP_dz2.0_dm0.3',	$
			'zerodSFQ_IP_dz1.0_singleMass',	$
			'dataIP_T',			$
			'dataIP_T',			$
			'zerodSFQ_IP_dz2.0',		$
			'zerodSFQ_IP_dz2.0'],		$

	iterateOver:   ['z',	'z',	'z',	'z',	'thirds', 'thirds', 'z', 'z'], $

	figNameSuffix: ['variable_mass_0.5dex',		$	
			'conservative_mass_dz2.0',	$
			'conservative+0.3dex_dz2.0',	$
			'single_mass', 			$
			'mass_bins_0.8dex',		$
		 	'mass_bins_1.0dex',		$
			'variable_mass_0.8dex',		$
			'variable_mass_1.0dex'], 	$	

	title:	       [textoidl('Moustakas13 Mass Limit + 0.5 (0.0) dex for SF (Q) IP (\Deltaz=2.0)'),	$
			textoidl('Moustakas13 IP Mass Limit (\Deltaz = 2.0)'), $
			textoidl('Moustakas13 IP Mass Limit + 0.3 dex (\Deltaz=2.0)'), $
			textoidl('Log (M_{IP} / M_{\odot}) > 10'), $
			textoidl('IP Stellar Mass Thirds; variable IP mass diff. 0.8 dex; \Deltaz=2.0'), $
			textoidl('IP Stellar Mass Thirds; variable IP mass diff. 1.0 dex; \Deltaz=2.0'), $
			textoidl('Moustakas13 Mass Limit + 0.8 (0.0) dex for SF (Q) IP (\Deltaz=2.0)'),	$
			textoidl('Moustakas13 Mass Limit + 1.0 (0.0) dex for SF (Q) IP (\Deltaz=2.0)')]}

	IF (binning eq 1) THEN BEGIN
		binSuffix = '_zbins'
		binType   = 'z bins'
 	ENDIF ELSE BEGIN
		binSuffix = '_zall'
		binType   = 'none'
	ENDELSE
	IF specs.paramTypeDir[dataSet] eq 'mass_bins' THEN BEGIN
		binSuffix = ''
		binType   = 'IP mass thirds'
	ENDIF

	PRINT, 'Data selected: ' + specs.paramType[dataSet]
	PRINT, 'Binning: ' + binType

	IPdataDir	= '~/conformity/results/' + specs.IPdataDir[dataSet] + '/IP_data/'
	latefracdataDir = '~/conformity/results/' + specs.paramTypeDir[dataSet] + '/latefrac_data_hist/'
	latefracSuffix 	= specs.latefracSuffix[dataSet]

	IF (specs.figNameSuffix[dataSet] EQ 'mass_bins_0.8dex') THEN BEGIN
		IPdataDir 	= IPdataDir + 'variable_mass+0.8dex/'
		latefracdataDir = latefracdataDir + 'variable_mass+0.8dex/'
	ENDIF ELSE BEGIN
	IF (specs.figNameSuffix[dataSet] EQ 'mass_bins_1.0dex') THEN BEGIN
		IPdataDir 	= IPdataDir + 'variable_mass+1.0dex/'
		latefracdataDir = latefracdataDir + 'variable_mass+1.0dex/'
	ENDIF
	ENDELSE
	
	outputFileName = 'latefracplot_' + specs.figNameSuffix[dataSet] + binSuffix
	outputFileLocation = '~/figures/' + outputFileName

	ymin	= 0.35
	ymax	= 1.0
;	IF (specs.paramTypeDir[dataSet] eq 'dR_60Mpc_largeBins') OR (specs.paramTypeDir[dataSet] eq 'mass_bins_largeBins') THEN BEGIN
;		ymin=0.625
;		ymax=0.825
;	ENDIF	
	yr	= ymax - ymin

	!p.multi	= 0
	!p.charsize	= 1.15
	xyoutscharsize	= 1.5
	
	nx 	= 1.
	ny 	= 1.
	zArray 	= [0.2, 1.0]
	zbinwidth = 0.8
	xm 	= 0.08
;	IF dataSet eq 4 THEN xm = 0.10
;	IF (specs.paramTypeDir[dataSet] eq 'dR_60Mpc_largeBins') OR (specs.paramTypeDir[dataSet] eq 'mass_bins_largeBins') THEN xm = 0.10
	ym 	= 0.09

	IF (binType eq 'IP mass thirds') THEN BEGIN
		!p.multi = 3
		!p.charsize = 2
		xyoutscharsize = 1.1

	        nx 	= 1.
	        ny 	= 3.
		zArray  = [0.6, 0.8, 0.2, 0.4, 1.0]
		zbinwidth = 0.2
	        xm	= 0.08
	        ym	= 0.08
	ENDIF

	IF specs.iterateOver[dataSet] eq 'z' THEN iterateList=zArray ELSE $
	IF specs.iterateOver[dataSet] eq 'thirds' THEN iterateList=['3','2','1','','']	

        dx = (1 - 1.5*xm)/nx
        dy = (1 - 2.0*ym)/ny

	IF (string(outputFormat) eq 'ps') THEN BEGIN
		PS_OPEN, outputFileLocation, THICK=5, /ENCAP
		DEVICE, /INCH, XS=8.5, YS=11
	ENDIF ELSE BEGIN
		SET_PLOT, 'X'
	ENDELSE

	ERASE
	FOR i=0,n_elements(iterateList)-2 DO BEGIN
                IF (i le (nx-1)) THEN BEGIN
                        xtitle = 'Projected Radius (Mpc)'
                        xtickformat = ''
                ENDIF ELSE BEGIN
                        xtitle = ''
                        xtickformat = '(A1)'
                ENDELSE
                IF (i mod nx eq 0) THEN BEGIN
                        ytitle = 'Late-type Fraction'
                        ytickformat = ''
                ENDIF ELSE BEGIN
                        ytitle = ''
                        ytickformat = '(A1)'
                ENDELSE

		IF (i eq 2) THEN title = specs.title[dataSet] ELSE title = ''

		; define each plot position
		pos = [ xm + dx*float(i mod nx), ym + dy*float(floor(float(i)/nx)), $
			xm + dx*(1. + float(i mod nx)), ym + dy*(1. + float(floor(float(i)/nx))) ]

		IF specs.iterateOver[dataSet] eq 'z' THEN BEGIN
			zlow 	= strtrim(string(iterateList[i], format='(f20.1)'),1)
			zhigh	= strtrim(string(iterateList[i] + zbinwidth, format='(f20.1)'),1)
			third = ''
		ENDIF
		IF specs.iterateOver[dataSet] eq 'thirds' THEN BEGIN
			zlow 	= '0.2'
			zhigh	= '1.0'
			third = iterateList[i]
		ENDIF

		latefracInputFile = 'latefrac_' + zlow + '_' + zhigh + '_' + latefracSuffix + third + '.fits'
		data 	= mrdfits(latefracdataDir + latefracInputFile, 1)

		IPdataInputFile = latefracSuffix + third + '.fits'

		IF (specs.paramTypeDir[dataSet] eq 'mass_bins') THEN $
;		  IPdataInputFile = 'zerodSFQ_IP_dz2.0.fits'
		  IPdataInputFile = 'dataIP_T' + third + '.fits'
		IF (specs.paramTypeDir[dataSet] eq 'variable_mass+0.5dex') OR $
		   (specs.paramTypeDir[dataSet] eq 'variable_mass+0.8dex') OR $
		   (specs.paramTypeDir[dataSet] eq 'variable_mass+1.0dex') THEN IPdataInputFile = 'zerodSFQ_IP_dz2.0.fits'

		dataAll  = mrdfits(IPdataDir + IPdataInputFile, 1)

		PRINT, 'Late fraction data: ', latefracdataDir + latefracInputFile
		PRINT, 'IP data: ', IPdataDir + IPdataInputFile		

		xmin	= 0.0
		xmax	= max(data.Rmax) + min(data.Rmax)
		xmax 	= input_xmax
		xr 	= xmax - xmin

		frac_IPSF = data.n_late_IPSF/data.n_tot_IPSF
		frac_IPQ  = data.n_late_IPQ/data.n_tot_IPQ

		dfrac_IPSF = poissonError(data.n_late_IPSF, data.n_tot_IPSF)
		dfrac_IPQ  = poissonError(data.n_late_IPQ, data.n_tot_IPQ)

		PRINT, data.subtotal_neighbors_IPSF
		PRINT, data.subtotal_neighbors_IPQ
		Nneighbors = bigint(max(data.subtotal_neighbors_IPSF) + max(data.subtotal_neighbors_IPQ))
		
		dataAll_zcut = dataAll[where((dataAll.zprimus ge float(zlow)) and (dataAll.zprimus le float(zhigh)))]
		Ntotal	= bigint(n_elements(dataAll_zcut))
		NIPint 	= n_elements(dataAll_zcut[where(dataAll_zcut.IP eq 1)])
		NIP 	= bigint(NIPint)
		IPfrac  = string(100.*NIPint/float(n_elements(dataAll_zcut)),format='(f20.1)')

		medianMassAll	= median(dataAll_zcut.mstar)
		medianIPmass	= median(dataAll_zcut[where(dataAll_zcut.IP eq 1)].mstar)
		meanMassAll	= mean(dataAll_zcut.mstar)
		dataIP		= dataAll_zcut[where(dataAll_zcut.IP eq 1)]

		medianIPmass	= median(dataIP.mstar)
		medianIPmassSF	= median(dataIP[where(dataIP.SFQ eq 1)].mstar)
		medianIPmassQ	= median(dataIP[where(dataIP.SFQ eq 0)].mstar)
		
		plot, data.Rmax, frac_IPSF, xrange=[xmin,xmax], yrange=[ymin,ymax], xtickformat=xtickformat, ytickformat=ytickformat, $
		  xtitle=xtitle, ytitle=ytitle, title=title, position=pos
		errplot, data.Rmax, frac_IPSF-dfrac_IPSF, frac_IPSF+dfrac_IPSF

		oplot, data.Rmax, frac_IPQ, linestyle=2
		errplot, data.Rmax, frac_IPQ-dfrac_IPQ, frac_IPQ+dfrac_IPQ

		IF i eq nx*ny-1 THEN LEGEND, ["Late-type IP","Early-type IP"], linestyle=[0,2], box=0, /TOP, /RIGHT, charsize=xyoutscharsize

		statsFile = latefracdataDir + latefracInputFile

		sigmaRange000_075 = sigmaRange(statsFile,0,2)
;		sigmaRange000_200 = sigmaRange(statsFile,0,7)
;		sigmaRange000_300 = sigmaRange(statsFile,0,11)
;		sigmaRange000_400 = sigmaRange(statsFile,0,15)

		sigmaRange075_200 = sigmaRange(statsFile,3,7)
;		sigmaRange075_300 = sigmaRange(statsFile,3,11)
		sigmaRange075_400 = sigmaRange(statsFile,3,15)

                xyouts, xmin+0.05*xr, ymin+0.13*yr, zlow + ' < z < ' + zhigh, ALIGNMENT=0.0, charsize=xyoutscharsize

                xyouts, xmin+0.95*xr, ymin+0.45*yr, textoidl("\sigma_{R < 0.75 Mpc}: ") + strtrim(string(sigmaRange000_075,format='(f20.2)'),1), ALIGNMENT=1.0, charsize=xyoutscharsize
                xyouts, xmin+0.95*xr, ymin+0.37*yr, textoidl("\sigma_{0.75 < R < 2 Mpc}: ") + strtrim(string(sigmaRange075_200,format='(f20.2)'),1), ALIGNMENT=1.0, charsize=xyoutscharsize
                xyouts, xmin+0.95*xr, ymin+0.29*yr, textoidl("\sigma_{0.75 < R < 4 Mpc}: ") + strtrim(string(sigmaRange075_400,format='(f20.2)'),1), ALIGNMENT=1.0, charsize=xyoutscharsize

;               xyouts, xmin+0.95*xr, ymin+0.21*yr, textoidl("Total galaxies: ") + strtrim(Ntotal, 1), ALIGNMENT=1.0, charsize=xyoutscharsize
;		xyouts, xmin+0.95*xr, ymin+0.21*yr, textoidl("N_{neighbors}: ") + Nneighbors, ALIGNMENT=1.0, charsize=xyoutscharsize
;		xyouts, xmin+0.95*xr, ymin+0.13*yr, textoidl("N_{IP}: ") + NIP + " ("+strtrim(IPfrac, 1)+"%)", ALIGNMENT=1.0, charsize=xyoutscharsize
                xyouts, xmin+0.95*xr, ymin+0.05*yr, textoidl("Median M_{IP}: ") + strtrim(string(medianIPmassSF, format='(f20.2)'), 1) + ' (SF) ' $
										+ strtrim(string(medianIPmassQ, format='(f20.2)'), 1) + ' (Q)', ALIGNMENT=1.0, charsize=xyoutscharsize
	ENDFOR

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

