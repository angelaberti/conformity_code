; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; latefrac:	calculates late-type fractions
; BSE:		uses bootstrap error method

PRO latefracplot_BSE_allz_v02, outputFormat;, zmin, zmax;, dz_coeff, printEvery
	!P.FONT=0

	SFcolor = cgcolor('blue')
	Qcolor  = cgcolor('red')

  IPdatafiles = [ $
;	'~/results/default_parameters/IP_data/zerodSFQ_IP_dz2.0_dm0.5.fits', $
	'~/results/single_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_singleMass.fits', $
	'~/results/match_IP_sample_rigorous/unmatchedExample/matchMassOnly/matchedIPsampleMassOnlyFBF_PHI3.7.fits', $
;	'~/results/conservative+0.3dex/IP_data/zerodSFQ_IP_dz2.0_dm0.3.fits', $
;	'~/results/variable_mass+1.0dex/IP_data/zerodSFQ_IP_dz2.0.fits', $
	'~/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', $
	'~/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits']

  titles = [$;'Default Parameters', $
	'(c) Matched Redshift Only', $
;	'Conservative Mass Cutoff +0.3 dex', $
;	'Variable Mass Cutoff (+1.0 dex for SF)', $
	textoidl('(d) Matched M_{*} Only'), $
	textoidl('(a) Matched M_{*} & Redshift'), $
	'(b) No Matching']

;  tags = ['default', 'singleMass', 'var1.0dex']
;  tags = ['default', 'singleMass', 'conservative', 'conservative0.3dex', 'var1.0dex']
  tags = ['singleMass', 'matchMassOnly', '', 'conservative']

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	Rmax 		= 15.
	dRproj		= 1.
	n_annuli 	= CEIL(Rmax/dRproj)
	printEvery	= 100

        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

	ERASE
;	!P.MULTI=[4,2,2]
	!P.MULTI=4
	!P.CHARSIZE=2.6
	charsz=1.4

	Rproj_array = FLOAT(dRproj)*(findgen(n_annuli+1))	
;	PRINT, Rproj_array

	xmin = 0
	xmax = n_annuli*dRproj
	xr = xmax-xmin

	ymin = 0.71
	ymax = 0.88
	yr = ymax-ymin

	z_low = 0.2
	z_high = 1.0

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/latex/figures/unmatchedIPsampleCompare', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=8, YS=8, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

	dataAll_allz = MRDFITS('~/results/zerodSFQ_all_cart.fits', 1)
	; eliminate data with targ_weight < 1
        dataAll_allz = dataAll_allz[where(dataAll_allz.targ_weight GE 1.)]

	nx=2.
	ny=2.

	xm=0.10
	ym=0.10

	dx = (1 - 1.5*xm)/nx
	dy = (1 - 2.0*ym)/ny

  FOR i=0,nx*ny-1 DO BEGIN
        IF (i LE (nx-1)) THEN BEGIN
                xtitle = 'Projected Radius (Mpc)'
                xtickformat = ''
        ENDIF ELSE BEGIN
                xtitle = ''
                xtickformat = '(A1)'
        ENDELSE
        IF (i mod nx EQ 0) THEN BEGIN
                ytitle = 'Late-type Fraction'
                ytickformat = ''
        ENDIF ELSE BEGIN
                ytitle = ''
                ytickformat = '(A1)'
        ENDELSE

        pos = [ xm + dx*FLOAT(i mod nx), ym + dy*FLOAT(FLOOR(FLOAT(i)/nx)), $
                xm + dx*(1. + FLOAT(i mod nx)), ym + dy*(1. + FLOAT(FLOOR(FLOAT(i)/nx))) ]
;	PRINT, pos

;  FOR i=0,n_elements(IPdatafiles)-1 DO BEGIN
	PRINT, ''
	PRINT, IPdatafiles[i]
	IPstats, IPdatafiles[i]

	SFcolor = cgcolor('blue')
	Qcolor  = cgcolor('red')

;	IF (i EQ 5) THEN dataIP_allz = MRDFITS('~/results/match_IP_sample_rigorous/latefrac_allz_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE.fits', 1) ELSE $
;	IF (i EQ 5) THEN dataIP_allz = MRDFITS('~/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1) ELSE $
	  dataIP_allz  = MRDFITS(IPdatafiles[i], 1)

	; eliminate data with targ_weight < 1
        dataIP_allz  = dataIP_allz[where((dataIP_allz.targ_weight GE 1.) AND (dataIP_allz.IP EQ 1))]

	dataAll = dataAll_allz
	dataIP = dataIP_allz

	zlabel = 'z=[' + decimal(z_low,2) + ', ' + decimal(z_high,2) + ']'

	IF (i EQ 2) THEN data = MRDFITS('~/results/match_IP_sample_rigorous/latefrac_allz_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE.fits', 1) ELSE $
	  data = MRDFITS('~/results/match_IP_sample_rigorous/unmatchedExample/latefrac_' + tags[i] + '_allz_dR1Mpc_BSE.fits', 1)

	n_tot_IPSF  = data.n_tot_IPSF
	n_late_IPSF = data.n_late_IPSF
	errors_IPSF = data.errors_IPSF
	frac_IPSF   = n_late_IPSF/n_tot_IPSF

	n_tot_IPQ  = data.n_tot_IPQ
	n_late_IPQ = data.n_late_IPQ
	errors_IPQ = data.errors_IPQ
	frac_IPQ   = n_late_IPQ/n_tot_IPQ

;	PLOT, Rproj_array+0.5*dRproj, frac_IPSF, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', POSITION=pos, $
	PLOT, Rproj_array+0.5*dRproj, frac_IPSF, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle=xtitle, ytitle=ytitle, xtickformat=xtickformat, ytickformat=ytickformat, POSITION=pos, $
	  /NODATA
;	  title = 'Field-by-Field Matched SF & Q IP Samples ' + textoidl('(\Deltaz=2.0)'), /NODATA
	XYOUTS, xmin+0.05*xr, ymin+0.9125*yr, titles[i], ALIGNMENT=0.0, CHARSIZE=charsz

	OPLOT, Rproj_array+0.5*dRproj, frac_IPSF, LINESTYLE=0, COLOR=SFcolor
	ERRPLOT, Rproj_array+0.5*dRproj, frac_IPSF-errors_IPSF, frac_IPSF+errors_IPSF, COLOR=SFcolor

	OPLOT, Rproj_array+0.5*dRproj, frac_IPQ, LINESTYLE=2, COLOR=Qcolor
	ERRPLOT, Rproj_array+0.5*dRproj, frac_IPQ-errors_IPQ, frac_IPQ+errors_IPQ, COLOR=Qcolor

;	XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'SF IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 1)])), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'Q IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 0)])), ALIGNMENT=1.0

;	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, textoidl('Median M_{*} = ') + decimal(MEDIAN(dataIP[WHERE(dataIP.SFQ EQ 1)].mstar),2) + ' (SF) ' + $
;		decimal(MEDIAN(dataIP[WHERE(dataIP.SFQ EQ 0)].mstar),2) + ' (Q)', ALIGNMENT=1.0, CHARSIZE=charsz
;	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, textoidl('Median z = ') + decimal(MEDIAN(dataIP[WHERE(dataIP.SFQ EQ 1)].zprimus),2) + ' (SF) ' + $
;		decimal(MEDIAN(dataIP[WHERE(dataIP.SFQ EQ 0)].zprimus),2) + ' (Q)', ALIGNMENT=1.0, CHARSIZE=charsz
;	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0, CHARSIZE=charsz
	IF (i EQ 2) THEN BEGIN
;		XYOUTS, xmin+0.15*xr, ymin+0.825*yr, zlabel, ALIGNMENT=0.0, CHARSIZE=charsz
		XYOUTS, xmin+0.05*xr, ymin+0.05*yr, zlabel, ALIGNMENT=0.0, CHARSIZE=charsz
		LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], COLOR=[SFcolor, Qcolor], CHARSIZE=1.25, BOX=0, /BOTTOM, /RIGHT
	ENDIF

;	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0
  ENDFOR

	IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
