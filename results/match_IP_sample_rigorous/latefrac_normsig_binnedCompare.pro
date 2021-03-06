; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; normsigplot:	plots normalized conformity signal
; JKE:		uses jackknife error method
; IPmatchFBF:	SF and Q IP populations have been selected to have same mass and redshift distribution WITHIN EACH FIELD
; PHI37:	IP matching procedure uses high mass cutoff threshold based on Moustakas13 stellar mass function results

PRO latefrac_normsig_binnedCompare, n_bins, COSMOScomp, outputFormat
;	n_bins=2

	!P.FONT=0

	PHI = 3.7
	stringPHI = strtrim(string(PHI, format='(f20.1)'),2)

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	; read in data files
;	IPdataAllz  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1)
	IPdataAllz  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI' + stringPHI + '.fits', 1)
	dataAllallz = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)

	Rmax 		= 5.
	dRproj		= 1.
	n_annuli 	= CEIL(Rmax/dRproj)
	plotMax 	= Rmax-1

        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

; eliminate data with targ_weight < 1
        dataAllallz = dataAllallz[where(dataAllallz.targ_weight GE 1.)]
        IPdataAllz  = IPdataAllz[where(IPdataAllz.targ_weight GE 1.)]

	!P.CHARSIZE=2
	charsz = 0.9

	Rproj_array = FINDGEN(Rmax+1)
	Rplot_array = FINDGEN(Rmax) + 0.5*dRproj

	xmin = 0
	xmax = 5
	xr = xmax-xmin

	ymin = -0.02
	ymax = 0.09
	yr = ymax-ymin

	IF (string(outputFormat) EQ 'ps') THEN BEGIN
;		PS_OPEN, '~/latex/figures/normsigplot_JKE_IPmatchFBF_PHI37_panels', THICK=5, /ENCAP
;		PS_OPEN, '~/latex/figures/latefrac_normsig_binnedCompare', THICK=5, /ENCAP
		PS_OPEN, '~/latex/figures/latefrac_normsig_binnedCompare_noCosmos', THICK=5, /ENCAP
		DEVICE, /INCH, XS=8, YS=8, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
		SET_PLOT, 'X'
		THICK=1
	ENDELSE

;	colors	= ['purple', 'magenta', 'orange']
	colors	= ['PUR8','PUR6','PUR4']
	colorsSF= ['BLU8','BLU6','BLU4']
	colorsQ	= ['RED8','RED6','RED4']
	COSMOS_ls = 2
	ls	= [0,3,4]
	lsSF	= [0,3,4]
	lsQ	= [2,5,6]

	deltaR	= 2*[-0.05, 0, 0.05]

	IF (n_bins EQ 2) THEN BEGIN
		zArray		= [0.2, MEDIAN(IPdataAllz.zprimus), 1.0]
		zSuffixArray	= ['H1', 'H2']
		zBinArray	= zSuffixArray

		massDir		= 'massHalves'
		massArray	= [MIN(IPdataAllz.mstar), MEDIAN(IPdataAllz.mstar), MAX(IPdataAllz.mstar)]
		massSuffixArray = ['H1', 'H2']
		massBinArray	= ['M1equalIP', 'M2equalIP']

		colors		= colors[0:n_elements(colors)-2]
		colorsSF	= colorsSF[0:n_elements(colorsSF)-2]
		colorsQ		= colorsQ[0:n_elements(colorsQ)-2]
		ls		= ls[0:n_elements(ls)-2]
		lsSF		= lsSF[0:n_elements(lsSF)-2]
		lsQ		= lsQ[0:n_elements(lsQ)-2]
;		deltaR		= deltaR[1:n_elements(deltaR)-1]

	ENDIF ELSE BEGIN ; n_bins=3
		lower_index	= CEIL(n_elements(IPdataAllz)/3.)
		upper_index	= 2*FLOOR(n_elements(IPdataAllz)/3.)

		orderedRedshifts= IPdataAllz[SORT(IPdataAllz.zprimus)].zprimus
		zArray		= [0.2, orderedRedshifts[lower_index], orderedRedshifts[upper_index], 1.0]
		zSuffixArray	= ['T1', 'T2', 'T3']
		zBinArray	= zSuffixArray

		massDir		= 'massThirds' 
		orderedMasses	= IPdataAllz[SORT(IPdataAllz.mstar)].mstar
		massArray	= [MIN(IPdataAllz.mstar), orderedMasses[lower_index], orderedMasses[upper_index], MAX(IPdataAllz.mstar)]
		massSuffixArray	= ['T1', 'T2', 'T3']
		massBinArray	= ['M1equalIP', 'M2equalIP', 'M3equalIP']
	ENDELSE

	PRINT, 'N bins: ', n_bins
	PRINT, 'z array: ', zArray
	PRINT, 'Mass array: ', massArray

	NX = 2
	NY = 2
	ERASE
	!P.MULTI=NX*NY
  	posArray = grid_array_deluxe(NX,NY,0.0,0.0)

; ___________________
; |   z    |  mass  |
; |________|________|
; |latefrac|latefrac|
; |___2____|___3____|
; |normsig |normsig |
; |___0____|___1____|

;;====================================================================================================
;;====================================================================================================
;; (0) LOWER LEFT - NORMSIG Z
;;====================================================================================================
;;====================================================================================================
	legend = []
  FOR n=0,n_elements(zArray)-2 DO BEGIN
	z_low   = zArray[n]
	z_high  = zArray[n+1]

    IF (n EQ 0) THEN ( zlabel = 'z = [' + decimal(z_low,2) + ', ' + decimal(z_high,2) + ']' ) $
	ELSE ( zlabel = '[' + decimal(z_low,2) + ', ' + decimal(z_high,2) + ']' )
	legend = [legend, zlabel]

  ; INITIALIZE PLOT AREA
    IF (n EQ 0) THEN BEGIN
	PLOT, [xmin,xmax], [0,0], xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle=textoidl('R_{proj} (Mpc)'), ytitle=textoidl('\xi_{norm}'), $
	  LINESTYLE=1, THICK=2, POSITION=posArray[*,0]
;	XYOUTS, xmin+0.50*xr, ymin+0.90*yr, 'Redshift Bins', ALIGNMENT=0.5, CHARSIZE=charsz
    ENDIF

    IF (COSMOScomp EQ 1) AND (n_bins EQ 2) THEN BEGIN
	COSMOSbin = 'H2'
	dataNC = MRDFITS('~/conformity/results/match_IP_sample_rigorous/jackknife_error/noCosmos/zBins/normsig_'+COSMOSbin+'_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_JKE_noCosmos.fits', 1)
	normsigNC	= dataNC.normsig
	normsigNC_errors= dataNC.normsig_errors
	OPLOT, Rplot_array[0:plotMax], normsigNC[0:plotMax], COLOR=cgColor('BLK5'), LINESTYLE=COSMOS_ls
	ERRPLOT, Rplot_array+deltaR[2], normsigNC-normsigNC_errors, normsigNC+normsigNC_errors, COLOR=cgColor('BLK5')
    ENDIF

	data = MRDFITS('~/conformity/results/match_IP_sample_rigorous/jackknife_error/zBins/normsig_' + zSuffixArray[n] + '_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_JKE.fits', 1)
	normsig		= data.normsig
	normsig_errors	= data.normsig_errors
	OPLOT, Rplot_array[0:plotMax], normsig[0:plotMax], COLOR=cgcolor(colors[n]), LINESTYLE=ls[n]
	ERRPLOT, Rplot_array+deltaR[n], normsig-normsig_errors, normsig+normsig_errors, COLOR=cgcolor(colors[n])
  ENDFOR
  IF (COSMOScomp EQ 1) THEN BEGIN
	LEGEND, [legend, 'w/o COSMOS, '+legend[n_bins-1]], LINESTYLE=[ls,COSMOS_ls], COLOR=[cgColor(colors),cgColor('BLK5')], BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3, CHARSIZE=charsz
  ENDIF ELSE BEGIN
	LEGEND, legend, LINESTYLE=ls, COLOR=[cgColor(colors)], BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3, CHARSIZE=charsz
  ENDELSE

;;====================================================================================================
;;====================================================================================================
;; (1) LOWER RIGHT - NORMSIG MASS
;;====================================================================================================
;;====================================================================================================
	legend	= []
  FOR n=0,n_elements(massArray)-2 DO BEGIN
	mass_low   = massArray[n]
	mass_high  = massArray[n+1]

	IF (n EQ 0) THEN ( masslabel = textoidl('log ( M_{*}/M') + sunsymbol() + ' ) = [' + decimal(mass_low,2) + ', ' + decimal(mass_high,2) + ']' ) $
	  ELSE ( masslabel = '[' + decimal(mass_low,2) + ', ' + decimal(mass_high,2) + ']' )
	legend = [legend, masslabel]

  ; INITIALIZE PLOT AREA
    IF (n EQ 0) THEN BEGIN
	PLOT, [xmin,xmax], [0,0], xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle=textoidl('R_{proj} (Mpc)'), ytickformat='(A1)', LINESTYLE=1, THICK=2, POSITION=posArray[*,1]
;	XYOUTS, xmin+0.50*xr, ymin+0.90*yr, 'Stellar Mass Bins', ALIGNMENT=0.5, CHARSIZE=charsz
    ENDIF

    IF (COSMOScomp EQ 1) AND (n_bins EQ 2) THEN BEGIN
	COSMOSbin = 'H2'
	dataNC = MRDFITS('~/conformity/results/match_IP_sample_rigorous/jackknife_error/noCosmos/massBins/normsig_'+COSMOSbin+'_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_JKE_noCosmos.fits', 1)
	normsigNC	= dataNC.normsig
	normsigNC_errors= dataNC.normsig_errors
	OPLOT, Rplot_array[0:plotMax], normsigNC[0:plotMax], COLOR=cgColor('BLK5'), LINESTYLE=COSMOS_ls
	ERRPLOT, Rplot_array+deltaR[2], normsigNC-normsigNC_errors, normsigNC+normsigNC_errors, COLOR=cgColor('BLK5')
    ENDIF

	data = MRDFITS('~/conformity/results/match_IP_sample_rigorous/jackknife_error/massBins/normsig_' + massSuffixArray[n] + '_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_JKE.fits', 1)
	normsig		= data.normsig
	normsig_errors	= data.normsig_errors
	OPLOT, Rplot_array[0:plotMax], normsig[0:plotMax], COLOR=cgcolor(colors[n]), LINESTYLE=ls[n]
	ERRPLOT, Rplot_array+deltaR[n], normsig-normsig_errors, normsig+normsig_errors, COLOR=cgcolor(colors[n])
  ENDFOR

  IF (COSMOScomp EQ 1) THEN BEGIN
	LEGEND, [legend, 'w/o COSMOS, '+legend[n_bins-1]], LINESTYLE=[ls,COSMOS_ls], COLOR=[cgColor(colors),cgColor('BLK5')], BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3, CHARSIZE=charsz
  ENDIF ELSE BEGIN
	LEGEND, legend, LINESTYLE=ls, COLOR=cgColor(colors), BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3, CHARSIZE=charsz
  ENDELSE

;;====================================================================================================
;;====================================================================================================

	ymin = 0.71
	ymax = 0.88
	yr = ymax-ymin

;;====================================================================================================
;;====================================================================================================
;; (2) UPPER LEFT - LATEFRAC Z
;;====================================================================================================
;;====================================================================================================
	legend = []
  FOR n=0,n_elements(zArray)-2 DO BEGIN
	z_low   = zArray[n]
	z_high  = zArray[n+1]

    IF (n EQ 0) THEN ( zlabel = 'z = [' + decimal(z_low,2) + ', ' + decimal(z_high,2) + ']' ) $
	ELSE ( zlabel = '[' + decimal(z_low,2) + ', ' + decimal(z_high,2) + ']' )
	legend = [legend, zlabel, 'SF IP', 'Q IP']

  ; INITIALIZE PLOT AREA
    IF (n EQ 0) THEN BEGIN
	PLOT, [xmin,xmax], [0,0], xrange=[xmin,xmax], yrange=[ymin,ymax], xtickformat='(A1)', ytitle=textoidl('f_{late}'), LINESTYLE=1, THICK=2, POSITION=posArray[*,2]
	XYOUTS, xmin+0.05*xr, ymin+0.925*yr, 'Redshift Bins', ALIGNMENT=0.0, CHARSIZE=charsz
    ENDIF

;    IF (COSMOScomp EQ 1) THEN BEGIN
;	dataNC = MRDFITS('~/conformity/results/match_IP_sample_rigorous/jackknife_error/noCosmos/zBins/normsig_H2_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_BSE_noCosmos.fits', 1)
;	normsigNC	= dataNC.normsig
;	normsigNC_errors= dataNC.normsig_errors
;	OPLOT, Rplot_array[0:plotMax], normsigNC[0:plotMax], COLOR=cgColor('BLK5'), LINESTYLE=COSMOS_ls
;	ERRPLOT, Rplot_array+deltaR[1], normsigNC-normsigNC_errors, normsigNC+normsigNC_errors, COLOR=cgColor('BLK5')
;    ENDIF

	data = MRDFITS('~/conformity/results/match_IP_sample_rigorous/latefrac_' + zBinArray[n] + '_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_BSE.fits', 1)
        n_tot_IPSF	= data.n_tot_IPSF
        n_late_IPSF	= data.n_late_IPSF
        errors_IPSF	= data.errors_IPSF
        n_tot_IPQ 	= data.n_tot_IPQ
        n_late_IPQ	= data.n_late_IPQ
        errors_IPQ	= data.errors_IPQ

        frac_IPSF	= n_late_IPSF/n_tot_IPSF
        frac_IPQ	= n_late_IPQ/n_tot_IPQ

        OPLOT, Rplot_array, frac_IPSF[0:plotMax], LINESTYLE=lsSF[n], color=cgColor(colorsSF[n])
        ERRPLOT, Rplot_array+0.5*deltaR[0], frac_IPSF-errors_IPSF, frac_IPSF+errors_IPSF, color=cgColor(colorsSF[n])

        OPLOT, Rplot_array, frac_IPQ[0:plotMax], LINESTYLE=lsQ[n], color=cgColor(colorsQ[n])
        ERRPLOT, Rplot_array+0.5*deltaR[2], frac_IPQ-errors_IPQ, frac_IPQ+errors_IPQ, color=cgColor(colorsQ[n])
  ENDFOR

;  IF (COSMOScomp EQ 1) THEN BEGIN
;	LEGEND, ['', legend, 'w/o COSMOS, '+legend[n_bins-1]], LINESTYLE=[1,lsSF,lsQ,COSMOS_ls], $
;	  COLOR=[cgColor('white'),cgColor(colorsSF),cgColor(colorsQ),cgColor('BLK5')], $
;	  BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3, CHARSIZE=charsz
;  ENDIF ELSE BEGIN
	LINESTYLEarray = [1, 1, lsSF[0], lsQ[0], 1, lsSF[1], lsQ[1]]
	COLORarray = [cgColor('white'), cgColor('white'), cgColor(colorsSF[0]), cgColor(colorsQ[0]), $
					cgColor('white'), cgColor(colorsSF[1]), cgColor(colorsQ[1])]
	IF (n_bins EQ 3) THEN BEGIN
		LINESTYLEarray = [LINESTYLEarray, 1, lsSF[2], lsQ[2]]
		COLORarray = [COLORarray, cgColor('white'), cgColor(colorsSF[2]), cgColor(colorsQ[2])]
	ENDIF
	LEGEND, ['', legend], LINESTYLE=LINESTYLEarray, $
	  COLOR=COLORarray, $
	  BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3, CHARSIZE=charsz
;  ENDELSE

;  IF (COSMOScomp EQ 1) THEN BEGIN
;	LEGEND, ['', legend, 'w/o COSMOS, '+legend[n_bins-1]], LINESTYLE=[1,ls,COSMOS_ls], COLOR=[cgColor('white'),cgColor(colors),cgColor('BLK5')], BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3, CHARSIZE=charsz
;  ENDIF ELSE BEGIN
;	LEGEND, ['', legend], LINESTYLE=[1, ls], COLOR=[cgColor('white'),cgColor(colors)], BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3, CHARSIZE=charsz
;  ENDELSE

;;====================================================================================================
;;====================================================================================================
;; (3) UPPER RIGHT - LATEFRAC MASS
;;====================================================================================================
;;====================================================================================================
	legend	= []
  FOR n=0,n_elements(massArray)-2 DO BEGIN
	mass_low   = massArray[n]
	mass_high  = massArray[n+1]

	IF (n EQ 0) THEN ( masslabel = textoidl('log ( M_{*}/M') + sunsymbol() + ' ) = [' + decimal(mass_low,2) + ', ' + decimal(mass_high,2) + ']' ) $
	  ELSE ( masslabel = '[' + decimal(mass_low,2) + ', ' + decimal(mass_high,2) + ']' )
	legend = [legend, masslabel, 'SF IP', 'Q IP']

  ; INITIALIZE PLOT AREA
    IF (n EQ 0) THEN BEGIN
	PLOT, [xmin,xmax], [0,0], xrange=[xmin,xmax], yrange=[ymin,ymax], xtickformat='(A1)', ytickformat='(A1)', LINESTYLE=1, THICK=2, POSITION=posArray[*,3]
	XYOUTS, xmin+0.05*xr, ymin+0.925*yr, 'Stellar Mass Bins', ALIGNMENT=0.0, CHARSIZE=charsz
    ENDIF

;    IF (COSMOScomp EQ 1) THEN BEGIN
;	dataNC = MRDFITS('~/conformity/results/match_IP_sample_rigorous/jackknife_error/noCosmos/massBins/normsig_H2_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_BSE_noCosmos.fits', 1)
;	normsigNC	= dataNC.normsig
;	normsigNC_errors= dataNC.normsig_errors
;	OPLOT, Rplot_array[0:plotMax], normsigNC[0:plotMax], COLOR=cgColor('BLK5'), LINESTYLE=COSMOS_ls
;	ERRPLOT, Rplot_array+deltaR[1], normsigNC-normsigNC_errors, normsigNC+normsigNC_errors, COLOR=cgColor('BLK5')
;    ENDIF

	data = MRDFITS('~/conformity/results/match_IP_sample_rigorous/'+massDir+'/latefrac_' + massBinArray[n] + '_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_BSE.fits', 1)
        n_tot_IPSF	= data.n_tot_IPSF
        n_late_IPSF	= data.n_late_IPSF
        errors_IPSF	= data.errors_IPSF
        n_tot_IPQ 	= data.n_tot_IPQ
        n_late_IPQ	= data.n_late_IPQ
        errors_IPQ	= data.errors_IPQ

        frac_IPSF	= n_late_IPSF/n_tot_IPSF
        frac_IPQ	= n_late_IPQ/n_tot_IPQ

        OPLOT, Rplot_array, frac_IPSF[0:plotMax], LINESTYLE=lsSF[n], color=cgColor(colorsSF[n])
        ERRPLOT, Rplot_array+0.5*deltaR[0], frac_IPSF-errors_IPSF, frac_IPSF+errors_IPSF, color=cgColor(colorsSF[n])

        OPLOT, Rplot_array, frac_IPQ[0:plotMax], LINESTYLE=lsQ[n], color=cgColor(colorsQ[n])
        ERRPLOT, Rplot_array+0.5*deltaR[2], frac_IPQ-errors_IPQ, frac_IPQ+errors_IPQ, color=cgColor(colorsQ[n])
  ENDFOR

;  IF (COSMOScomp EQ 1) THEN BEGIN
;	LEGEND, ['', legend, 'w/o COSMOS, '+legend[n_bins-1]], LINESTYLE=[1,lsSF,lsQ,COSMOS_ls], $
;	  COLOR=[cgColor('white'),cgColor(colorsSF),cgColor(colorsQ),cgColor('BLK5')], $
;	  BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3, CHARSIZE=charsz
;  ENDIF ELSE BEGIN
	LINESTYLEarray = [1, 1, lsSF[0], lsQ[0], 1, lsSF[1], lsQ[1]]
	COLORarray = [cgColor('white'), cgColor('white'), cgColor(colorsSF[0]), cgColor(colorsQ[0]), $
					cgColor('white'), cgColor(colorsSF[1]), cgColor(colorsQ[1])]
	IF (n_bins EQ 3) THEN BEGIN
		LINESTYLEarray = [LINESTYLEarray, 1, lsSF[2], lsQ[2]]
		COLORarray = [COLORarray, cgColor('white'), cgColor(colorsSF[2]), cgColor(colorsQ[2])]
	ENDIF
	LEGEND, ['', legend], LINESTYLE=LINESTYLEarray, $
	  COLOR=COLORarray, $
	  BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3, CHARSIZE=charsz
;  ENDELSE

;  IF (COSMOScomp EQ 1) THEN BEGIN
;	LEGEND, ['', legend, 'w/o COSMOS, '+legend[n_bins-1]], LINESTYLE=[1,ls,COSMOS_ls], COLOR=[cgColor('white'),cgColor(colors),cgColor('BLK5')], BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3, CHARSIZE=charsz
;  ENDIF ELSE BEGIN
;	LEGEND, legend;, LINESTYLE=ls, COLOR=cgColor(colors), BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3, CHARSIZE=charsz
;	LEGEND, ['',legend], LINESTYLE=[1,ls], COLOR=[cgColor('white'),cgColor(colors)], BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3, CHARSIZE=charsz
;  ENDELSE

	IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
