; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; latefrac:	calculates late-type fractions
; BSE:		uses bootstrap error method
; IPmatchFBF:	SF and Q IP populations have been selected to have same mass and redshift distribution WITHIN EACH FIELD
; PHI_._:	IP matching procedure uses high mass cutoff threshold based on Moustakas13 stellar mass function results

PRO latefracplot, outputFormat;, zmin, zmax;, dz_coeff, printEvery
	PHIarray = 3.0 + 0.1*findgen(10)

	W = 1.5
	A = findgen(17) * (!PI*2/16.) 
	USERSYM, W*COS(A), W*SIN(A), /FILL

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	Rmax 		= 10.
	dRproj		= 1.
	n_annuli 	= Ceil(Rmax/dRproj)
	printEvery	= 100

	dataAllallz = MRDFITS('~/results/zerodSFQ_all_cart.fits', 1)

	!p.multi=0
	!p.charsize=1.25

	Rproj_array = float(dRproj)*(findgen(n_annuli+1))	
;	PRINT, Rproj_array

	xmin = 0
	xmax = n_annuli*dRproj
	xr = xmax-xmin

	ymin = -0.04
	ymax = 0.06
	yr = ymax-ymin

	colors = ['red', 'orange', 'yellow', 'green', 'cyan', 'blue', 'purple', 'magenta', 'pink', 'gray']
;	colors = ['r', 'o', 'y', 'g', 'c', 'b', 'p', 'm', 'r', 'o', 'y']
;	colors = []

  ; loop through z thirds
  FOR i=0,2 DO BEGIN

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/matchedIPsampleFBF/PHItest_T' + strtrim(string(i+1, format='(i20)'),2), THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

    ; run through all test PHI
    FOR p=0,n_elements(PHIarray)-1 DO BEGIN
	stringPhi = strtrim(string(PHIarray[p], format='(f20.1)'), 2)
        dataIPallz  = MRDFITS('~/results/match_IP_sample_rigorous/PHItest/matchedIPsampleFBF_PHI' + stringPhi + '.fits', 1)

	data = MRDFITS('~/results/match_IP_sample_rigorous/PHItest/latefrac_T' + strtrim(string(i+1, format='(i20)'),2) + '_targ_weight_IPmatchFBF_PHI' + stringPhi + '_dR1Mpc_BSE.fits', 1)

	n_tot_IPSF  = data.n_tot_IPSF
	n_late_IPSF = data.n_late_IPSF
	errors_IPSF = data.errors_IPSF
	frac_IPSF   = n_late_IPSF/n_tot_IPSF

	IF (p EQ 0) THEN PLOT, [1], [1], xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Normalized Conformity Signal', $
	  title='Norm. Conformity Signal (T' + strtrim(string(i+1, format='(i20)'),2) + ')', /NODATA

	n_tot_IPQ  = data.n_tot_IPQ
	n_late_IPQ = data.n_late_IPQ
	errors_IPQ = data.errors_IPQ
	frac_IPQ   = n_late_IPQ/n_tot_IPQ

	OPLOT, [0,Rmax], [0,0], LINESTYLE=2
	OPLOT, Rproj_array+0.5*dRproj, (frac_IPSF-frac_IPQ)/((frac_IPSF+frac_IPQ)/2), COLOR=cgcolor(colors[p]), THICK=3, PSYM=8; LINESTYLE=0
;	OPLOT, Rproj_array+0.5*dRproj, SQRT(errors_IPSF^2+errors_IPQ^2), COLOR=cgcolor(colors[p]), THICK=3, LINESTYLE=p
;	OPLOT, Rproj_array+0.5*dRproj, errors_IPQ, COLOR=cgcolor(colors[p]), PSYM=2, THICK=3;, LINESTYLE=p
;	ERRPLOT, Rproj_array+0.5*dRproj, frac_IPQ-errors_IPQ, frac_IPQ+errors_IPQ

;	XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'SF IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 1)])), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'Q IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 0)])), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, zlabel, ALIGNMENT=1.0

;	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0

    ENDFOR
	LEGEND, [decimal(-1.*PHIarray, 1)], COLOR=cgcolor(colors[indgen(n_elements(PHIarray))]), PSYM=8, THICK=3, BOX=0, /TOP, /RIGHT, CHARSIZE=1.0

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
  ENDFOR
END
