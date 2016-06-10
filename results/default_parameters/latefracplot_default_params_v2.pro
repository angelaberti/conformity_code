PRO latefracplot_default_params_v2, outputFormat
  IF (string(outputFormat) eq 'ps') THEN BEGIN
	PS_OPEN, '../../figures/latefracplot_default_params_v2', THICK=5, /ENCAP
	DEVICE, /INCH, XS=8, YS=6
  ENDIF ELSE BEGIN
	SET_PLOT, 'X'
  ENDELSE

;	!p.font=1

	ymin	= 0.65
	ymax	= 0.825
	yr	= ymax - ymin

	!p.charsize = 1.5
	linethickness = 5
	xyoutscharsize = 1.5
	legendcharsize = 1.5

;	zArray  = [0.2, 0.4, 0.6, 0.8, 1.0]
	zArray  = [0.2, 1.0]
	dz_coeffList = [1.0]

	!p.multi=0
;	!p.multi=[0,2,2]
;	!x.margin=[6,1.5] ; left, right
;	!y.margin=[3.5,1] ; bottom, top	

	FOR z_index=0,n_elements(zArray)-2 DO BEGIN
;	FOR z_index=0,0 DO BEGIN
                FOR coeff_index=0,n_elements(dz_coeffList)-1 DO BEGIN
                        inputFile = "latefrac_data_hist/latefrac_" + strtrim(string(zArray[z_index], format='(f20.1)'),1) + "_" $
                                                + strtrim(string(zArray[z_index+1], format='(f20.1)'),1) + "_zerodSFQ_IP_dz" $
                                                + strtrim(string(dz_coeffList[coeff_index], format='(f20.1)'),1) + "_dm0.5.fits"
			print, "Late fraction data: ", inputFile
			data = mrdfits(inputFile, 1)
		
			xmin	= 0.0
			xmax	= max(data.Rmax)+min(data.Rmax)
			xr 	= xmax - xmin

			Nneighbors = strcompress(string(total(data.n_tot_IPSF + data.n_tot_IPQ), format='(i20)'))

			frac_IPSF = data.n_late_IPSF/data.n_tot_IPSF
			frac_IPQ  = data.n_late_IPQ/data.n_tot_IPQ

			dfrac_IPSF = poissonError(data.n_late_IPSF, data.n_tot_IPSF)
			dfrac_IPQ  = poissonError(data.n_late_IPQ, data.n_tot_IPQ)

	                dataIP  = mrdfits('IP_data/'+strmid(inputFile, strpos(inputFile, 'zerodSFQ_IP')), 1)
		
			dataIP_zcut = dataIP[where((dataIP.zprimus ge zArray[z_index]) and (dataIP.zprimus le zArray[z_index+1]))]
			Ntotal	= strcompress(n_elements(dataIP_zcut))
			NIP 	= float(n_elements(dataIP_zcut[where(dataIP_zcut.IP eq 1)]))
			IPfrac  = string(100.*NIP/float(n_elements(dataIP_zcut)),format='(f20.1)')

			medianMassAll	= median(dataIP_zcut.mstar)
			medianMassIP	= median(dataIP_zcut[where(dataIP_zcut.IP eq 1)].mstar)
			meanMassAll	= mean(dataIP_zcut.mstar)
			meanMassIP	= mean(dataIP_zcut[where(dataIP_zcut.IP eq 1)].mstar)
			
			plot, data.Rmax, frac_IPSF, xrange=[xmin,xmax], yrange=[ymin,ymax], title='Default Sample Parameters', $
			  xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', thick=linethickness
			errplot, data.Rmax, frac_IPSF-dfrac_IPSF, frac_IPSF+dfrac_IPSF

			oplot, data.Rmax, frac_IPQ, linestyle=2, thick=linethickness

			errplot, data.Rmax, frac_IPQ-dfrac_IPQ, frac_IPQ+dfrac_IPQ

			IF z_index eq 0 THEN BEGIN
				LEGEND, ["Late-type IP","Early-type IP"], linestyle=[0,2], box=0, thick=linethickness, /TOP, /RIGHT
			ENDIF

			sigmaRange000_075 = sigmaRange(inputFile,0,2)
;			sigmaRange000_200 = sigmaRange(inputFile,0,7)
;			sigmaRange000_300 = sigmaRange(inputFile,0,11)
;			sigmaRange000_400 = sigmaRange(inputFile,0,15)

			sigmaRange075_200 = sigmaRange(inputFile,3,7)
;			sigmaRange075_300 = sigmaRange(inputFile,3,11)
			sigmaRange075_400 = sigmaRange(inputFile,3,15)
	
			xyouts, xmin+0.05*xr, ymin+0.05*yr, strtrim(string(zArray[z_index],format='(f20.1)'),1) + textoidl(" < z < ") + $
						  	    strtrim(string(zArray[z_index+1],format='(f20.1)'),1), ALIGNMENT=0.0

			xyouts, xmin+0.95*xr, ymin+0.45*yr, textoidl("\sigma_{R < 0.75 Mpc}    ") + strtrim(string(sigmaRange000_075,format='(f20.2)'),1), ALIGNMENT=1.0
			xyouts, xmin+0.95*xr, ymin+0.37*yr, textoidl("\sigma_{0.75 < R < 2 Mpc} ") + strtrim(string(sigmaRange075_200,format='(f20.2)'),1), ALIGNMENT=1.0
			xyouts, xmin+0.95*xr, ymin+0.29*yr, textoidl("\sigma_{0.75 < R < 4 Mpc} ") + strtrim(string(sigmaRange075_400,format='(f20.2)'),1), ALIGNMENT=1.0

;			xyouts, xmin+0.95*xr, ymin+0.21*yr, textoidl("Total galaxies ") + strtrim(Ntotal, 1), ALIGNMENT=1.0
			xyouts, xmin+0.95*xr, ymin+0.21*yr, textoidl("N_{neighbors} ") + strtrim(Nneighbors, 1), ALIGNMENT=1.0
			xyouts, xmin+0.95*xr, ymin+0.13*yr, textoidl("N_{IP} ") + strtrim(string(NIP, format='(i20)'), 1) + " ("+strtrim(IPfrac, 1)+"%)", ALIGNMENT=1.0
			xyouts, xmin+0.95*xr, ymin+0.05*yr, textoidl("Median IP mass ") + strtrim(string(medianMassIP, format='(f20.2)'), 1), ALIGNMENT=1.0
		ENDFOR
	ENDFOR

  IF (string(outputFormat) eq 'ps') THEN BEGIN
	PS_CLOSE
  ENDIF
	SET_PLOT, 'X'
END
