PRO latefracplot_single_mass_zbins, outputFormat
  IF (string(outputFormat) eq 'ps') THEN BEGIN
	PS_OPEN, '../../figures/latefracplot_single_mass_zbins', THICK=5, /ENCAP
	DEVICE, /INCH, XS=8, YS=6
  ENDIF ELSE BEGIN
	SET_PLOT, 'X'
  ENDELSE

;	!p.font=1

	ymin	= 0.35
	ymax	= 1.0
	yr	= ymax - ymin

	erase
	!p.charsize = 2
	legend_size=1
	xyouts_size=1
;	linethickness = 2

;	zArray  = [0.2, 0.4, 0.6, 0.8, 1.0]
	zArray  = [0.6, 0.8, 0.2, 0.4, 1.0]
;	zArray  = [0.2, 1.0]
	dz_coeffList = [1.0]

	!p.multi=4
;	!p.multi=[0,2,2]
;	!x.margin=[6,1.5] ; left, right
;	!y.margin=[3.5,2.5] ; bottom, top	

        nx=2.
        ny=2.

        xm=0.08
        ym=0.08

        dx = (1 - 1.5*xm)/nx
        dy = (1 - 2.0*ym)/ny

        FOR z_index=0,n_elements(zArray)-2 DO BEGIN
                IF (z_index le (nx-1)) THEN BEGIN
                        xtitle = 'Projected Radius (Mpc)'
                        xtickformat = ''
                ENDIF ELSE BEGIN
                        xtitle = ''
                        xtickformat = '(A1)'
                ENDELSE
                IF (z_index mod nx eq 0) THEN BEGIN
                        ytitle = 'Late-type Fraction'
                        ytickformat = ''
                ENDIF ELSE BEGIN
                        ytitle = ''
                        ytickformat = '(A1)'
                ENDELSE

;		IF z_index ge nx*(ny-1) THEN title=textoidl('Log (M_{IP} / M_{\odot}) > 10') ELSE title=''
                IF z_index eq 3 THEN BEGIN
                        XYOUTS, 4.25, 1.025, textoidl('Log (M_{IP} / M_{\odot}) > 10'), ALIGNMENT=0.5
                ENDIF

                pos = [ xm + dx*float(z_index mod nx), ym + dy*float(floor(float(z_index)/nx)), $
                        xm + dx*(1. + float(z_index mod nx)), ym + dy*(1. + float(floor(float(z_index)/nx))) ]
                PRINT, pos
			
                FOR coeff_index=0,n_elements(dz_coeffList)-1 DO BEGIN
                        inputFile = "latefrac_data_hist/latefrac_" + strtrim(string(zArray[z_index], format='(f20.1)'),1) + "_" $
                                                + strtrim(string(zArray[z_index]+0.2, format='(f20.1)'),1) + "_zerodSFQ_IP_dz" $
                                                + strtrim(string(dz_coeffList[coeff_index], format='(f20.1)'),1) + "_singleMass.fits"
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
		
			dataIP_zcut = dataIP[where((dataIP.zprimus ge zArray[z_index]) and (dataIP.zprimus le (zArray[z_index]+0.2)))]
			Ntotal	= strcompress(n_elements(dataIP_zcut))
			NIP 	= float(n_elements(dataIP_zcut[where(dataIP_zcut.IP eq 1)]))
			IPfrac  = string(100.*NIP/float(n_elements(dataIP_zcut)),format='(f20.1)')

			medianMassAll	= median(dataIP_zcut.mstar)
			medianMassIP	= median(dataIP_zcut[where(dataIP_zcut.IP eq 1)].mstar)
			meanMassAll	= mean(dataIP_zcut.mstar)
			meanMassIP	= mean(dataIP_zcut[where(dataIP_zcut.IP eq 1)].mstar)
			
			plot, data.Rmax, frac_IPSF, xrange=[xmin,xmax], yrange=[ymin,ymax], xtickformat=xtickformat, ytickformat=ytickformat, $
			  xtitle=xtitle, ytitle=ytitle, title=title, thick=linethickness, position=pos
			errplot, data.Rmax, frac_IPSF-dfrac_IPSF, frac_IPSF+dfrac_IPSF

			oplot, data.Rmax, frac_IPQ, linestyle=2, thick=linethickness

			errplot, data.Rmax, frac_IPQ-dfrac_IPQ, frac_IPQ+dfrac_IPQ

			IF z_index eq nx*ny-1 THEN BEGIN
				LEGEND, ["Late-type IP","Early-type IP"], linestyle=[0,2], box=0, /TOP, /RIGHT, charsize=legend_size
			ENDIF

			sigmaRange000_075 = sigmaRange(inputFile,0,2)
;			sigmaRange000_200 = sigmaRange(inputFile,0,7)
;			sigmaRange000_300 = sigmaRange(inputFile,0,11)
;			sigmaRange000_400 = sigmaRange(inputFile,0,15)

			sigmaRange075_200 = sigmaRange(inputFile,3,7)
;			sigmaRange075_300 = sigmaRange(inputFile,3,11)
			sigmaRange075_400 = sigmaRange(inputFile,3,15)
	
			xyouts, xmin+0.05*xr, ymin+0.05*yr, strtrim(string(zArray[z_index],format='(f20.1)'),1) + textoidl(" < z < ") + $
						  	    strtrim(string(zArray[z_index]+0.2,format='(f20.1)'),1), ALIGNMENT=0.0, charsize=xyouts_size

			xyouts, xmin+0.95*xr, ymin+0.45*yr, textoidl("\sigma_{R < 0.75 Mpc}    ") + strtrim(string(sigmaRange000_075,format='(f20.2)'),1), ALIGNMENT=1.0, charsize=xyouts_size
			xyouts, xmin+0.95*xr, ymin+0.37*yr, textoidl("\sigma_{0.75 < R < 2 Mpc} ") + strtrim(string(sigmaRange075_200,format='(f20.2)'),1), ALIGNMENT=1.0, charsize=xyouts_size
			xyouts, xmin+0.95*xr, ymin+0.29*yr, textoidl("\sigma_{0.75 < R < 4 Mpc} ") + strtrim(string(sigmaRange075_400,format='(f20.2)'),1), ALIGNMENT=1.0, charsize=xyouts_size

;			xyouts, xmin+0.95*xr, ymin+0.21*yr, textoidl("Total galaxies ") + strtrim(Ntotal, 1), ALIGNMENT=1.0, charsize=xyouts_size
			xyouts, xmin+0.95*xr, ymin+0.21*yr, textoidl("N_{neighbors} ") + strtrim(Nneighbors, 1), ALIGNMENT=1.0, charsize=xyouts_size
			xyouts, xmin+0.95*xr, ymin+0.13*yr, textoidl("N_{IP} ") + strtrim(string(NIP, format='(i20)'), 1) + " ("+strtrim(IPfrac, 1)+"%)", ALIGNMENT=1.0, charsize=xyouts_size
			xyouts, xmin+0.95*xr, ymin+0.05*yr, textoidl("Median IP mass ") + strtrim(string(medianMassIP, format='(f20.2)'), 1), ALIGNMENT=1.0, charsize=xyouts_size
		ENDFOR
	ENDFOR

  IF (string(outputFormat) eq 'ps') THEN BEGIN
	PS_CLOSE
  ENDIF
	SET_PLOT, 'X'
END
