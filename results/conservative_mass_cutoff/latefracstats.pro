PRO latefracstats, dm

;	SET_PLOT, 'ps'
;	DEVICE, file='latefracplot_zrange.ps', /landscape;, encapsulated=1
;	!p.font=1

	zArray  = [0.2, 0.4, 0.6, 0.8, 1.0]
	dvArray = [500, 1000, 1500, 2000]

	OPENW, 1, 'latefracstats.txt', width=1000

	FOR z_index=0,n_elements(zArray)-2 DO BEGIN
		FOR dv_index=0,n_elements(dvArray)-1 DO BEGIN
			
			inputFile = "latefrac_" + strtrim(string(zArray[z_index], format='(f20.1)'),1) + "_" $
						+ strtrim(string(zArray[z_index+1], format='(f20.1)'),1) + "_zerodSFQ_IP_dv" $
						+ strtrim(string(dvArray[dv_index], format='(i20)'),1) + "kms_dm" $
						+ strtrim(string(dm, format='(f20.1)'),1) + ".fits"
			print, inputFile
			data = mrdfits(inputFile, 1)
		
			frac_IPSF = data.n_late_IPSF/data.n_tot_IPSF
			frac_IPQ  = data.n_late_IPQ/data.n_tot_IPQ

			dfrac_IPSF = poissonError(data.n_late_IPSF, data.n_tot_IPSF)
			dfrac_IPQ  = poissonError(data.n_late_IPQ, data.n_tot_IPQ)

		        IF strmatch(inputFile, "_dm0.5") THEN BEGIN
		                dataIPinputFile = strmid(inputFile, strpos(inputFile, 'zerodSFQ_IP'))
		                dataIP  = mrdfits(dataIPinputFile, 1)
		                dm      = 0.5
		        ENDIF ELSE BEGIN
		                dataIPinputFile = strmid(inputFile, strpos(inputFile, 'zerodSFQ_IP'))
		                dataIP  = mrdfits(dataIPinputFile, 1)
		                dm      = 0.0
		        ENDELSE
		
			dataIPcut = dataIP[where((dataIP.zprimus ge zArray[z_index]) and (dataIP.zprimus le zArray[z_index+1]))]
			Ntotal	= strcompress(n_elements(dataIPcut))
			NIP 	= float(n_elements(dataIPcut[where(dataIPcut.IP eq 1)]))
			IPfrac  = string(NIP/float(n_elements(dataIPcut)),format='(f20.3)')
			
			fields  = dataIP[uniq(dataIP.field, sort(dataIP.field))].field

			sigmaRange000_075 = sigmaRange(inputFile,0,2)
;			sigmaRange000_200 = sigmaRange(inputFile,0,7)
;			sigmaRange000_300 = sigmaRange(inputFile,0,11)
;			sigmaRange000_400 = sigmaRange(inputFile,0,15)

			sigmaRange075_200 = sigmaRange(inputFile,3,7)
			sigmaRange075_300 = sigmaRange(inputFile,3,11)
			sigmaRange075_400 = sigmaRange(inputFile,3,15)

;			IF (strmatch(inputFile, "*_dm0.5*") eq 1) THEN (massOffset="-0.5") ELSE (massOffset="")
;			ipfracArray = ipfrac_by_field(dataIPinputField, fields)

			allStats = [strtrim(dvArray[dv_index]), strtrim(zArray[z_index]), strtrim(zArray[z_index+1]), strtrim(Ntotal), strtrim(IPfrac), strtrim(sigmaRange000_075), strtrim(sigmaRange075_200), strtrim(sigmaRange075_300), strtrim(sigmaRange075_400)]
;			allStats = [strtrim(dvArray[dv_index]), strtrim(Ntotal), strtrim(IPfrac), ipfracArray, strtrim(sigmaRange000_075), strtrim(sigmaRange075_200), strtrim(sigmaRange075_300), strtrim(sigmaRange075_400)]
;			print, allStats

			printf, 1, allStats
		ENDFOR
	ENDFOR

	CLOSE, 1

;	DEVICE, /CLOSE
;	SET_PLOT, 'X'

END
