PRO makeIPplots, fields

	SET_PLOT, 'ps'
	DEVICE, file="IPplots_all.ps", /landscape

	inputFile = "zerodSFQ_IP_dv1500kms_dm0.5.fits"

	FOR f=0,n_elements(fields)-1 DO BEGIN
	  FOR SFQ=0,1 DO BEGIN
	    FOR IP=0,1 DO BEGIN
		ipplot, inputFile, fields[f], 'dz', 'dy', SFQ, IP
	    ENDFOR
	  ENDFOR
	ENDFOR

	DEVICE, /CLOSE
	SET_PLOT, 'X'
END
