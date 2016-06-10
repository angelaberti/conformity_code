PRO latefrac_vs_mstar, outputFormat
        IF (string(outputFormat) EQ 'ps') THEN BEGIN
		PS_OPEN, '~/figures/latefrac_vs_mstar', THICK=5, /ENCAP
		DEVICE, /INCH, XS=8, YS=6
        ENDIF ELSE BEGIN
                SET_PLOT, 'X'
        ENDELSE

	binsize = 0.1

	A = FINDGEN(17) * (!PI*2/16.)
	USERSYM, 0.5*COS(A), 0.5*SIN(A), /FILL
	ptsymbol = 8

	ERASE
	!p.multi=1
	!p.charsize=1.5

	zmin = 0.2
	zmax = 1.0

	ymin = 0.2
	ymax = 1.05
	
	dataAll = mrdfits('~/conformity/results/zerodSFQ_all_cart.fits', 1)
	dataAllcomp = mrdfits('~/conformity/results/conservative_mass_cutoff/allAboveMassCompLim-0.0.fits', 1)
	dataIP	= mrdfits('~/conformity/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits', 1)
	dataIP = dataIP[WHERE(dataIP.IP EQ 1)]

	mRange = MAX(dataAll.mstar)-MIN(dataAll.mstar)
	dm = 0.3

	mArray = FLOOR(MIN(dataAll.mstar)) + dm*INDGEN(CEIL(mRange/dm)+1)
	PRINT, mArray

	colors = ['Magenta', 'Red', 'Blue', 'green', 'purple']

	fracSF	= []
;	massBin	= []

	FOR m=1,n_elements(mArray)-2 DO BEGIN
		mMin = mArray[m]
		mMax = mArray[m+1]

		data = dataAll[WHERE((dataAll.mstar GE mMin) AND (dataAll.mstar LT mMax), /NULL)]
;		zbin = mean([zmin+z_index*binsize, zmin+(z_index+1)*binsize])

		fracSF_mbin = TOTAL(data[where(data.SFQ EQ 1, /NULL)].targ_weight)/TOTAL(data.targ_weight)

		fracSF = [fracSF, fracSF_mbin]
	ENDFOR

	PLOT, mArray[0:-1]+0.5*dm, fracSF, xrange=MINMAX(mArray), yrange=[ymin,ymax], LINESTYLE=0, $
	xtitle=textoidl('log (M_{stellar}/M_{\odot})'), ytitle='Late-type Fraction'

        IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
       	SET_PLOT, 'X'
END
