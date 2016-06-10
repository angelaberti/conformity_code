PRO zhist_IPcompare, outputFormat
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '../figures/zhist_IPcompare', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	!p.multi=0

	!p.charsize = 1.25

	binsize = 0.02
	linethickness = 5
	linestyle = 0

	mydefault	= mrdfits('../results/default_parameters/IP_data/zerodSFQ_IP_dz1.0_dm0.5.fits',1)
	NIPdefault	= n_elements(mydefault[where(mydefault.IP eq 1)])
	mycons		= mrdfits('../results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz1.0_dm0.0.fits',1)
	NIPcons		= n_elements(mycons[where(mycons.IP eq 1)])
	single		= mrdfits('../results/single_mass_cutoff/IP_data/zerodSFQ_IP_dz1.0_singleMass.fits',1)	
	NIPsingle	= n_elements(single[where(single.IP eq 1)])

	PLOTHIST, mydefault[where(mydefault.IP eq 1)].zprimus, bin=binsize, thick=linethickness, linestyle=linestyle, xtitle='Redshift', yrange=[0,1450], title='Redshift Distributions of IP Samples', ytitle='Number', /NODATA
	PLOTHIST, mycons[where(mycons.IP eq 1)].zprimus, bin=binsize, thick=linethickness, linestyle=linestyle, color=cgColor('Red'), /overplot
	PLOTHIST, single[where(single.IP eq 1)].zprimus, bin=binsize, thick=linethickness, linestyle=linestyle, color=cgColor('Black'), /overplot
	PLOTHIST, mydefault[where(mydefault.IP eq 1)].zprimus, bin=binsize, thick=linethickness, linestyle=linestyle, color=cgColor('Blue'), /overplot

	LEGEND, [textoidl('Default IP M_{*} limit (N_{IP}: ') + bigint(NIPdefault) + ')', $
		 textoidl('Conservative IP M_{*} limit (N_{IP}: ') + bigint(NIPcons) + ')', $
		 textoidl('Log (M_{IP}/M_{\odot}) > 10 (N_{IP}: ') + bigint(NIPsingle) + ')'], $
		linestyle=[0,0,0], thick=linethickness, color=[cgColor('Blue'),cgColor('Red'),cgColor('Black')], $ ; cgColor('Green'),cgColor('Purple')], $
		box=0, /BOTTOM, /RIGHT, charsize=1.1

        IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
       	SET_PLOT, 'X'
END
