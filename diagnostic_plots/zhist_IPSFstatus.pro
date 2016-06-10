PRO zhist_IPSFstatus, outputFormat
	IF (string(outputFormat) eq 'ps') THEN BEGIN
		PS_OPEN, '~/figures/zhist_IP_SFstatus', THICK=5, /ENCAP
		DEVICE, /INCH, XS=8, YS=6
        ENDIF ELSE BEGIN
                SET_PLOT, 'X'
        ENDELSE

	!p.multi=0
	!p.charsize=1.25

;	!x.margin=[7,2.5] ; left, right
;	!y.margin=[3.5,1] ; bottom, top

	binsize = 0.02

	path = '../results/conservative_mass_cutoff/IP_data/'

	data = mrdfits(path+'zerodSFQ_IP_dz2.0_dm0.0.fits', 1)
	
	dataIP  = data[where(data.IP eq 1)]
	dataNIP = data[where(data.IP eq 0)]
	dataIPSF= dataIP[where(dataIP.SFQ eq 1)]
	dataIPQ = dataIP[where(dataIP.SFQ eq 0)]

	ymin=0
	ymax=900

	plothist, dataIP.zprimus,   bin=binsize, linestyle=0, yrange=[ymin,ymax], xtitle="redshift", ytitle='number', title='Redshift distributions of IP samples'
	plothist, dataIPSF.zprimus, bin=binsize, linestyle=0, color=cgColor('Blue'), /overplot
	plothist, dataIPQ.zprimus,  bin=binsize, linestyle=0, color=cgColor('Red'), /overplot

	LEGEND, ["All IP", "Late-type IP", "Early-type IP"], color=[cgColor('Black'),cgColor('Blue'),cgColor('Red')], linestyle=[0,0,0], BOX=0, /RIGHT, /TOP
	LEGEND, ['Conservative mass limit', textoidl('\Deltaz = 2.0')], BOX=0, /LEFT, /TOP
	
;	zmin = 0.2
;	zmax = 1.0
;	dz = binsize

;	ymin = 0
;	ymax = 10000

;	zArray = [zmin]
;	FOR z_index=0,(zmax-zmin)/dz DO BEGIN
;		zArray = [zArray, zmin+float(z_index)*dz]
;	ENDFOR
;	print, zArray		

;	IPfrac_zbinned = []
;	n_IParray = []
;	FOR z_index=1,(zmax-zmin)/dz DO BEGIN
;		data_zbin = data[where( (data.zprimus ge zArray[z_index]) AND $
;					(data.zprimus lt zArray[z_index+1]) )]
;		data_zbinIP  = data_zbin[where(data_zbin.IP eq 1)]
;		data_zbinNIP = data_zbin[where(data_zbin.IP eq 0)]

;		n_IP 	= n_elements(data_zbinIP)
;		n_total = n_elements(data_zbin)
		
;		n_IParray = [n_IParray, float(n_IP)]
;		IPfrac_zbinned = [IPfrac_zbinned, float(n_IP)/float(n_total)]

;		PRINT, zArray[z_index], zArray[z_index+1], n_IP, n_total, float(n_IP)/float(n_total)
;	ENDFOR

;	plot, zArray+dz/2., n_IParray, xtitle="redshift", ytitle=textoidl("N_{IP}"), xrange=[zmin,zmax], yrange=[ymin,ymax]
;	plot, zArray+dz/2, IPfrac_zbinned, xtitle="redshift", ytitle="IP fraction", xrange=[zmin,zmax], yrange=[0,1]

        IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
       	SET_PLOT, 'X'
END
