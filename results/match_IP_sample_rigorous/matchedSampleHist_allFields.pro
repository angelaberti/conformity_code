PRO matchedSampleHist_allFields, outputFormat
	!P.FONT=0
 
	targPhi = -3.7
	stringPHI = strtrim(string(-1.*targPhi, format='(f20.1)'), 2)

	dataAll = mrdfits('~/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits',1)
        data = dataAll[where((dataAll.IP eq 1) AND (dataAll.targ_weight GE 1))]
	dataSF = data[where(data.SFQ EQ 1)]
	dataQ  = data[where(data.SFQ EQ 0)]

	IPselect = MRDFITS('~/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1)
	SFselect = IPselect[where((IPselect.SFQ EQ 1) AND (IPselect.targ_weight GE 1) AND (IPselect.mstar GE 9.1))]
	Qselect  = IPselect[where((IPselect.SFQ EQ 0) AND (IPselect.targ_weight GE 1) AND (IPselect.mstar GE 9.1))]

	PRINT, MIN(SFselect.mstar)
	PRINT, MIN(Qselect.mstar)

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/latex/figures/matchedSampleHist_allFields', /ENCAP, THICK=5
	        DEVICE, /INCH, XS=8, YS=6, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	!P.MULTI = 0
	!P.CHARSIZE = 1.25

	xminPlot = 9.1
	xmaxPlot = 12.0
	xrPlot = xmaxPlot-xminPlot

	IF outputFormat eq 'ps' THEN SFselectColor = cgColor('Black') ELSE SFselectColor = cgColor('White')
	Qcolor  = cgColor('Red')
	SFcolor = cgColor('Blue')
	textColor = cgColor('Black')

	mass_binwidth	= 0.2

	xtitle = textoidl('log ( M_{stellar} / M') + sunsymbol() + ' )'
	ytitle = 'Number'

	; SET UP PLOT AREA
	b=mass_binwidth
	PLOTHIST, dataSF.mstar, bin=b, MIN=9.1, xrange=[xminPlot, xmaxPlot], xtitle=xtitle, ytitle=ytitle, COLOR=SFcolor
;	PLOTHIST, dataSF.mstar, bin=b, xtitle=xtitle, ytitle=ytitle, COLOR=SFcolor
	PLOTHIST, dataQ.mstar, bin=b, MIN=9.1, COLOR=Qcolor, /OVERPLOT
	PLOTHIST, SFselect.mstar, bin=b, MIN=9.1, COLOR=SFcolor, LINESTYLE=2, /OVERPLOT
	PLOTHIST, Qselect.mstar, bin=b, MIN=9.1, COLOR=Qcolor, LINESTYLE=2, /OVERPLOT

	dm = 0.2
	marray = 9.1+dm*FINDGEN(21)

	FOREACH m,marray DO BEGIN
		mmin=m
		mmax=m+dm
		PRINT, mmin, mmax, $
			n_elements(SFselect[where((SFselect.mstar GE mmin) AND (SFselect.mstar LT mmax), /NULL)]), $
			n_elements(Qselect[where((Qselect.mstar GE mmin) AND (Qselect.mstar LT mmax), /NULL)])
	ENDFOREACH

;	LEGEND, ['Late-type IPs','Early-type IPs'], LINESTYLE=[0,2], COLOR=[cgColor('blue'),cgColor('red')], BOX=0, /BOTTOM, /RIGHT
	LEGEND, ['All SF', 'All Q', 'Selected SF', 'Selected Q'], LINESTYLE=[0,0,2,2], COLOR=[SFcolor,Qcolor,SFcolor,Qcolor], BOX=0, /TOP, /LEFT, THICK=5, NUMBER=2, PSPACING=3

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

FUNCTION twoDecimal, input
	RETURN, strtrim(string(input, format='(f20.2)'),1)
END
