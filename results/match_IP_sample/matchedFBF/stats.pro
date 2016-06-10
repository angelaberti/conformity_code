PRO STATS
;	xmin = 0.2
;	xmax = 1.0

;	ymin = 9.0
;	ymax = 11.5

;	bin = 0.05

	dataAll = mrdfits('~/results/zerodSFQ_all_cart.fits', 1)
	dataAll = dataAll[where((dataAll.field EQ 'cosmos    ') AND (dataAll.targ_weight GE 1) AND (dataAll.zprimus GE 0.7))]
	dataAllSF = dataAll[where(dataAll.SFQ EQ 1)]
	dataAllQ  = dataAll[where(dataAll.SFQ EQ 0)]

;	temp = SUBDIV('../matchedIPsample.fits', dataIPQ, dataIPSF)
;	PLOTHIST, dataIPQ.zprimus, bin=bin, color=cGcolor('Red'), xrange=[xmin,xmax];, yrange=[ymin,ymax]
;	PLOTHIST, dataIPSF.zprimus, bin=bin, color=cGcolor('Blue'), /OVERPLOT
;	PLOT, dataIPSF.zprimus, dataIPSF.mstar, psym=3, color=cGcolor('Blue'), xrange=[xmin,xmax], yrange=[ymin,ymax]
;	PLOT, dataIPQ.zprimus, dataIPQ.mstar, psym=3, color=cGcolor('Red')

	psym = 3
	IPsym = 4

	temp = SUBDIV('matchedIPsample_FBF.fits', dataIPQ, dataIPSF)
	PLOT, dataAll.dz, dataAll.yprop, /NODATA
	OPLOT, dataAllQ.dz, dataAllQ.yprop, psym=psym, color=cgcolor('Orange')
	OPLOT, dataAllSF.dz, dataAllSF.yprop, psym=psym, color=cgcolor('Cyan')
	OPLOT, dataIPQ.dz, dataIPQ.yprop, psym=IPsym, color=cgcolor('Red')
	OPLOT, dataIPSF.dz, dataIPSF.yprop, psym=IPsym, color=cgcolor('Blue')
;	PLOTHIST, dataIPSF.zprimus, bin=bin, color=cGcolor('Cyan'), /OVERPLOT
;	PLOTHIST, dataIPQ.zprimus, bin=bin, color=cGcolor('Orange'), /OVERPLOT
;	OPLOT, dataIPSF.zprimus, dataIPSF.mstar, psym=3, color=cGcolor('Cyan')
;	OPLOT, dataIPQ.zprimus, dataIPQ.mstar, psym=3, color=cGcolor('Green')
END

FUNCTION SUBDIV, datafile, dataIPQ, dataIPSF

	data	= mrdfits(datafile, 1)
	dataIP  = data[where((data.IP EQ 1) AND (data.field EQ 'cosmos    ') AND (data.zprimus GE 0.7))]
	dataIPQ = dataIP[where(dataIP.SFQ EQ 0)]
	dataIPSF = dataIP[where(dataIP.SFQ EQ 1)]
	
	RETURN, dataIP
END
