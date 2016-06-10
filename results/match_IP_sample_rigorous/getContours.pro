FUNCTION getContours, data, xbinwidth, ybinwidth, xpts, ypts
;	xbinwidth = 0.05
;	ybinwidth = 0.2

;	data = MRDFITS('~/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1)

	xdata = data.zprimus
	ydata = data.mstar

;	xmin = MIN(xdata)
;	xmax = MAX(xdata)
	xmin = 0.2
	xmax = 1.0

;	ymin = MIN(ydata)
;	ymax = MAX(ydata)
;	yrange = ymax-ymin

	ymin = FLOOR(10.*MIN(ydata))/10.
	ymax = CEIL(10.*MAX(ydata))/10.
;	ymin = 9.1
;	ymax = 11.4

	n_xbins = FLOOR( (xmax-xmin)/xbinwidth )
	n_ybins = CEIL( (ymax-ymin)/ybinwidth )

	hist2d = []
	FOR i=0,n_xbins DO BEGIN
		x_low  = xmin + xbinwidth*FLOAT(i)
		x_high = xmin + xbinwidth*FLOAT(i+1)
;		PRINT, 'x_low: ', x_low, ' x_high: ', x_high

		in_x_range = data[WHERE((xdata GE x_low) AND (xdata LE x_high), /NULL)]
;		PRINT, 'in_x_range: ', n_elements(in_x_range)

                x_row = []
                FOR j=0,n_ybins DO BEGIN
			y_low   = ymin + ybinwidth*FLOAT(j)
			y_high  = ymin + ybinwidth*FLOAT(j+1)
;			PRINT, 'y_low: ', y_low, ' y_high: ', y_high

			IF (in_x_range NE !NULL) THEN $
;			in_unit = n_elements(in_x_range[WHERE((ydata GE y_low) AND (ydata LE y_high), /NULL)]) $
			in_unit = n_elements(in_x_range[WHERE((in_x_range.mstar GE y_low) AND (in_x_range.mstar LE y_high), /NULL)]) $
			ELSE in_unit = 0

;			PRINT, 'in_unit: ', in_unit
                        x_row = [x_row, in_unit]
                ENDFOR

		hist2d = [[hist2d], [x_row]]
	ENDFOR

	xpts	= xmin + xbinwidth*(FINDGEN(n_xbins+1))
	ypts	= ymin + ybinwidth*(FINDGEN(n_ybins+1))
	
;	PRINT, 'x bins: ', n_xbins
;	PRINT, 'x pts:  ', n_elements(xpts)
;	PRINT, 'y bins: ', n_ybins
;	PRINT, 'y pts:  ', n_elements(ypts)

;	levels = 0.9*MAX(hist2d)*FINDGEN(n_levels)/(n_levels-1)
;	levels = levels[1:*]
;	PRINT, 'Levels: ', levels

	hist = TRANSPOSE(hist2d)

;	PRINT, hist
	RETURN, hist
END
