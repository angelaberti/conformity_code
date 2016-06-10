PRO ipfrac_all

;zlist = [0.4, 0.6, 0.8, 1.0]
;dvlist = [500, 750, 1000, 1250] ; +/- km/s

;FOR i=0,n_elements(zlist)-1 DO BEGIN
;	FOR j=0,n_elements(dvlist)-1 DO BEGIN
;		iptest_v4, "zerodSFQ_all_cart.fits", zlist[i], 0.5, dvlist[j], 1000
;	ENDFOR
;ENDFOR

FOREACH element, findfile("MoustakasMassLimits/zerodSFQ_IP*.fits") DO BEGIN
	print, element
	IF (strmid(element, strpos(element, 'dv')+2,1) ne '1') $
	  THEN dv = strmid(element, strpos(element, 'dv')+2,3) ELSE dv = strmid(element, strpos(element, 'dv')+2, 4)
	print, "deltav: " + dv
	print, "latefrac, " + element + ", 12, 0.25, " + dv + ", 5000"
;	latefrac, element, 12, 0.25, float(dv), 1000
	ipfrac, element
ENDFOREACH

FOREACH element, findfile("MoustakasMassLimits/zerodSFQ_IP*.fits") DO BEGIN
	print, element
	IF (strmid(element, strpos(element, 'dv')+2,1) ne '1') $
	  THEN dv = strmid(element, strpos(element, 'dv')+2,3) ELSE dv = strmid(element, strpos(element, 'dv')+2, 4)
	print, "deltav: " + dv
	print, "latefrac, " + element + ", 12, 0.25, " + dv + ", 5000"
;	latefrac, element, 12, 0.25, float(dv), 1000
	ipfrac, element
ENDFOREACH

;OPENW, 1, "ipfrac.txt", width=600

;FOREACH element, findfile("MoustakasMassLimits-0.5/zerodSFQ_IP*.fits") DO ipfrac, element
;FOREACH element, findfile("MoustakasMassLimits/zerodSFQ_IP*.fits") DO ipfrac, element

;CLOSE, 1

END
