PRO ipplot3D, zerodInput, f, IPindex
; f = field
; IPindex = index of IP galaxy in field
; plots all galaxies in 3D area centered at location of IP galaxy specified by IPindex

data	= zerodInput[where(zerodInput.field eq f)]

dataIP	= data[where(data.IP eq 1)]
IPmain	= dataIP[IPindex] ; first IP galaxy of field

x0 = redh100()*IPmain.dx
y0 = redh100()*IPmain.dy
z0 = redh100()*IPmain.dz

print,"IP coordinates: ", x0, y0, z0

delta_x = 1 ; Mpc
delta_y = 1 ; Mpc
delta_v = (500.0^2 + dvrad(IPmain.zprimus)^2)^0.5	
delta_z = delta_v/(100.0*redh100())

print, "x width (Mpc): " + string(strcompress(2.0*delta_x))
print, "y width (Mpc): " + string(strcompress(2.0*delta_y))
print, "z width (Mpc): " + string(strcompress(2.0*delta_z))
print, "z width (km/s): " + string(strcompress(2.0*delta_v))

dataNIP	= data[where(data.IP eq 0)] ; all NIP galaxies of field

p = plot3d(data.dx, data.dy, data.dz, '*', $
	xrange=[x0-delta_x, x0+delta_x], $
	yrange=[y0-delta_y, y0+delta_y], $
	zrange=[z0-delta_z, z0+delta_z])

END
