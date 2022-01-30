# -*- coding: utf-8 -*-

def array_positions(nx, ny, pixel_scale):
	import numpy as np
	index_pixel_x = np.zeros((ny,nx))
	index_pixel_y = np.zeros((ny,nx))
	for iy in range(ny):
		index_pixel_x[iy,:] = np.arange(nx)
	for ix in range(nx):
		index_pixel_y[:,ix] = np.arange(ny)
	x_arcsec = pos_pixel2arcsec(index_pixel_x, pixel_scale)
	y_arcsec = pos_pixel2arcsec(index_pixel_y, pixel_scale)
	return index_pixel_x, index_pixel_y, x_arcsec, y_arcsec

def rotate(x, y, theta):
	'''
	REVISION HISTORY:
		Written		A. Enia		September  2015
	'''
	import numpy as np
	# 2D rotation, given (x,y) and the rotation angle
	try:
		if theta.unit == 'deg': theta = np.radians(theta)
	except: pass
	x_rot = +np.cos(theta)*x + np.sin(theta)*y
	y_rot = -np.sin(theta)*x + np.cos(theta)*y
	return x_rot, y_rot;

def rotate_back(x_rot, y_rot, theta):
	'''
	REVISION HISTORY:
		Written		A. Enia		September  2015
	'''
	import numpy as np
	# 2D de-rotation, given (x_rot, y_rot) and the rotation angle
	try:
		if theta.unit == 'deg': theta = np.radians(theta)
	except: pass
	x = + np.cos(theta)*x_rot - np.sin(theta)*y_rot
	y = + np.sin(theta)*x_rot + np.cos(theta)*y_rot
	return x, y;

def rotate_derivative_back(xx_rot, yy_rot, xy_rot, theta):
	'''
	REVISION HISTORY:
		derivative	M. Negrello	January 2016
	'''
	import numpy as np
	# 2D de-rotation for derivatives, given (xx_rot, yy_rot, xy_rot) and the rotation angle
	try:
		if theta.unit == 'deg': theta = np.radians(theta)
	except: pass
	costheta = np.cos(theta)
	sintheta = np.sin(theta)
	xx = xx_rot*costheta**2.0 - 2.0*xy_rot*sintheta*costheta + yy_rot*sintheta**2.0
	yy = xx_rot*sintheta**2.0 + 2.0*xy_rot*sintheta*costheta + yy_rot*costheta**2.0
	xy = xx_rot*sintheta*costheta + xy_rot*(costheta**2.0 - sintheta**2.0) - yy_rot*sintheta*costheta
	return xx, yy, xy

def pos_arcsec2pixel(arcsec, pixel_scale):
	'''
	NAME:
		pos_arcsec2pixel and pos_pixel2arcsec
	PURPOSE:
		Convert from pixels to arcsecond, given a pixel scale.
	REVISION HISTORY:
		Written			M. Negrello	
		Translated 2 Py		A. Enia		September  2015
	'''
	# From arcsec to pixel
	pixel = arcsec/(pixel_scale) - 0.5
	return pixel;

def pos_pixel2arcsec(pixel, pixel_scale):
	# From pixel to arcsec
	arcsec = (pixel + 0.5)*pixel_scale
	return arcsec;

def nint(x):
	'''
	NAME:
		nint
	PURPOSE:
		Round x to the next integer and change its type to integer
	REVISION HISTORY:
		Written		A. Enia		Semptember  2015
	'''
	import numpy as np
	x_int = (np.rint(x)).astype(int)
	return x_int;

def make_header(image, pixel_scale, ra_ref, dec_ref, observer, telescope):
	'''
	NAME:
		make_header
	PURPOSE:
		Pretty self explanatory
	INPUT:
		image			2D array
		pixel_scale		well...
		ra_ref, dec_ref 	well...^2
	REVISION HISTORY:
		Written		A. Enia		September  2015
	'''
	from astropy.io import fits
	import numpy as np
	# Get the header
	hdu = fits.PrimaryHDU(image)
	hdr = hdu.header
	nx, ny = image.shape[1], image.shape[0]
	hdr['naxis1'] = nx
	hdr['naxis2'] = ny
	# Then start modify it.
	if (np.mod(nx, 2) != 0): hdr['crpix1'] = (nx-1)/2 + 1	# Odd
	else: hdr['crpix1'] = nx/2 + 0.5	# Even	
	if (np.mod(ny, 2) != 0): hdr['crpix2'] = (ny-1)/2 + 1	# Odd
	else: hdr['crpix2'] = ny/2 + 0.5	# Even
	hdr['equinox'] = 2000
	hdr['crval1'] = ra_ref.value
	hdr['crval2'] = dec_ref.value
	hdr['ctype1'] = 'RA---TAN'
	hdr['ctype2'] = 'DEC--TAN'
	hdr['RADECSYS'] = 'FK5'
	hdr['CD1_1'] = -pixel_scale.to('deg').value
	hdr['CD1_2'] = 0.0
	hdr['CD2_1'] = 0.0
	hdr['CD2_2'] = pixel_scale.to('deg').value
	hdr['observer'] = observer
	hdr['telescop'] = telescope
	return hdr

def make_ellipse(x_cen_arcsec, y_cen_arcsec, x_arcsec, y_arcsec, dx, dy, rot_angle, ell, mjaxis_min, mjaxis_max):

	'''
	NAME:
		make_ellipse
	PURPOSE:
		Create the ellipses enclosing the pixel actually taking part of the lens modeling
	INPUT:
		x_cen_arcsec, y_cen_arcsec	Center of the image
		x_arcsec, y_arcsec		Meshgrid arrays with the pixel positions in arcsecons
		dx, dy				Position of the ellipse with respect to the center
		rot_angle			Position angle of the ellipse
		ell				Ellipticity
		mjaxis_min, mjaxis_max		Radius, in arcsecond, of the inner and outer part of the ellipse
	REVISION HISTORY:
		Written		A. Enia		July 2016
	'''
	import numpy as np
	# picking a reference point for the ellipses
	x_cen_arcsec_ellipse = x_cen_arcsec + dx
	y_cen_arcsec_ellipse = y_cen_arcsec + dy
	# array of distances from the reference point
	dist_x_arcsec_ellipse = (x_arcsec - x_cen_arcsec_ellipse)    # arcsec
	dist_y_arcsec_ellipse = (y_arcsec - y_cen_arcsec_ellipse)    # arcsec
	# apply rotation
	(dist_x_rot_arcsec, dist_y_rot_arcsec) = rotate(dist_x_arcsec_ellipse, dist_y_arcsec_ellipse, rot_angle)
	# building the ellipses: ellipticity and major axes, other things follows
	a_min, a_max = mjaxis_min*ell, mjaxis_max*ell 
	b_min, b_max = mjaxis_min/ell, mjaxis_max/ell
	distance_min = dist_x_rot_arcsec**2.0/a_min**2.0 + dist_y_rot_arcsec**2.0/b_min**2.0
	distance_max = dist_x_rot_arcsec**2.0/a_max**2.0 + dist_y_rot_arcsec**2.0/b_max**2.0
	return (distance_min >= 1.0) & (distance_max <= 1.0)


def savefits(array, file_name, pixel_scale, hdr, ra, dec, observer, telescope):
	'''
	NAME:
		savefits
	PURPOSE:
		Save a fits file
	INPUT:
		array					2D array to save as a fits
		file_name				name of the fits file
		pixel_scale				pixel scale (for header)
		ra, dec					ra/dec of the center (for header)
		observer, telescope		observer and telescope (for header)
	REVISION HISTORY:
		Written		A. Enia		October 2017
	'''
	from image_tools import make_header
	from astropy.io import fits
	try: hdu = fits.PrimaryHDU(array.value)
	except: hdu = fits.PrimaryHDU(array)
	try: hdu.header = hdr
	except: hdu.header = make_header(array, pixel_scale, ra, dec, observer, telescope)
	hdu.writeto(file_name, overwrite = 'True')
	return
