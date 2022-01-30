def Singular_Isothermal_Ellipsoid(r_Ein, q_lens, x_arcsec, y_arcsec, x_lens_arcsec, y_lens_arcsec, pos_angle, gamma, theta_shear):
	'''
	NAME:
		Singular_Isothermal_Ellipsoid (SIE)
	PURPOSE:
		Compute SPLE model according to Oguri GLAFIC manual
	INPUTS:
		r_Ein:			the Einstein radius for a SIE profile
		q_lens:			elongation of the lens profile (1 - b/a)
		pos_angle:		rotation angle from the x_axis, anti-clockwise
		x/y_lens_arcsec:	position of the lens along the x/y axis, in arcsec
		x/y_arcsec:		array/matrix of pixel positions along the x/y axis, in arcsec
		gamma:			eventual external shear intensity
		theta_shear:		eventual external shera position angle, (FROM THE X_AXIS, ANTICLOCKWISE)
	OUTPUTS:
		alphax_map, alphay_map:	the deflection angle map in arcsec
		k_map:			array/matrix of the convergence
		shear:			array/matrix of the shear
		mu:			array/matrix of the amplifications
	REVISION HISTORY:
		Written		M. Negrello
		Translated	A. Enia		Semptember  2015
	'''
	import numpy as np
	import image_tools as ImagTools
	from astropy import units as u

	if q_lens == 1.0:
		q_lens = 0.9999999
	if theta_shear.unit == 'deg':
		theta_shear = theta_shear.to('rad')
#	mass model normalization
	b_SIE = r_Ein/np.sqrt(q_lens) # (Eq.B27, glafic manual) <--- IN RADIANS
	b_SIE = b_SIE.to('arcsec') # converted in arcsec
	s_SIE = 0.0*u.arcsec/np.sqrt(q_lens) # not 0 if there's a flat core
#	array of distances from the lens
	dist_lens_x_arcsec = (x_arcsec - x_lens_arcsec)
	dist_lens_y_arcsec = (y_arcsec - y_lens_arcsec)
#	apply rotation
	(dist_lens_x_rot_arcsec, dist_lens_y_rot_arcsec) = ImagTools.rotate(dist_lens_x_arcsec, dist_lens_y_arcsec, pos_angle)
#	make calculations (following the GLAFIC manual, see appendix B.5)
	phi = np.sqrt(q_lens**2.0*(s_SIE**2.0 + dist_lens_x_rot_arcsec**2.0) + dist_lens_y_rot_arcsec**2.0)
	alphax_rot_map = (b_SIE*q_lens/np.sqrt(1.0 - q_lens**2.0))*np.arctan((np.sqrt(1.0 - q_lens**2.0)*dist_lens_x_rot_arcsec)/(phi + s_SIE))/u.rad #!!!!!!
	alphay_rot_map = b_SIE*q_lens/np.sqrt(1.0 - q_lens**2.0)*np.arctanh((np.sqrt(1.0 - q_lens**2.0)*dist_lens_y_rot_arcsec)/(phi + q_lens**2.0*s_SIE))/u.rad
#	rotate back and add external shear
	(alphax_rot_map_BACK, alphay_rot_map_BACK) = ImagTools.rotate_back(alphax_rot_map, alphay_rot_map, pos_angle)
#	adding external shear (if there's any)
	if theta_shear.unit == 'deg':
		theta_shear = theta_shear.to('rad')
	gamma_x = -gamma * np.cos(2.0*(theta_shear))  
	gamma_y = -gamma * np.sin(2.0*(theta_shear))
	shear_x_arcsec = gamma_x*dist_lens_x_arcsec + gamma_y*dist_lens_y_arcsec 
	shear_y_arcsec = - gamma_x*dist_lens_y_arcsec + gamma_y*dist_lens_x_arcsec
#	Final result: the deflection angle maps
	alphax_map = alphax_rot_map_BACK + shear_x_arcsec
	alphay_map = alphay_rot_map_BACK + shear_y_arcsec

#	Now, the lens shear, the magnification map and the convergence map
	phi_xx_rot = b_SIE*q_lens/phi*(q_lens**2.0*s_SIE**2.0 + dist_lens_y_rot_arcsec**2.0 + s_SIE*phi)/ \
						((1.0 + q_lens**2.0)*s_SIE**2.0 + 2.0*phi*s_SIE + dist_lens_x_rot_arcsec**2.0 + dist_lens_y_rot_arcsec**2.0)
	phi_yy_rot = b_SIE*q_lens/phi*(s_SIE**2.0 + dist_lens_x_rot_arcsec**2.0 + s_SIE*phi)/ \
						((1.0 + q_lens**2.0)*s_SIE**2.0 + 2.0*phi*s_SIE + dist_lens_x_rot_arcsec**2.0 + dist_lens_y_rot_arcsec**2.0)
	phi_xy_rot = - b_SIE*q_lens/phi*(dist_lens_x_rot_arcsec*dist_lens_y_rot_arcsec)/ \
						((1.0 + q_lens**2.0)*s_SIE**2.0 + 2.0*phi*s_SIE + dist_lens_x_rot_arcsec**2.0 + dist_lens_y_rot_arcsec**2.0)
	(phi_xx, phi_yy, phi_xy) = ImagTools.rotate_derivative_back(phi_xx_rot, phi_yy_rot, phi_xy_rot, pos_angle)
	phi_xx_map = phi_xx + gamma_x
	phi_yy_map = phi_yy - gamma_x
	phi_xy_map = phi_xy + gamma_y
#	convergence, shear & magnification map
	k_map = 0.5*(phi_xx_map + phi_yy_map)
	shear_x = 0.5*(phi_xx_map - phi_yy_map)
	shear_y = phi_xy_map 
	shear = np.sqrt(shear_x**2.0 + shear_y**2.0)
	mu_radial = 1.0/(1.0 - k_map + shear)
	mu_tangential = 1.0/(1.0 - k_map - shear)
	mu = 1.0/((1.0 - phi_xx_map)*(1.0 - phi_yy_map) - phi_xy_map**2.0)

	return alphax_map, alphay_map, k_map, shear, mu

def Single_Power_Law_Ellipsoid(r_Ein, q_lens, slope, x_arcsec, y_arcsec, x_lens_arcsec, y_lens_arcsec, pos_angle, gamma, theta_shear):
	'''
	NAME:
		Singular_Power_Law_Ellipsoid (SPLE)
	PURPOSE:
		Compute SPLE model according to Tessore & Metcalf (2015, A&A, 580, 79)
	INPUTS:	
		r_Ein:				the Einstein radius for a SIE profile
		slope:				slope of the 2D profile (0 < slope < 2); if slope = 1 ==> SIE
		q_lens:				elongation of the lens profile (1 - b/a)
		pos_angle:			rotation angle from the x_axis, anti-clockwise
		x/y_lens_arcsec:	position of the lens along the x/y axis, in arcsec
		x/y_arcsec:			array/matrix of pixel positions along the x/y axis, in arcsec
		gamma:				eventual external shear intensity
		theta_shear:		eventual external shera position angle, (FROM THE X_AXIS, ANTICLOCKWISE)
	OUTPUTS:
		alphax_map, alphay_map:	the deflection angle map in arcsec
		k_map:					array/matrix of the convergence
		shear:					array/matrix of the shear
		mu, mu r, mu t:			array/matrix of the amplifications
	REVISION HISTORY:
		Written		A. Enia		January  2016
	'''
	import numpy as np
	import image_tools as ImagTools
	from astropy import units as u

	if q_lens == 1.0:
		q_lens = 0.9999999
#	mass model normalization (the same as before: in SIE is b_SIE*q = r_Ein*q^(-0.5)*q = r_Ein*q^(0.5)
	b = r_Ein.to('arcsec') * np.sqrt(q_lens)
#	array of distances from the lens
	dist_lens_x_arcsec = (x_arcsec - x_lens_arcsec)
	dist_lens_y_arcsec = (y_arcsec - y_lens_arcsec)
#	Transform into the lens f.o.r..
	(dist_lens_x_rot_arcsec, dist_lens_y_rot_arcsec) = ImagTools.rotate(dist_lens_x_arcsec, dist_lens_y_arcsec, pos_angle)
#	Radial Part
	norm = 2*b/(1+q_lens)
	R = np.sqrt((q_lens**2 * dist_lens_x_rot_arcsec**2 + dist_lens_y_rot_arcsec**2))/(r_Ein.to('arcsec')*np.sqrt(q_lens))
	radial_part = norm*(R)**(1 - slope)
#	Angular Part up to the 20th order
	order = 20
	phi = np.arctan2(dist_lens_y_rot_arcsec.value, q_lens*dist_lens_x_rot_arcsec.value) # The numpy arctan2 function is identical to atan2 in C
	cos = np.cos(2*phi)
	sin = np.sin(2*phi)
	f = (1-q_lens)/(1+q_lens)
	T = 2-slope
	nx = np.shape(x_arcsec)[0]
	ny = np.shape(y_arcsec)[0]
	omega_x = np.zeros((ny, nx, order))
	omega_y = np.zeros((ny, nx, order))
	omega_x[:, :, 0] = np.cos(phi)
	omega_y[:, :, 0] = np.sin(phi)
	for n in range(1, order):
		omega_x[:, :, n] = - f*(2*n-T)/(2*n+T)*(cos*omega_x[:, :, n-1] - sin*omega_y[:, :, n-1])
		omega_y[:, :, n] = - f*(2*n-T)/(2*n+T)*(sin*omega_x[:, :, n-1] + cos*omega_y[:, :, n-1])
	angular_part_x = np.sum(omega_x, axis = 2)
	angular_part_y = np.sum(omega_y, axis = 2)
#	Deflection angle = radial part * angular part
	alphax_rot_map_noshear = radial_part * angular_part_x
	alphay_rot_map_noshear = radial_part * angular_part_y
#	De-rotate
	(alphax_map_noshear, alphay_map_noshear) = ImagTools.rotate_back(alphax_rot_map_noshear, alphay_rot_map_noshear, pos_angle)
#	adding external shear (if there's any)
	if theta_shear.unit == 'deg':
		theta_shear = theta_shear.to('rad')
	gamma_x = -gamma*np.cos(2.0*(theta_shear))  
	gamma_y = -gamma*np.sin(2.0*(theta_shear))
	shear_x_arcsec = gamma_x*dist_lens_x_arcsec + gamma_y*dist_lens_y_arcsec 
	shear_y_arcsec = - gamma_x*dist_lens_y_arcsec + gamma_y*dist_lens_x_arcsec
#	Final result: the deflection angle maps
	alphax_map = alphax_map_noshear + shear_x_arcsec
	alphay_map = alphay_map_noshear + shear_y_arcsec

#	Convergence map
	k_map = ((2.0 - slope)/2.0)*R**(-slope)
#	Now, the lens shear (plus eventual external shear)
	R_forshear = np.sqrt((dist_lens_x_rot_arcsec**2 + dist_lens_y_rot_arcsec**2))
	phi_forshear = np.arctan2(dist_lens_y_rot_arcsec.value, dist_lens_x_rot_arcsec.value)
	cos_forshear = np.cos(phi_forshear)
	sin_forshear = np.sin(phi_forshear)
	cos2_forshear = np.cos(2.0*phi_forshear)
	sin2_forshear = np.sin(2.0*phi_forshear) 
	phi_xx_rot = k_map - k_map*cos2_forshear + (1.0 - slope)/R_forshear * (alphax_rot_map_noshear*cos_forshear - alphay_rot_map_noshear*sin_forshear)
	phi_yy_rot = k_map + k_map*cos2_forshear + (1.0 - slope)/R_forshear * (alphax_rot_map_noshear*cos_forshear - alphay_rot_map_noshear*sin_forshear)
	phi_xy_rot = -k_map*sin2_forshear + (1.0 - slope)/R_forshear * (alphax_rot_map_noshear*sin_forshear + alphay_rot_map_noshear*cos_forshear)
	(phi_xx, phi_yy, phi_xy) = ImagTools.rotate_derivative_back(phi_xx_rot, phi_yy_rot, phi_xy_rot, pos_angle)
	shear_x = 0.5*(phi_xx - phi_yy) + gamma_x
	shear_y = phi_xy + gamma_y
	shear = np.sqrt(shear_x**2.0 + shear_y**2.0)

#	Amplification map
	mu_radial = 1.0/(1.0 - k_map + shear)
	mu_tangential = 1.0/(1.0 - k_map - shear)
	mu = mu_radial * mu_tangential

#	Output
	return alphax_map, alphay_map, k_map, shear, mu
