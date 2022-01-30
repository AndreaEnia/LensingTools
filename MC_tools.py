# -*- coding: utf-8 -*-

def bayesian_evidence_uvplane(lg_regterm):
	"""
	NAME:
		bayesian_evidence_uvplane
	PURPOSE:
		The Bayesian approach to the RSI method suggests the use of this quantity, the evidence, to find the best regularization term given a certain lens
		model, and to rank the models.
	INPUTS:
		reg_term: 		regularization term
		H_matrix:		regularization matrix
		noise_map:		the noise map (pixel per pixel or a estimated mean level)
		F_matrix:		as in eq. 3 of Dye et al., 2008, MNRAS
		f_real, f_imag:	as in eq. # of Bussman et al., 2013,
		D_matrix:		as in eq. 3 of Dye et al., 2008, MNRAS
		term3_part:		as in eq. 5 of Dye et al., 2008, MNRAS
		term5:			as in eq. 5 of Dye et al., 2008. MNRAS
	OUTPUTS:
		evidence:		as in Eq. 5 of Dye et al., 2008, MNRAS
	"""
	import numpy as np

	reg_term = 10.0**(lg_regterm)
	if np.shape(reg_term) == (1,):
		reg_term = reg_term[0]
	dimens = np.shape(H_matrix)[0]
	aux = F_matrix + reg_term*H_matrix
	S_matrix = np.linalg.solve(aux, D_matrix) # Never ever invert a matrix
	realFT_image_plane_matrix, imagFT_image_plane_matrix = np.zeros(n_vis), np.zeros(n_vis)
	realFT_image_plane_matrix, imagFT_image_plane_matrix = np.sum(f_real_aux*S_matrix, axis = 1), np.sum(f_imag_aux*S_matrix, axis = 1)
	chi_sq_min = np.sum(wgt*(realFT_image_plane_matrix - real_tab)**2.0) + np.sum(wgt*(imagFT_image_plane_matrix - imag_tab)**2.0)
	term2 = np.linalg.slogdet(aux)[1]
	term3 = term3_part + dimens*np.log(reg_term)
	term4_part = np.dot(S_matrix, np.dot(H_matrix, np.transpose(S_matrix)))
	term4 = reg_term*term4_part
	evidence = .5*(chi_sq_min + term2 - term3 + term4 + term5)
	# Necessary to stop SLSQP/CCOBYLA when they fuck up the evidence, giving np.inf as a result.
	if evidence == np.inf:
		return 67101551;
	return evidence

def bayesian_evidence(lg_regterm):
	"""
	NAME:
		bayesian_evidence
	PURPOSE:
		From the Bayesian approach to the RSI by Suyu et. al (2006), this quantity is fundamental to search in the parameters space for the best-fit
		regularization term and the best-fit lens model parameters
	INPUTS:
		lamb: 		regularization term
		H_matrix:	regularization matrix
		noise_map:	the noise map (pixel per pixel or a estimated mean level)
		F_matrix:	as in eq. 3 of Dye et al., 2008, MNRAS
		D_matrix:	as in eq. 3 of Dye et al., 2008, MNRAS
		term3_part:	as in eq. 5 of Dye et al., 2008, MNRAS
		term5:		as in eq. 5 of Dye et al., 2008. MNRAS
	OUTPUTS:
		evidence:	as in Eq. 5 of Dye et al., 2008, MNRAS
	"""
	import numpy as np
	import subprocess

	reg_term = 10.0**(lg_regterm)
	if np.shape(reg_term) == (1,):
		reg_term = reg_term[0]
	dimens = np.shape(H_matrix)[0]
	nx, ny = np.shape(image_plane)[1], np.shape(image_plane)[0]
	aux = F_matrix + reg_term*H_matrix
	S_matrix = np.linalg.solve(aux, D_matrix)
	Image_array = np.zeros(ndim_image)
	Image_array = np.sum(f_aux*S_matrix, axis = 1)
	image_plane_output = np.zeros((ny,nx))
	image_plane_output[iy_image_cut[:], ix_image_cut[:]] = Image_array[:]
	chi_sq_min = np.sum((image_plane_output[GOOD_image] - image_plane[GOOD_image])**2.0/noise_map[GOOD_image]**2.0)
	term2 = np.linalg.slogdet(aux)[1]
	term3 = term3_part + dimens*np.log(reg_term)
	term4_part = np.dot(S_matrix, np.dot(H_matrix, np.transpose(S_matrix)))
	term4 = reg_term*term4_part
	evidence = .5*(chi_sq_min + term2 - term3 + term4 + term5)
	# Necessary to stop SLSQP/CCOBYLA when they fuck up the evidence, giving np.inf as a result.
	if evidence == np.inf:
		return 67101551
			
	return evidence





def ln_epsilon(array_param, nx, ny, x_arcsec, y_arcsec, pixel_scale, x_cen_arcsec, y_cen_arcsec, ny_SP, nx_SP, pixel_scale_SP, tass_number, \
		area_image, db_image_plane, db_noise_map, db_psf_map):
	global image_plane
	global noise_map
	global psf_map
	image_plane, noise_map, psf_map = db_image_plane, db_noise_map, db_psf_map

	"""
	NAME:
		ln_epsilon
	PURPOSE:
		Run multinest, in order to find the lens parameters.
	INPUTS:
		array_param:					The lens parameters.
		nx, ny: 						Lens plane x and y size.
		nx_SP, ny_SP: 					Source plane x and y size (might be different from nx and ny)
		x_arcsec, y_arcsec:				Array of pixels (in arcsecond) in the LP.
		x_cen_arcsec, y_cen_arcsec:		Center of LP.
		pixel_scale, pixel_scale_SP:	LP and SP pixel scales (might be different).
		image_plane:					input lens plane
		area_image:						LP region that enters in the lensing process
		noise_map:						the noise map (pixel per pixel or a estimated mean level)
		psf_map:						the observed PSF
	OUTPUTS:
		ln_epsilon:						as in Eq. 5 of Dye et al., 2008, MNRAS
	"""

	# --------
	# PACKAGES
	import numpy as np
	from scipy.optimize import minimize
	from scipy.signal import fftconvolve
	from scipy.spatial import Voronoi, voronoi_plot_2d
	from Voronoi_tools import voronoi_finite_polygons_2d
	from sklearn.cluster import KMeans
	from astropy import units as u
	from tqdm import tqdm
	from RSI_tools import H_matrix_tesserae
	from image_tools import pos_arcsec2pixel, pos_arcsec2pixel, nint
	from lens_models import Singular_Isothermal_Ellipsoid as SIE
	# --------

#	=============
#	PARAMETER SET
	r_Ein = array_param[0]*u.arcsec
	q_lens, pos_angle = array_param[1], array_param[2]*u.deg
	dx_lens, dy_lens = array_param[3]*u.arcsec, array_param[4]*u.arcsec
	gamma, theta_shear = 0.0, 0.0*u.deg
#	=============

#	================
#	DEFLECTION ANGLE
#	coordinates of the lens
	x_lens_arcsec, y_lens_arcsec = x_cen_arcsec + dx_lens, y_cen_arcsec + dy_lens
#	Single Isothermal ELlipsoid
	(alphax_map, alphay_map, k_map, shear, mu) = SIE(r_Ein, q_lens, x_arcsec, y_arcsec, x_lens_arcsec, y_lens_arcsec, pos_angle, gamma, theta_shear)
#	==================

# 	===============
# 	PERFORM LENSING
	LENSED_x_arcsec, LENSED_y_arcsec = (x_arcsec - alphax_map), (y_arcsec - alphay_map)	# arcsec - position along the i-axis: THE WAY THE (FULL) LP IS DEFORMED
	LENSED_pixel_x, LENSED_pixel_y = pos_arcsec2pixel(LENSED_x_arcsec, pixel_scale), pos_arcsec2pixel(LENSED_y_arcsec, pixel_scale)	# pixel positions in the LENS PLANE
	LENSED_pixel_x_SP, LENSED_pixel_y_SP = pos_arcsec2pixel(LENSED_x_arcsec, pixel_scale_SP), pos_arcsec2pixel(LENSED_y_arcsec, pixel_scale_SP)	# pixel positions in the SOURCE PLANE
	# make sure to keep only the pixels in the source plane that are lensed into a visible region of the lens plane 
	ok = ((LENSED_pixel_x >= 0) & (LENSED_pixel_x <= nx-1) & (LENSED_pixel_y >= 0) & (LENSED_pixel_y <= ny-1)).nonzero()
	i1, i2 = ok[1], ok[0]	# index position in the LENS PLANE along the i-axis
# 	===============

#	===============================
#	SELECT AREA IN THE SOURCE PLANE
	area_source = np.zeros((ny_SP, nx_SP))
	for i in range(len(ok[0])):
		if area_image[i1[i],i2[i]] != 0.0:
			area_source[nint(LENSED_pixel_y_SP[i1[i],i2[i]]), nint(LENSED_pixel_x_SP[i1[i],i2[i]])] = 1.0

	GOOD_source = np.where(area_source == 1)
	ndim_source = len(GOOD_source[0])
	iy_source_cut, ix_source_cut = GOOD_source[0], GOOD_source[1]
	global GOOD_image
	GOOD_image = np.where(area_image == 1)
	global ndim_image
	ndim_image = len(GOOD_image[0])
	global iy_image_cut
	global ix_image_cut
	iy_image_cut, ix_image_cut = GOOD_image[0], GOOD_image[1]
#	===============================

#	========================
#	SOURCE PLANE TESSELATION
#	Number of Voronoi tesserae:
	if tass_number == 0: n_tass_SP = ndim_source
	else: n_tass_SP = tass_number
# 	Define the array to feed into the Voronoi routine:
	ok_LENSED_x_arcsec, ok_LENSED_y_arcsec = LENSED_x_arcsec[GOOD_image], LENSED_y_arcsec[GOOD_image]
	X, Y = (ok_LENSED_x_arcsec - x_cen_arcsec).value, (ok_LENSED_y_arcsec - y_cen_arcsec).value
	bad = np.isnan(X)
	if np.sum(bad) >= 1: X[bad] = 0.0
	bad = np.isnan(Y)
	if np.sum(bad) >= 1: Y[bad] = 0.0
#	Number of pixels actually used in the tesselation
	n_pt = np.size(X)
# 	Now identify the LP pixels associated with each Voronoi tesserae
	XY_SP_arcsec = np.zeros((n_pt,2))	# This will contain the pixel positions (in arcsec, because X & Y are in arcsec) in the source plane
	XY_SP_arcsec[:,0], XY_SP_arcsec[:,1] = X[:], Y[:]

	try:
		kmeans = KMeans(init='random', n_clusters=n_tass_SP, n_init=20)	# 'random’: choose k observations (rows) at random from data for the initial centroids.
		kmeans.fit(XY_SP_arcsec)
		centroids = kmeans.cluster_centers_
		labels = kmeans.labels_
		x_SP_centroids, y_SP_centroids = centroids[:,0], centroids[:,1]
		vor = Voronoi(centroids)
		regions_SP, vertices_SP = voronoi_finite_polygons_2d(vor)
	except: return -np.inf
# 	========================

#	==================
#	BUILDING F MATRIX
	global f_aux
	f, f_aux, sigma = np.zeros((ndim_image, n_tass_SP)), np.zeros((ndim_image, n_tass_SP)), np.zeros(ndim_image)
	d = image_plane[GOOD_image]/noise_map[GOOD_image]
	for i_vor in tqdm(range(n_tass_SP)):
		ok = labels == i_vor
		n_ok = np.sum(ok)
		image_plane_aux = np.zeros((ny,nx))
		image_plane_aux[iy_image_cut[ok], ix_image_cut[ok]] = 1.0
#		CONVOLUTION WITH THE PSF
#		To simulate an observation, we need to convolve with the instrumental PSF.
		image_plane_conv = fftconvolve(image_plane_aux, psf_map, mode = 'same')
		f[:, i_vor] = image_plane_conv[GOOD_image]/noise_map[GOOD_image]
		f_aux[:, i_vor] = image_plane_conv[GOOD_image]

#	matrix F (eq 3: F_ik = sum on J of (f_ij*f_kj)/sigma_j^2; sigma is implied in the f lower case matrix), and D (called c on eq 3)
	global F_matrix
	global D_matrix
	F_matrix = np.dot(np.transpose(f), f)
	D_matrix = np.dot(d,f)
#	==================

#	=====================
#	REGULARIZATION MATRIX (as in Eq.4 of Nightingale & Dye 2015)
	global H_matrix
	H_matrix = H_matrix_tesserae(n_tass_SP, regions_SP)
#	=====================

#	=================
#	BAYESIAN EVIDENCE
#	Third term (non lambda dependent part)) of Eq.(5) of Dye et al. (2008)
	global term3_part
	term3_part = np.linalg.slogdet(H_matrix)[1]
#	Fifth term in Eq.(5) of Dye et al. (2008)
	global term5
	term5 = np.sum(np.log(2.0*np.pi*noise_map[GOOD_image]**2.0))

# #	Test to see if the evidence has the expected shape
# 	import matplotlib.pyplot as plt
# 	o = 300
# 	lam = np.linspace(-15, 10, o)
# 	ev = np.zeros(o)
# 	for i in range(o):
# 		lam_aux = lam[i]
# 		ev[i] = bayesian_evidence(lam_aux)
# 	test = plt.figure()
# 	plt.plot(lam, -ev, 'rx')
# 	plt.plot(lam, -ev)
# 	test.savefig('ev.png')
# 	sys.exit()

#	Reg. term search through a minimization with the SLSQP routine
	x0 = np.array([-5.0]) # Always check with some iterations where, approximately, the x0 should be.
	bnds = [(-15.0, 15.0)] # Bnds keeps the solution inside this boundaries (it should, at least...)
	res = minimize(bayesian_evidence, x0, method='SLSQP', bounds = bnds, tol = 1e-8, options={'disp': False, 'maxiter': 80})
#	=================

#	Ok, but tomorrow never knows.
	if res.fun == 67101551:
		return -np.inf
	ln_epsilon = -res.fun
	return ln_epsilon



def ln_epsilon_uvplane(array_param, off_ra, off_dec, db_n_vis, x_arcsec, y_arcsec, pixel_scale, x_cen_arcsec, y_cen_arcsec, nx, ny, ny_SP, nx_SP, pixel_scale_SP, tass_number, \
		area_image, db_real_tab, db_imag_tab, db_wgt, db_uwave_tab, db_vwave_tab):
	global n_vis	
	global real_tab
	global imag_tab
	global wgt
	global uwave_tab
	global vwave_tab
	n_vis, real_tab, imag_tab, wgt, uwave_tab, vwave_tab = db_n_vis, db_real_tab, db_imag_tab, db_wgt, db_uwave_tab, db_vwave_tab
	"""
	NAME:
		ln_epsilon_uvplane
	PURPOSE:
		Run multinest
	INPUTS:
		array_param:					The lens parameters.
		nx, ny: 						Lens plane x and y size.
		nx_SP, ny_SP: 					Source plane x and y size (might be different from nx and ny)
		x_arcsec, y_arcsec:				Array of pixels (in arcsecond) in the LP.
		x_cen_arcsec, y_cen_arcsec:		Center of LP.
		pixel_scale, pixel_scale_SP:	LP and SP pixel scales (might be different).
		real/imag tab/wgt:				Visibilities
		uwave/vwave:					idem
		area_image:						IP region that enters in the lensing process
		primary_beam:					the image primary beam
	OUTPUTS:
		ln_epsilon:						as in Eq. 5 of Dye et al., 2008, MNRAS
	"""
	# --------
	# PACKAGES
	# --------
	import subprocess
	import numpy as np
	from scipy.optimize import minimize
	from scipy.spatial import Voronoi, voronoi_plot_2d
	from Voronoi_tools import voronoi_finite_polygons_2d
	from sklearn.cluster import KMeans
	from astropy import units as u
	from tqdm import tqdm
	from RSI_tools import H_matrix_tesserae
	from image_tools import pos_arcsec2pixel, pos_pixel2arcsec, nint
	from lens_models import Singular_Isothermal_Ellipsoid as SIE

#	==============
#	PARAMETERS SET
	r_Ein = array_param[0]*u.arcsec
	q_lens, pos_angle = array_param[1], array_param[2]*u.deg
	dx_lens, dy_lens = array_param[3]*u.arcsec - off_ra, array_param[4]*u.arcsec - off_dec # Eventual phase offset
	gamma, theta_shear = 0.0, 0.0*u.deg
#	==============

#	================
#	DEFLECTION ANGLE
#	coordinates of the lens
	x_lens_arcsec, y_lens_arcsec = x_cen_arcsec + dx_lens, y_cen_arcsec + dy_lens
#	SIE
	(alphax_map, alphay_map, k_map, shear, mu) = SIE(r_Ein, q_lens, x_arcsec, y_arcsec, x_lens_arcsec, y_lens_arcsec, pos_angle, gamma, theta_shear)
#	==================

#	===============
#	PERFORM LENSING
	LENSED_x_arcsec, LENSED_y_arcsec = (x_arcsec - alphax_map), (y_arcsec - alphay_map)
	LENSED_pixel_x, LENSED_pixel_y = pos_arcsec2pixel(LENSED_x_arcsec, pixel_scale), pos_arcsec2pixel(LENSED_y_arcsec, pixel_scale)
	LENSED_pixel_x_SP, LENSED_pixel_y_SP = pos_arcsec2pixel(LENSED_x_arcsec, pixel_scale_SP), pos_arcsec2pixel(LENSED_y_arcsec, pixel_scale_SP)
#	make sure to keep only the pixels in the source plane that are lensed into a visible region of the lens plane 
	ok = ((LENSED_pixel_x >= 0) & (LENSED_pixel_x <= nx-1) & (LENSED_pixel_y >= 0) & (LENSED_pixel_y <= ny-1)).nonzero()
	i1, i2 = ok[1], ok[0]	# index position in the LENS PLANE
#	===============

#	===============================
#	SELECT AREA IN THE SOURCE PLANE
	area_source = np.zeros((ny_SP, nx_SP))
	for i in range(np.size(ok)/2):
		if area_image[i1[i],i2[i]] != 0.0:
			area_source[nint(LENSED_pixel_y_SP[i1[i],i2[i]]), nint(LENSED_pixel_x_SP[i1[i],i2[i]])] = 1.0

	GOOD_source = (area_source == 1).nonzero()
	ndim_source = np.count_nonzero(area_source == 1)
	iy_source_cut, ix_source_cut = GOOD_source[0], GOOD_source[1]
	global GOOD_image
	GOOD_image = (area_image == 1).nonzero()
	global ndim_image
	global iy_image_cut
	global ix_image_cut
	ndim_image = np.count_nonzero(area_image == 1)
	iy_image_cut, ix_image_cut = GOOD_image[0], GOOD_image[1]
#	===============================

#	========================
#	SOURCE PLANE TESSELATION
#	Number of Voronoi tesserae:
	if tass_number == 0: n_tass_SP = ndim_source
	else: n_tass_SP = tass_number
# 	Define the array to feed into the Voronoi routine:
	ok_LENSED_x_arcsec, ok_LENSED_y_arcsec = LENSED_x_arcsec[GOOD_image], LENSED_y_arcsec[GOOD_image]
	X, Y = (ok_LENSED_x_arcsec - x_cen_arcsec).value, (ok_LENSED_y_arcsec - y_cen_arcsec).value
	bad = np.isnan(X)
	if np.sum(bad) >= 1: X[bad] = 0.0
	bad = np.isnan(Y)
	if np.sum(bad) >= 1: Y[bad] = 0.0
#	Number of pixels actually used in the tesselation
	n_pt = np.size(X)
# 	Now identify the LP pixels associated with each Voronoi tesserae
	XY_SP_arcsec = np.zeros((n_pt,2))	# This will contain the pixel positions (in arcsec, because X & Y are in arcsec) in the source plane
	XY_SP_arcsec[:,0], XY_SP_arcsec[:,1] = X[:], Y[:]

	try:
		kmeans = KMeans(init='random', n_clusters=n_tass_SP, n_init=20)	# 'random’: choose k observations (rows) at random from data for the initial centroids.
		kmeans.fit(XY_SP_arcsec)
		centroids = kmeans.cluster_centers_
		labels = kmeans.labels_
		x_SP_centroids, y_SP_centroids = centroids[:,0], centroids[:,1]
		vor = Voronoi(centroids)
		regions_SP, vertices_SP = voronoi_finite_polygons_2d(vor)
	except: return -np.inf
# 	========================

#	===================
#	BUILDING F MATRICES
	f_real = np.zeros((n_vis, n_tass_SP))
	f_imag = np.zeros((n_vis, n_tass_SP))
	global f_real_aux
	global f_imag_aux
	f_real_aux = np.zeros((n_vis, n_tass_SP))
	f_imag_aux = np.zeros((n_vis, n_tass_SP))

#	Mo ja faccio tassello per tassello. Quindi i_vor per i_vor, non i_source per i_source.
	for i_vor in tqdm(range(n_tass_SP)):
#		Since SP is unitary, and SB is conserved, I directly give value 1 to the LP intersted pixels.
		ok = labels == i_vor # Cioè: ok è dove labels coincide con i_vor.
		n_ok = np.sum(ok)
		image_plane_aux= np.zeros((ny,nx))
		image_plane_aux[iy_image_cut[ok], ix_image_cut[ok]] = 1.0#*primary_beam[iy_image_cut[ok], ix_image_cut[ok]]
		good = np.where(image_plane_aux > 0.0)
		n_good = len(good[0])
		if (n_good == 0):
			print('WARNING!: n_good = '), n_good
			return -np.inf
#		FOURIER transform
		x_delta_rad, y_delta_rad = (pos_pixel2arcsec(good[1], pixel_scale) - x_cen_arcsec).to('rad'), (pos_pixel2arcsec(good[0], pixel_scale) - y_cen_arcsec).to('rad')
		realFT_image_plane_aux, imagFT_image_plane_aux = np.zeros(n_vis), np.zeros(n_vis)
		for i_good in range(n_good):
			argument = -2.0*np.pi*(x_delta_rad[i_good]*uwave_tab + y_delta_rad[i_good]*vwave_tab)
			realFT_image_plane_aux += np.cos(argument)
			imagFT_image_plane_aux += np.sin(argument)

		f_real[:, i_vor] = realFT_image_plane_aux*np.sqrt(wgt)
		f_imag[:, i_vor] = imagFT_image_plane_aux*np.sqrt(wgt)
		f_real_aux[:, i_vor] = realFT_image_plane_aux
		f_imag_aux[:, i_vor] = imagFT_image_plane_aux

#	matrix F (eq 3: F_ik = sum on J of (f_ij*f_kj)/sigma_j^2, wich means a product between matrixes; sigma is implied in the f lower case matrix)
	F_real_matrix = np.dot(np.transpose(f_real), f_real)
	F_imag_matrix = np.dot(np.transpose(f_imag), f_imag)
	global F_matrix
	F_matrix = F_real_matrix + F_imag_matrix
#	matrix D  (called c on eq 3)
	d_real = real_tab*np.sqrt(wgt)
	d_imag = imag_tab*np.sqrt(wgt)
	D_real_matrix = np.dot(d_real,f_real)
	D_imag_matrix = np.dot(d_imag,f_imag)
	global D_matrix
	D_matrix = D_real_matrix + D_imag_matrix
#	===================

#	=====================
#	REGULARIZATION MATRIX (as in Eq.4 of Nightingale & Dye 2015)
	global H_matrix
	H_matrix = H_matrix_tesserae(n_tass_SP, regions_SP)
#	=====================

#	=================
#	BAYESIAN EVIDENCE
#	Third term (non lambda dependent part)) of Eq.(5) of Dye et al. (2008)
	global term3_part
	term3_part = np.linalg.slogdet(H_matrix)[1]
#	Fifth term in Eq.(5) of Dye et al. (2008)
	global term5
	term5 = np.sum(np.log(2.0*np.pi*(1/wgt)))

##	Test to see if the evidence has the expected shape
#	import matplotlib.pyplot as plt
#	o = 300
#	lam = np.linspace(-5, 15, o)
#	ev = np.zeros(o)
#	for i in tqdm(range(o)):
#		lam_aux = lam[i]
#		ev[i] = bayesian_evidence(lam_aux)
#	test = plt.figure()
#	plt.plot(lam, -ev, 'rx')
#	plt.plot(lam, -ev)
#	test.savefig('ev.png')
#	sys.exit()

#	Reg. term search through a minimization with the SLSQP routine
	x0 = np.array([5.0]) # Always check with some iterations where, approximately, the x0 should be.
	bnds = [(-3.0, 10.0)] # Bnds keeps the solution inside this boundaries (it should, at least...)
	res = minimize(bayesian_evidence_uvplane, x0, method='SLSQP', bounds = bnds, tol = 1e-8, options={'disp': False, 'maxiter': 80})
#	=================

#	Ok, but tomorrow never knows.
	if res.fun == 67101551:
		return -np.inf
	ln_epsilon = -res.fun
	return ln_epsilon

