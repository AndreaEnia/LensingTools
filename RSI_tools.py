#!/usr/bin/python
# -*- coding: latin-1 -*-

"""
NAME:
	H_matrix_tesserae
PURPOSE:
	Compute the H matrix in the RSI formalism, on Voronoi tesselation, as in Eq.4 of Nightingale & Dye 2015
REVISION HISTORY:
	Written		A. Enia		November  2017
"""
def H_matrix_tesserae(n_tass_SP, regions_SP):
	import numpy as np

	# =====================
	# REGULARIZATION MATRIX (as in Eq.4 of Nightingale & Dye 2015)
	# =====================
	H_matrix = np.zeros((n_tass_SP,n_tass_SP))
	# ---------------------
	# LOOP over the TASSELS to find those with at least a vertex in common
	# ---------------------
	index_vor = np.arange(n_tass_SP)
	for i_vor_ref in np.arange(n_tass_SP):
		region_ref = regions_SP[i_vor_ref]
		incommon = [filter(lambda x: x in region_ref, sublist) for sublist in regions_SP]
		dim = np.array([[np.shape(incommon[i])[0]] for i in np.arange(n_tass_SP)])
		dim = np.array(list(zip(*dim)[0]))
		ok = np.logical_and(dim > 0, index_vor <> i_vor_ref)
		n_neighboors = np.sum(ok)
		H_matrix[i_vor_ref, i_vor_ref] = H_matrix[i_vor_ref, i_vor_ref] + n_neighboors
		H_matrix[ok, ok] = H_matrix[ok, ok] + 1.0 
		H_matrix[ok,i_vor_ref] = H_matrix[ok,i_vor_ref] - 1.0
		H_matrix[i_vor_ref,ok] = H_matrix[i_vor_ref,ok] - 1.0

	return H_matrix

"""
NAME:
	Inversion
PURPOSE:
	Excute the regularized semilinear inversion (for square grid)
REVISION HISTORY:
	Written		A. Enia		October  2015
	Modified	A. Enia		November 2017
"""
def inversion(image_plane, F_matrix, lam_best, H_matrix, D_matrix, f_aux, ndim_image, ndim_source, ix_source_cut, iy_source_cut, ix_image_cut, iy_image_cut, nx, ny, nx_SP, ny_SP):
	import numpy as np
	import matplotlib.pyplot as plt

	# Linear inversion for lambda = lam_best
	Finv_matrix = np.linalg.inv(F_matrix + lam_best*H_matrix)
	# Result: matrix S 
	S_matrix = np.dot(D_matrix, Finv_matrix)
	# Matrix describing the lens plane
	#	this is the reconstructed image using the most likely reconstructed source (S_matrix), the f matrix (the one that mapped
	#	a unitary source plane inside the ellipses to the lens plane and is linked to the lensing potential) and PSF models
	#	confined to the annular region that maps on to the source plane
	Image_matrix = np.zeros(ndim_image)
	for i_image in range(ndim_image):
		Image_matrix[i_image] = np.sum(S_matrix*f_aux[i_image, :])
	# Covariance matrix (see Warren & Dye (2003), Eq.(17)).
       	error_matrix = np.zeros(ndim_source)
       	R_matrix = Finv_matrix
       	RH_matrix = np.dot(R_matrix, H_matrix)
       	for i_source in range(ndim_source):
	       	error_matrix[i_source] = R_matrix[i_source,i_source] - lam_best*np.sum(R_matrix[:, i_source]*RH_matrix[:, i_source])
	# Reconstructed source plane and the minimum chi_sq lens plane (here labeled as "reconstructed")
	source_plane_reconstructed = np.zeros((ny_SP, nx_SP))
	error_source_plane_reconstructed = np.zeros((ny_SP, nx_SP))
	for i_source in range(ndim_source):
		source_plane_reconstructed[iy_source_cut[i_source], ix_source_cut[i_source]] = S_matrix[i_source]
		error_source_plane_reconstructed[iy_source_cut[i_source], ix_source_cut[i_source]] = error_matrix[i_source]
	image_plane_reconstructed = np.zeros((ny, nx))
	residuals = np.zeros((ny,nx))
	for i_image in range(ndim_image):
		image_plane_reconstructed[iy_image_cut[i_image], ix_image_cut[i_image]] = Image_matrix[i_image]
		residuals[iy_image_cut[i_image], ix_image_cut[i_image]] = \
					image_plane_reconstructed[iy_image_cut[i_image], ix_image_cut[i_image]] - image_plane[iy_image_cut[i_image], ix_image_cut[i_image]]

	return source_plane_reconstructed, error_source_plane_reconstructed, image_plane_reconstructed, residuals

"""
NAME:
	adgrid_inversion
PURPOSE:
	Excute the regularized semilinear inversion (for adaptive grid)
REVISION HISTORY:
	Written		A. Enia		October  2015
	Modified	A. Enia		November 2017
"""
def adgrid_inversion(image_plane, F_matrix, lam_best, H_matrix, D_matrix, f_aux, ndim_image, ndim_source, ix_image_cut, iy_image_cut, nx, ny):
	import numpy as np

	# Linear inversion for lambda = lam_best
	Finv_matrix = np.linalg.inv(F_matrix + lam_best*H_matrix)
	# Result: matrix S 
	S_matrix = np.dot(D_matrix, Finv_matrix)
	# Covariance matrix (see Warren & Dye (2003), Eq.(17)).
	error_matrix = np.zeros(ndim_source)
	R_matrix = Finv_matrix
	RH_matrix = np.dot(R_matrix,H_matrix)
	for i_source in range(ndim_source):
		error_matrix[i_source] = R_matrix[i_source,i_source] - lam_best*np.sum(R_matrix[i_source, :]*RH_matrix[i_source, :])
	error_matrix = np.sqrt(error_matrix)
	# Matrix describing the image plane
	#	this is the reconstructed image using the most likely reconstructed source (S_matrix), the f matrix (the one that mapped
	#	a unitary source plane inside the ellipses to the image plane and is linked to the lensing potential) and PSF models
	#	confined to the annular region that maps on to the source plane
	Image_matrix = np.zeros(ndim_image)
	for i_image in range(ndim_image):
		Image_matrix[i_image] = np.sum(S_matrix*f_aux[i_image, :])

	image_plane_reconstructed = np.zeros((ny, nx))
	residuals = np.zeros((ny,nx))
	for i_image in range(ndim_image):
		image_plane_reconstructed[iy_image_cut[i_image], ix_image_cut[i_image]] = Image_matrix[i_image]
		residuals[iy_image_cut[i_image], ix_image_cut[i_image]] = \
					image_plane_reconstructed[iy_image_cut[i_image], ix_image_cut[i_image]] - image_plane[iy_image_cut[i_image], ix_image_cut[i_image]]

	norm_S_matrix = 1.0 - (S_matrix - S_matrix.min())/(S_matrix.max() - S_matrix.min())
	return norm_S_matrix, S_matrix, error_matrix, image_plane_reconstructed, residuals

"""
NAME:
	adgrid_inversion_interf
PURPOSE:
	Excute the regularized semilinear inversion (for adaptive grid, interferometric code)
REVISION HISTORY:
	Written		A. Enia		Febraury 2016
	Modified	A. Enia		November 2017
"""
def adgrid_inversion_uvplane(F_matrix, lam_best, H_matrix, D_matrix, n_tass_SP):
	import numpy as np

	# Linear inversion for lambda = lam_best
	Finv_matrix = np.linalg.inv(F_matrix + lam_best*H_matrix)
	# Result: matrix S 
	S_matrix = np.dot(D_matrix, Finv_matrix)
	norm_S_matrix = 1.0 - (S_matrix - S_matrix.min())/(S_matrix.max() - S_matrix.min())
	# Covariance matrix (see Warren & Dye (2003), Eq.(17)).
	error_matrix = np.zeros(S_matrix.shape[0])
	R_matrix = np.linalg.inv(F_matrix + lam_best*H_matrix)
	RH_matrix = np.dot(R_matrix,H_matrix)
	for i_source in range(n_tass_SP):
		error_matrix[i_source] = R_matrix[i_source,i_source] - lam_best*np.sum(R_matrix[i_source, :]*RH_matrix[i_source, :])
	error_matrix = np.sqrt(error_matrix)
	# Signal-to-noise counts on SP tesserae
	S2N_matrix = S_matrix/error_matrix

	return norm_S_matrix, S_matrix, error_matrix, S2N_matrix
