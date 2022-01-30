def SIE_glafic(folder, image_params, output_lens_params):
	'''
	NAME:
		SIE_glafic
	INPUT:
		ACHTUNG: pass only the values (pippo.value)

		folder			name of the output folder
		output_lens_params 	as follows:
						0 is the lens velocity dispersion
						1 is the ellipticity
						2 is the position angle
						3 is the x position of the lens with respect to the image center
						4 is the y position of the lens with respect to the image center
						5 is the external shear intensity
						6 is the external shear position angle

		image_params 		as follows:
						0 is the number of image pixels on the x-axis
						1 is the number of image pixels on the y-axis
						2 is the pixel scale
						3 is the source redshift
						4 is the lens redshift

	PURPOSE:
		Find the caustics and critics of a lensing phenomena via glafic, for a power-law profile
	REVISION HISTORY:
		Written		A. Enia		October  2015
	'''
	import numpy as np
	import subprocess
	from astropy.cosmology import Planck13 as cosmo

	# ===============
	# LENS PARAMETERS
	type_lens             = 'sie'
	sigma_lens_glafic     = output_lens_params[0]
	dx_lens_arcsec_glafic = output_lens_params[3]
	dy_lens_arcsec_glafic = output_lens_params[4]
	ell_lens_glafic       = 1.0 - output_lens_params[1]
	theta_lens_glafic     = - (90.0 - output_lens_params[2])     # (NOTE THE SIGN!)
	core_lens_glafic      = 0.0
	# ===============

	# ===============
	# EXTERNAL SHEAR PARAMETERS
	ext_shear_glafic = output_lens_params[5]
	theta_shear_glafic = output_lens_params[6]
	# ===============

	# =================
	# SOURCE PARAMETERS
	type_source             = 'gauss' 
	flux_source_glafic      = 10.0
	dx_source_arcsec_glafic = 0.0
	dy_source_arcsec_glafic = 0.0
	ell_source_glafic       = 1.0 - 0.75 
	theta_source_glafic     = 0.0
	sigma_source_glafic     = 0.4 
	# =================

	# =====================
	# INPUT FILE FOR GLAFIC

	# LENSES
	n_lenses   = 1
	lens_type  = np.chararray(n_lenses, itemsize = 5)
	lens_sigma = np.zeros(n_lenses)
	lens_x     = np.zeros(n_lenses)
	lens_y     = np.zeros(n_lenses)
	lens_ell   = np.zeros(n_lenses)
	lens_rot   = np.zeros(n_lenses)
	lens_core  = np.zeros(n_lenses)
	#	- lens 1
	lens_type[0]  = type_lens
	lens_sigma[0] = sigma_lens_glafic 
	lens_x[0]     = dx_lens_arcsec_glafic      
	lens_y[0]     = dy_lens_arcsec_glafic     
	lens_ell[0]   = ell_lens_glafic
	lens_rot[0]   = theta_lens_glafic    
	lens_core[0]  = core_lens_glafic    

	# EXTENDED SOURCES 
	n_extsources = 1
	source_type  = np.chararray(n_extsources, itemsize = 5)
	source_z     = np.zeros(n_extsources)
	source_flux  = np.zeros(n_extsources)
	source_x     = np.zeros(n_extsources)
	source_y     = np.zeros(n_extsources)
	source_ell   = np.zeros(n_extsources)
	source_rot   = np.zeros(n_extsources)
	source_sigma = np.zeros(n_extsources)
	#	- source 1
	source_type[0]  = type_source
	source_z[0]     = image_params[3]
	source_flux[0]  = flux_source_glafic    # mJy
	source_x[0]     = dx_source_arcsec_glafic
	source_y[0]     = dy_source_arcsec_glafic
	source_ell[0]   = ell_source_glafic
	source_rot[0]   = theta_source_glafic     # deg anticlockwise from North
	source_sigma[0] = sigma_source_glafic

	# POINT SOURCES
	n_ptsources = 0

	# create the input file for glafic
	xmin = - image_params[0]*image_params[2]/2.0
	xmax = + image_params[0]*image_params[2]/2.0
	ymin = - image_params[1]*image_params[2]/2.0
	ymax = + image_params[1]*image_params[2]/2.0
	pix_ext = image_params[2] # arcsec (pixel scale of the image)
	pix_poi = image_params[2] # arcsec (used for point sources) DON'T BOTHER: I AM USING AN EXTENDED SOURCE ANYWAY
	maxlev = 6

	with open(folder+'/glafic.input', 'w') as f:
		f.write('## setting primary parameters\n')
		aux = 'omega ', cosmo.Om0, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'lambda ', 1 - cosmo.Om0, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'weos ', cosmo.w(0), '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'hubble ', cosmo.h, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# redshift of the lens
		aux = 'zl ', image_params[4], '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# prefix for output files
		aux = 'prefix ', folder+'/glafic', '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# size and pixel scale of the image
		aux = 'xmin ', xmin, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'ymin ', ymin, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'xmax ', xmax, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'ymax ', ymax, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'pix_ext ', pix_ext, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'pix_poi ', pix_poi, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'maxlev ', maxlev, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# define number of lenses, extsources, ptsources
		aux = 'startup', n_lenses+1, n_extsources, n_ptsources, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# define lenses
		for i_lens in range(n_lenses):
			aux = 'lens', lens_type[i_lens], lens_sigma[i_lens], lens_x[i_lens], lens_y[i_lens], lens_ell[i_lens], \
					lens_rot[i_lens], lens_core[i_lens], 0.0, '\n'
			aux = ' '.join((map(str, aux)))
			f.write(aux)
	# ext shear (only one)
		aux = 'lens pert', source_z[0], 0.0, 0.0, ext_shear_glafic, theta_shear_glafic, 0.0, 0.0, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# define extsources
		for i_source in range(n_extsources):
			aux = 'extend', source_type[i_source], source_z[i_source], source_flux[i_source], source_x[i_source], source_y[i_source], source_ell[i_source], \
					source_rot[i_source], source_sigma[i_source], 0.0, '\n'
			aux = ' '.join((map(str, aux)))
			f.write(aux)
	# define point source
		for i_source in range(n_ptsources-1):
			aux = 'point', source_z[i_source], source_flux[i_source], source_x[i_source], source_y[i_source], '\n'
			aux = ' '.join((map(str, aux)))
			f.write(aux)

		f.write('end_startup \n')
	# execute commands
		f.write('start_command \n')
		aux = 'writelens', source_z[0], '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'writeimage_ori', 0.0, 0.0, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'writeimage', 0.0, 0.0, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# write critical lines and caustics
		aux = 'writecrit', source_z[0], '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# terminate
		f.write('quit \n')

	# RUN GLAFIC
	print
	print(' run GLAFIC ... ')
	subprocess.call("glafic "+folder+"/glafic.input", shell=True)
	print(' ... done!')
	print

	return;


def SPLE_glafic(folder, image_params, output_lens_params):
	'''
	NAME:
		SPLE_glafic
	INPUT:
		ACHTUNG: pass only the values (pippo.value)

		folder			name of the output folder
		output_lens_params 	as follows:
						0 is the lens Einstein radius
						1 is the ellipticity
						2 is the position angle
						3 is the x position of the lens with respect to the image center
						4 is the y position of the lens with respect to the image center
						5 is the external shear intensity
						6 is the external shear position angle
						7 is the slope of the power law ellipsoid
		image_params 		as follows:
						0 is the number of image pixels on the x-axis
						1 is the number of image pixels on the y-axis
						2 is the pixel scale
						3 is the source redshift
						4 is the lens redshift

	PURPOSE:
		Find the caustics and critics of a lensing phenomena via glafic, for a power-law profile
	REVISION HISTORY:
		Written		A. Enia		October  2015
		PLE profile	M. Negrello	January  2016
	'''
	import numpy as np
	import subprocess
	from astropy.cosmology import Planck13 as cosmo

	# ===============
	# LENS PARAMETERS
	type_lens             = 'pow'
	norm_lens_glafic      = output_lens_params[0] # Einstein radius in arcsec
	dx_lens_arcsec_glafic = output_lens_params[3]
	dy_lens_arcsec_glafic = output_lens_params[4]
	ell_lens_glafic       = 1.0 - output_lens_params[1]
	slope_lens_glafic     = output_lens_params[7] + 1.0
	theta_lens_glafic     = output_lens_params[2] - 90.0	# (NOTE THE SIGN!)
	core_lens_glafic      = 0.0
	# ===============

	# ===============
	# EXTERNAL SHEAR PARAMETERS
	ext_shear_glafic = output_lens_params[5]
	theta_shear_glafic = output_lens_params[6] - 90.0 # in degrees
	# ===============

	# =================
	# SOURCE PARAMETERS
	type_source             = 'gauss' 
	flux_source_glafic      = 10.0
	dx_source_arcsec_glafic = 0.0
	dy_source_arcsec_glafic = 0.0
	ell_source_glafic       = 1.0 - 0.75 
	theta_source_glafic     = 0.0
	sigma_source_glafic     = 0.4 
	# =================

	# =====================
	# INPUT FILE FOR GLAFIC

	# LENSES
	n_lenses   = 1
	lens_type  = np.chararray(n_lenses, itemsize = 5)
	lens_norm = np.zeros(n_lenses)
	lens_x     = np.zeros(n_lenses)
	lens_y     = np.zeros(n_lenses)
	lens_slope = np.zeros(n_lenses)
	lens_ell   = np.zeros(n_lenses)
	lens_rot   = np.zeros(n_lenses)
	lens_core  = np.zeros(n_lenses)
	#	- lens 1
	lens_type[0]  = 'pow'
	lens_norm[0]  = norm_lens_glafic 
	lens_x[0]     = dx_lens_arcsec_glafic      
	lens_y[0]     = dy_lens_arcsec_glafic  
	lens_slope[0] = slope_lens_glafic   
	lens_ell[0]   = ell_lens_glafic
	lens_rot[0]   = theta_lens_glafic    
	lens_core[0]  = core_lens_glafic    

	# EXTENDED SOURCES 
	n_extsources = 1
	source_type  = np.chararray(n_extsources, itemsize = 5)
	source_z     = np.zeros(n_extsources)
	source_flux  = np.zeros(n_extsources)
	source_x     = np.zeros(n_extsources)
	source_y     = np.zeros(n_extsources)
	source_ell   = np.zeros(n_extsources)
	source_rot   = np.zeros(n_extsources)
	source_sigma = np.zeros(n_extsources)
	#	- source 1
	source_type[0]  = type_source
	source_z[0]     = image_params[3]
	source_flux[0]  = flux_source_glafic    # mJy
	source_x[0]     = dx_source_arcsec_glafic
	source_y[0]     = dy_source_arcsec_glafic
	source_ell[0]   = ell_source_glafic
	source_rot[0]   = theta_source_glafic     # deg anticlockwise from North
	source_sigma[0] = sigma_source_glafic

	# POINT SOURCES
	n_ptsources = 0

	# create the input file for glafic
	xmin = - image_params[0]*image_params[2]/2.0
	xmax = + image_params[0]*image_params[2]/2.0
	ymin = - image_params[1]*image_params[2]/2.0
	ymax = + image_params[1]*image_params[2]/2.0
	pix_ext = image_params[2] # arcsec (pixel scale of the image)
	pix_poi = image_params[2] # arcsec (used for point sources) DON'T BOTHER: I AM USING AN EXTENDED SOURCE ANYWAY
	maxlev = 5

	with open(folder+'/glafic.input', 'w') as f:
		f.write('## setting primary parameters\n')
		aux = 'omega ', cosmo.Om0, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'lambda ', 1 - cosmo.Om0, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'weos ', cosmo.w(0), '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'hubble ', cosmo.h, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# redshift of the lens
		aux = 'zl ', image_params[4], '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# prefix for output files
		aux = 'prefix ', folder+'/glafic', '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# size and pixel scale of the image
		aux = 'xmin ', xmin, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'ymin ', ymin, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'xmax ', xmax, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'ymax ', ymax, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'pix_ext ', pix_ext, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'pix_poi ', pix_poi, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'maxlev ', maxlev, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# define number of lenses, extsources, ptsources
		aux = 'startup', n_lenses+1, n_extsources, n_ptsources, '\n' # + 1 because of external shear!!!
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# define lenses
		for i_lens in range(n_lenses):
			aux = 'lens pow', source_z[0], lens_x[i_lens], lens_y[i_lens], lens_ell[i_lens], \
					lens_rot[i_lens], lens_norm[i_lens], lens_slope[i_lens], '\n'
			aux = ' '.join((map(str, aux)))
			f.write(aux)
	# ext shear (only one)
		aux = 'lens pert', source_z[0], 0.0, 0.0, ext_shear_glafic, theta_shear_glafic, 0.0, 0.0, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# define extsources
		for i_source in range(n_extsources):
			aux = 'extend gauss', source_z[i_source], source_flux[i_source], source_x[i_source], source_y[i_source], source_ell[i_source], \
					source_rot[i_source], source_sigma[i_source], 0.0, '\n'
			aux = ' '.join((map(str, aux)))
			f.write(aux)
	# define point source
		for i_source in range(n_ptsources-1):
			aux = 'point', source_z[i_source], source_flux[i_source], source_x[i_source], source_y[i_source], '\n'
			aux = ' '.join((map(str, aux)))
			f.write(aux)

		f.write('end_startup \n')
	# execute commands
		f.write('start_command \n')
		aux = 'writelens', source_z[0], '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'writeimage_ori', 0.0, 0.0, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'writeimage', 0.0, 0.0, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# write critical lines and caustics
		aux = 'writecrit', source_z[0], '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# terminate
		f.write('quit \n')

	# RUN GLAFIC
	print
	print(' run GLAFIC ... ')
	subprocess.call("glafic "+folder+"/glafic.input", shell=True)
	print(' ... done!')
	print

	return;

def SPLE_glafic_twolenses(folder, image_params, output_lens1_params, output_lens2_params):
	'''
	NAME:
		SPLE_glafic_twolenses
	INPUT:
		ACHTUNG: pass only the values (pippo.value)

		folder			name of the output folder
		output_lens_params 	as follows:
						0 is the lens Einstein radius
						1 is the ellipticity
						2 is the position angle
						3 is the x position of the lens with respect to the image center
						4 is the y position of the lens with respect to the image center
						5 is the external shear intensity
						6 is the external shear position angle
						7 is the slope of the power law ellipsoid
		image_params 		as follows:
						0 is the number of image pixels on the x-axis
						1 is the number of image pixels on the y-axis
						2 is the pixel scale
						3 is the source redshift
						4 is the lens redshift

	PURPOSE:
		Find the caustics and critics of a lensing phenomena via glafic, for two power-law profiles
	REVISION HISTORY:
		Written		A. Enia		October  2015
		PLE profile	M. Negrello	January  2016
	'''
	import numpy as np
	import subprocess
	from astropy.cosmology import Planck13 as cosmo

	# ===============
	# LENS PARAMETERS
	type_lens1             = 'pow'
	norm_lens1_glafic      = output_lens1_params[0] # Einstein radius in arcsec
	dx_lens1_arcsec_glafic = output_lens1_params[3]
	dy_lens1_arcsec_glafic = output_lens1_params[4]
	ell_lens1_glafic       = 1.0 - output_lens1_params[1]
	slope_lens1_glafic     = output_lens1_params[7] + 1.0
	theta_lens1_glafic     = output_lens1_params[2] - 90.0	# (NOTE THE SIGN!)
	core_lens1_glafic      = 0.0

	type_lens2             = 'pow'
	norm_lens2_glafic      = output_lens2_params[0] # Einstein radius in arcsec
	dx_lens2_arcsec_glafic = output_lens2_params[3]
	dy_lens2_arcsec_glafic = output_lens2_params[4]
	ell_lens2_glafic       = 1.0 - output_lens2_params[1]
	slope_lens2_glafic     = output_lens2_params[7] + 1.0
	theta_lens2_glafic     = output_lens2_params[2] - 90.0	# (NOTE THE SIGN!)
	core_lens2_glafic      = 0.0
	# ===============

	# ===============
	# EXTERNAL SHEAR PARAMETERS
	ext_shear_glafic = output_lens1_params[5]
	theta_shear_glafic = output_lens1_params[6] - 90.0 # in degrees
	# ===============

	# =================
	# SOURCE PARAMETERS
	type_source             = 'gauss' 
	flux_source_glafic      = 10.0
	dx_source_arcsec_glafic = 0.0
	dy_source_arcsec_glafic = 0.0
	ell_source_glafic       = 1.0 - 0.75 
	theta_source_glafic     = 0.0
	sigma_source_glafic     = 0.4 
	# =================

	# =====================
	# INPUT FILE FOR GLAFIC

	# LENSES
	n_lenses   = 2
	lens_type  = np.chararray(n_lenses, itemsize = 5)
	lens_norm = np.zeros(n_lenses)
	lens_x     = np.zeros(n_lenses)
	lens_y     = np.zeros(n_lenses)
	lens_slope = np.zeros(n_lenses)
	lens_ell   = np.zeros(n_lenses)
	lens_rot   = np.zeros(n_lenses)
	lens_core  = np.zeros(n_lenses)
	#	- lens 1
	lens_type[0]  = type_lens1
	lens_norm[0]  = norm_lens1_glafic 
	lens_x[0]     = dx_lens1_arcsec_glafic      
	lens_y[0]     = dy_lens1_arcsec_glafic  
	lens_slope[0] = slope_lens1_glafic   
	lens_ell[0]   = ell_lens1_glafic
	lens_rot[0]   = theta_lens1_glafic    
	lens_core[0]  = core_lens1_glafic    
	#	- lens 2
	lens_type[1]  = type_lens2
	lens_norm[1]  = norm_lens2_glafic 
	lens_x[1]     = dx_lens2_arcsec_glafic      
	lens_y[1]     = dy_lens2_arcsec_glafic  
	lens_slope[1] = slope_lens2_glafic   
	lens_ell[1]   = ell_lens2_glafic
	lens_rot[1]   = theta_lens2_glafic    
	lens_core[1]  = core_lens2_glafic    

	# EXTENDED SOURCES 
	n_extsources = 1
	source_type  = np.chararray(n_extsources, itemsize = 5)
	source_z     = np.zeros(n_extsources)
	source_flux  = np.zeros(n_extsources)
	source_x     = np.zeros(n_extsources)
	source_y     = np.zeros(n_extsources)
	source_ell   = np.zeros(n_extsources)
	source_rot   = np.zeros(n_extsources)
	source_sigma = np.zeros(n_extsources)
	#	- source 1
	source_type[0]  = type_source
	source_z[0]     = image_params[3]
	source_flux[0]  = flux_source_glafic    # mJy
	source_x[0]     = dx_source_arcsec_glafic
	source_y[0]     = dy_source_arcsec_glafic
	source_ell[0]   = ell_source_glafic
	source_rot[0]   = theta_source_glafic     # deg anticlockwise from North
	source_sigma[0] = sigma_source_glafic

	# POINT SOURCES
	n_ptsources = 0

	# create the input file for glafic
	xmin = - image_params[0]*image_params[2]/2.0
	xmax = + image_params[0]*image_params[2]/2.0
	ymin = - image_params[1]*image_params[2]/2.0
	ymax = + image_params[1]*image_params[2]/2.0
	pix_ext = image_params[2] # arcsec (pixel scale of the image)
	pix_poi = image_params[2] # arcsec (used for point sources) DON'T BOTHER: I AM USING AN EXTENDED SOURCE ANYWAY
	maxlev = 5

	with open(folder+'/glafic.input', 'w') as f:
		f.write('## setting primary parameters\n')
		aux = 'omega ', cosmo.Om0, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'lambda ', 1 - cosmo.Om0, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'weos ', cosmo.w(0), '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'hubble ', cosmo.h, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# redshift of the lens
		aux = 'zl ', image_params[4], '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# prefix for output files
		aux = 'prefix ', folder+'/glafic', '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# size and pixel scale of the image
		aux = 'xmin ', xmin, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'ymin ', ymin, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'xmax ', xmax, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'ymax ', ymax, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'pix_ext ', pix_ext, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'pix_poi ', pix_poi, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'maxlev ', maxlev, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# define number of lenses, extsources, ptsources
		aux = 'startup', n_lenses+1, n_extsources, n_ptsources, '\n' # + 1 because of external shear!!!
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# define lenses
		for i_lens in range(n_lenses):
			aux = 'lens', lens_type[i_lens], source_z[0], lens_x[i_lens], lens_y[i_lens], lens_ell[i_lens], \
					lens_rot[i_lens], lens_norm[i_lens], lens_slope[i_lens], '\n'
			aux = ' '.join((map(str, aux)))
			f.write(aux)
	# ext shear (only one)
		aux = 'lens pert', source_z[0], 0.0, 0.0, ext_shear_glafic, theta_shear_glafic, 0.0, 0.0, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# define extsources
		for i_source in range(n_extsources):
			aux = 'extend', source_type[i_source], source_z[i_source], source_flux[i_source], source_x[i_source], source_y[i_source], source_ell[i_source], \
					source_rot[i_source], source_sigma[i_source], 0.0, '\n'
			aux = ' '.join((map(str, aux)))
			f.write(aux)
	# define point source
		for i_source in range(n_ptsources-1):
			aux = 'point', source_z[i_source], source_flux[i_source], source_x[i_source], source_y[i_source], '\n'
			aux = ' '.join((map(str, aux)))
			f.write(aux)

		f.write('end_startup \n')
	# execute commands
		f.write('start_command \n')
		aux = 'writelens', source_z[0], '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'writeimage_ori', 0.0, 0.0, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
		aux = 'writeimage', 0.0, 0.0, '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# write critical lines and caustics
		aux = 'writecrit', source_z[0], '\n'
		aux = ' '.join((map(str, aux)))
		f.write(aux)
	# terminate
		f.write('quit \n')

	# RUN GLAFIC
	print
	print(' run GLAFIC ... ')
	subprocess.call("glafic "+folder+"/glafic.input", shell=True)
	print(' ... done!')
	print

	return

def Gauss2D(x_cen_pixel, y_cen_pixel, sigma_x_pixel, sigma_y_pixel, nx, ny):
    import numpy as np
#   array of pixel positions:
    pixel_x, pixel_y = np.zeros(shape=(nx,ny)), np.zeros(shape=(nx,ny))
    index_x, index_y= np.arange(nx), np.arange(ny)
    for iy in index_y: pixel_x[:,iy] = index_x*1.0    
    for ix in index_x: pixel_y[ix,:] = index_y*1.0
#   Build Gaussian:
    Gauss2D = 1.0/(2.0*np.pi*sigma_x_pixel*sigma_y_pixel)*np.exp(-0.5*((pixel_x - x_cen_pixel)**2.0/sigma_x_pixel**2.0 + (pixel_y - y_cen_pixel)**2.0/sigma_y_pixel**2.0))
    return Gauss2D

def simulate_source(FWHM_x_arcsec, FWHM_y_arcsec, FWHM_dx_cen, FWHM_dy_cen, nx_SP, ny_SP, pixel_scale_SP):
    import numpy as np
    x_cen_arcsec, y_cen_arcsec = nx_SP*pixel_scale_SP/2.0, ny_SP*pixel_scale_SP/2.0
    FWHM_x_cen, FWHM_y_cen = (x_cen_arcsec + FWHM_dx_cen)/pixel_scale_SP, (y_cen_arcsec + FWHM_dy_cen)/pixel_scale_SP
    FWHM_x_pixel, FWHM_y_pixel = FWHM_x_arcsec/pixel_scale_SP, FWHM_y_arcsec/pixel_scale_SP
    sigma_x_pixel, sigma_y_pixel = FWHM_x_pixel/(2.0*np.sqrt(2.0*np.log(2.0))), FWHM_y_pixel/(2.0*np.sqrt(2.0*np.log(2.0)))
    return Gauss2D(FWHM_y_cen, FWHM_x_cen, sigma_y_pixel, sigma_x_pixel, nx_SP, ny_SP)
