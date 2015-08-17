import healpy as hp
import pylab as pl
from myclass import asinh10, make_disc
import numpy as np
import sys
from collections import OrderedDict

# ====================================================================
# To make Table.4 (previously Tab.3) with two step
if True:
	from twostep_fit_flux import twostep_fit_flux
	twostep_fit_flux('12132h', inc_bubbles=False, inc_haze=False, fsmaskf='ebvh_mask_ns0128.fits',
		raw_sync=True)

	#onestep_fit_flux( write_png=False, inc_disc=False, inc_haze=True, inc_bubbles=True, raw_sync=False,
	#	mfile='ebvh_mask_top-half_disc_40_ns0128.fits', model_list=['12132h'], plot_moll=False )
	sys.exit()



dir = '/global/scratch2/sd/dpietrob/DM-haze/data/maps/ns0128/'
nside = 128
snside = '0128'

# ====================================================================
# Make bubbles template as two discs for the haze/gamma-rays counterpart
# Not used: back to Gaussian for consistency and much difference after masking
if False:
	dirn = hp.ang2vec( (90.-22.5)/180.*np.pi, 0. )
	dn = make_disc( nside, 20., direction=dirn)
	dirs = hp.ang2vec( (90.+22.5)/180.*np.pi, 0. )
	ds = make_disc( nside, 20., direction=dirs)

	bubbles = (ds + dn ) * 0.01
	bubbles = hp.smoothing( bubbles, fwhm=10./180.*np.pi)
	hp.mollview( bubbles )
	hp.write_map( dir+'fermi_bubbles_ns'+snside+'.fits', bubbles)
	pl.show()

# ====================================================================
# Residual maps for the best fit: Tab3col1 first and 11113
if False:
	from myclass import *
	from mcmc import *
	from data_model import *

	dir = '/global/scratch2/sd/dpietrob/DM-haze/data/maps/ns0128/'
	frequencies = [ 'K', '030', 'Ka', 'Q', '044', 'V', '070']
	model = '12132h' #'11113'
	mfile="ebvh_mask_ns0128.fits" # Actual mask used in the mcmc for monopole/dipole termination. chisq doesn't match mcmc one...
#	mfile='ebvh_mask_top-half_disc_40_ns0128.fits'
	raw_sync = True #False used to get large scales
	# Residuals for the default case w/ mono/dipole
	inc_disc = False
#	inc_haze = False
	inc_bubbles = False
	inc_haze = True
	inc_mono = True
	inc_dipo = True
	fix_cmb = False
	regression = 0
	write_maps = False

	mcmc_region = hp.read_map(dir+"south_haze_ptsrc_mask_ns0128.fits") 
	data = Data( DMmodel=model, maskfile=mfile, raw_sync=raw_sync, frequencies=frequencies )
	gmodel = data.gmodel

	## DM+Bubbles
	## fname = '/global/scratch2/sd/dpietrob/DM-haze/results/fixed_mcmc_chain_bubbles_FINAL_long_'+model+'_2stp_08x.likestats'
	## DM only
	DMrec = True
	fname = '/global/scratch2/sd/dpietrob/DM-haze/results/fixed_mcmc_chain_DM_FINAL_long_'+model+'_2stp_08x.likestats'

	f = open( fname, 'r')
	mcmc_chisq = np.float( f.readline().split("=")[1] )
	f.close()

	bf, partag = np.loadtxt( fname, skiprows=3, usecols=(1,6), unpack=True, dtype=str)
	bf = bf.astype(float)
	dbf = OrderedDict(zip(partag,bf))
	#print dbf

	# --- Find minimum for each frequency
	trc2 = 0.
	for ifreq, freq in enumerate( data.data.frequencies ):

		tfiles = [ gmodel.galaxy.dust.filename,
			gmodel.galaxy.freefree.filename, 
			gmodel.galaxy.sync.filename]
		idust = 0
		iff   = 1
		isync = 2
		# ihaze = 3
		tags = ['f'+freq+'_{dust}','f'+freq+'_{freefree}','f'+freq+'_{sync}' ]

		if inc_haze:
			h = getattr( gmodel.haze, 'f'+freq )
			tfiles.append( h.filename )
			tags.append('f'+freq+'_mono')
			ihaze = len(tfiles)-1

		if inc_mono:
			tfiles.append( dir+'map_monopole_ns0128.fits' )
			imono = len(tfiles)-1
			tags.append('f'+freq+'_mono')
		if inc_dipo:
			tfiles.append( dir+'map_xdipole_ns0128.fits' )
			idipox = len(tfiles)-1

			tfiles.append( dir+'map_ydipole_ns0128.fits' )
			idipoy = len(tfiles)-1

			tfiles.append( dir+'map_zdipole_ns0128.fits' )
			idipoz = len(tfiles)-1
			tags.append( 'f'+freq+'_dipole_x' )
			tags.append( 'f'+freq+'_dipole_y' )
			tags.append( 'f'+freq+'_dipole_z' )
		if inc_disc:
			tfiles.append( gmodel.galaxy.disc.filename )
			idisc = len(tfiles)-1
		if inc_bubbles:
			tfiles.append( gmodel.galaxy.bubbles.filename )
			ibubbles = len(tfiles)-1
			tags.append('f'+freq+'_{bubbles}')
		if not fix_cmb:
			tfiles.append( gmodel.galaxy.cmb.filename )
			icmb = len(tfiles)-1
			tags.append('f'+freq+'_cmb')


		map = getattr( data.data, 'f'+freq )
		regmap = map.imap.map
		if fix_cmb:
			regmap -= gmodel.galaxy.cmb.map
		
		r = map_regression( regmap, tfiles, rms=map.irms.map,
			maskfile=data.mask.filename, regression=regression, 
			return_res=True, return_chi2=True )

#		print r[0]
		if write_maps:
			hp.write_map( 'Tab3col1_residuals_'+freq+'_ns0128.fits', r[2]*data.mask.map )

		# Jul 2015 changed tags to matches new label in getdist latex friendly
# 		tags = ['f'+freq+'_dust','f'+freq+'_freefree','f'+freq+'_sync',
# 			'f'+freq+'_mono','f'+freq+'_dipole_x','f'+freq+'_dipole_y','f'+freq+'_dipole_z',
# 			'f'+freq+'_cmb','f'+freq+'_bubbles','global_haze']
#		if inc_haze:
#			tags.pop(ihaze)
#			tfiles.pop(ihaze)
#			imono -= 1
#			idipox -= 1
#			idipoy -= 1
#			idipoz -= 1
#			ibubbles -= 1
#			icmb -= 1

		dbf['f'+freq+'_mono'] = r[0][imono]
		dbf['f'+freq+'_dipole_x'] = r[0][idipox]
		dbf['f'+freq+'_dipole_y'] = r[0][idipoy]
		dbf['f'+freq+'_dipole_z'] = r[0][idipoz]
		dbf['f'+freq+'_cmb'] = r[0][icmb]

#		tfiles.append( gmodel.galaxy.bubbles.filename )

		tags.append('global_{haze}')
		tfiles.append( getattr( gmodel.haze, 'f'+freq ).filename )

		if DMrec:
			if inc_bubbles:
				xx = tags.pop(ibubbles)
				xx = tfiles.pop(ibubbles)
			if inc_haze:
				xx = tags.pop(ihaze)
				xx = tfiles.pop(ihaze)
			
		dtf = OrderedDict(zip(tags,tfiles))
		print tfiles
		print dbf
		print dtf

		mcmc_model = np.zeros( 12*128**2, dtype=float )

		for k in dtf.keys():
			print k
			print dtf[k]
			print dbf[k]
			tmp = hp.read_map( dtf[k] )
			mcmc_model += tmp * dbf[k]
	
	#	hp.mollview( (regmap-mcmc_model)*data.mask.map )
	#	pl.show()
		if write_maps:
			hp.write_map( 'Tab3col4_wDM-BB_residuals_'+freq+'_ns0128.fits', (regmap-mcmc_model)*data.mask.map )

		rc2 = np.sum( ( (regmap-mcmc_model)*mcmc_region/map.irms.map )**2 ) / np.sum( mcmc_region )
		trc2 += rc2
		print r"$\chi_r^2=$"+"%5.2f" %rc2
		print r"$N_{pix}+$"+"%5.2f" %np.sum( mcmc_region )

	print r"total $\chi_r^2=$"+"%5.2f" %trc2
	print mcmc_chisq
	print r"total $\chi_{r,MCMC}^2=$"+"%5.2f" % (2.*mcmc_chisq/np.sum( mcmc_region ) )


	sys.exit()


# ====================================================================
# To make Table.3 and associated fig.7; run it three times with inc_haze/inc_bubbles on/off
from onestep_fit_flux import onestep_fit_flux
onestep_fit_flux( write_png=False, inc_disc=False, inc_haze=True, inc_bubbles=True, raw_sync=False,
	mfile='ebvh_mask_top-half_disc_40_ns0128.fits', model_list=['12132h'], plot_moll=False )
sys.exit()

# ====================================================================
maskfile = 'ebvh_mask_ns0128.fits'
mask = hp.read_map( dir + maskfile )
mask[mask == 0.] = hp.UNSEEN
hp.mollview(mask, cbar=False, title='Galactic mask')
hp.graticule(dpar=20,dmer=20)

reg_maskfile = 'south_haze_ptsrc_mask_ns0128.fits'
reg_mask = hp.read_map( dir + reg_maskfile )
reg_mask[reg_mask == 0.] = hp.UNSEEN
hp.mollview(reg_mask, title='Southern region mask')
hp.graticule(dpar=20,dmer=20)

pl.show()

sys.exit()

# ====================================================================
raw_sync = True 
files = {}
files['cmb']      = dir + 'planck_hfi_ilc_ns'+snside+'_1deg_mK.fits'
files['dust']     = dir + 'map_fds_dust_94_ns'+snside+'_1deg.fits'
files['freefree'] = dir + 'map_halpha_ns'+snside+'_1deg_norm.fits'
if not raw_sync:
	files['sync']     = dir + 'map_haslam_dsds_ns'+snside+'_1deg_norm.fits'
else:
	files['sync']     = dir + 'map_haslam_nofilt_ns'+snside+'_1deg_norm.fits'
# files['dipole_x'] = dir + 'map_xdipole_ns'+snside+'.fits'
# files['dipole_y'] = dir + 'map_ydipole_ns'+snside+'.fits'
# files['dipole_z'] = dir + 'map_zdipole_ns'+snside+'.fits'
files['disc']     = dir + 'map_mydisc_ns'+snside+'_new.fits'
files['haze']     = dir + 'map_myhaze_ns'+snside+'_new.fits'

mx_ar = {}
mx_ar['freefree'] = 0.1
mx_ar['sync'] = 0.1
mx_ar['dust'] = 0.1
mx_ar['disc'] = 1.
mx_ar['haze'] = 1.

for key in files.keys():
	map = hp.read_map(files[key])
	if key == 'cmb':
		mn = -0.3
		mx = 0.3
	else:
		mn = 0.
		mx = mx_ar[ key ]
	hp.mollview( map, title=key, min=mn, max=mx, unit='mK' )
	hp.graticule(dpar=20,dmer=20)
	

pl.show()
sys.exit()

