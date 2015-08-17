import matplotlib.pyplot as plt
import pylab as pl
import pickle
import healpy as hp
import sys

from myclass import *
from mcmc import *
from data_model import *

# 
# # --------------------------------------
# class Component( MyClass ):
# 	def __init__(self, filename='', set_ones=False, scaling=1.):
# 		self.filename = filename
# 		if len(filename) == 0:
# 			if not set_ones:
# 				self.map = 0.
# 			else:
# 				self.map = 1. * scaling
# 		else:
# 			self.map = hp.read_map( filename ) * scaling
# 
# # --------------------------------------
# class FGData( MyClass ):
# 
# 	# ----------------------------------
# 	def __init__(self, selection=[], dir='/global/scratch2/sd/dpietrob/DM-haze/data/maps/',
# 		nside=128, DMmodel='3225'):
# 
# 		# ------------------------------
# 		class Galaxy( MyClass ):
# 			def __init__( self, galactic_templates, selection=[] ):
# 				if len(selection) == 0:
# 					selection = haze_templates.keys()
# 				for key in selection:
# 					file = galactic_templates[key]
# 					comp = Component( file )
# 					setattr( self, key, comp)
# 					
# 		# ------------------------------
# 		class Haze( MyClass ):
# 			def __init__( self, haze_templates, selection=[] ):
# 				if len(selection) == 0:
# 					selection = haze_templates.keys()
# 				for key in selection:
# 					file = haze_templates[key]
# 					comp = Component( file )
# 					setattr( self, key, comp)
# 					
# 		snside = str(nside).zfill(4)
# 		fgdir = dir+'ns'+snside+'/'
# 
# 		# --- Define templates
# 		galactic_templates = {}
# 		galactic_templates['cmb']      = fgdir + 'planck_hfi_ilc_ns'+snside+'_1deg_mK.fits'
# 		galactic_templates['dust']     = fgdir + 'map_fds_dust_94_ns'+snside+'_1deg.fits'
# 		galactic_templates['freefree'] = fgdir + 'map_halpha_ns'+snside+'_1deg_norm.fits'
# 		galactic_templates['sync']     = fgdir + 'map_haslam_dsds_ns'+snside+'_1deg_norm.fits'
# 		galactic_templates['dipole_x'] = fgdir + 'map_xdipole_ns'+snside+'.fits'
# 		galactic_templates['dipole_y'] = fgdir + 'map_ydipole_ns'+snside+'.fits'
# 		galactic_templates['dipole_z'] = fgdir + 'map_zdipole_ns'+snside+'.fits'
# 		galactic_templates['disc']     = fgdir + 'map_mydisc_ns'+snside+'_new.fits'
# 
# 		if ( (len(selection) == 0) or (selection[0].lower() == 'all') ):
# 			selection = ['cmb', 'dust', 'freefree', 'sync', 'dipole_x', 'dipole_y', 'dipole_z', 'disc']
# 
# # 		galaxy = Galaxy( selection )
# # 		galaxy.define( galactic_templates, selection )
# 		galaxy = Galaxy( galactic_templates, selection=selection )
# 		
# 		#print galaxy
# 
# 
# 		haze_templates = {}
# # 		haze_templates['single'] = fgdir + 'map_myhaze_ns'+snside+'.fits'
# 		#flist = ['fK', 'f030']
# # --- Generic template for all frequencies
# 		flist = ['generic']
# 		for item in flist:
# 			haze_templates[ item ]   = fgdir + 'map_myhaze_ns'+snside+'_new.fits'
# 
# # --- Tailored templateper frequency
# 		# DMmodel = '3225' 
# 		flist = ['K','030','Ka','Q','044','V','070','W']
# 		for item in flist:
# 			if DMmodel != 'myhaze':
# 				haze_templates[ 'f'+item ]   = fgdir + 'haze_model_54_'+DMmodel+'_ns'+snside+'_'+item+'_norm.fits'
# 			else:
# 				haze_templates[ 'f'+item ]   = fgdir + 'map_'+DMmodel+'_ns'+snside+'_new_'+item+'.fits'
# 
# # 		haze = Haze( frequency_list='haze_f'+flist)
# # 		haze.define( haze_templates )
# 		haze = Haze( haze_templates )
# 
# 		#print haze
# 		
# 		self.galaxy = galaxy
# 		self.haze = haze
# 
# # --------------------------------------
# class FreqData( MyClass ):
# 		
# 	# ----------------------------------
# 	def __init__(self, selection=['K','030','Ka','Q','044','V','070','W'],
# 		dir='/global/scratch2/sd/dpietrob/DM-haze/data/maps/',
# 		maskfile='ebvh_mask_top-half_disc_40_ns0128.fits',
# 		nside=128):
# 
# 		# maskfile='/global/scratch2/sd/dpietrob/DM-haze/data/maps/ns0128/ebvh_mask_ns0128.fits',
# 		# ------------------------------
# 		class Frequency( MyClass ):
# 			def __init__(self):
# 				self.imap = Component()
# 				self.irms = Component()
# 				self.cfreq = 0.
# 
# 			def define( self, datainfo ):
# 				# intensity
# 				imap = getattr( self, 'imap' )
# 				file = datainfo['ifilename']
# 				setattr( imap, 'filename', file)
# 				map = hp.read_map( file )
# 				setattr( imap, 'map', map)
# 				setattr( self, 'imap', imap)
# 
# 				# rms
# 				irms = getattr( self, 'irms' )
# 				file = datainfo['infilename']
# 				setattr( irms, 'filename', file)
# 				map = hp.read_map( file )
# 				setattr( irms, 'map', map)
# 				setattr( self, 'irms', irms)
# 
# 				# cfreq
# 				setattr( self, 'cfreq', datainfo['cfreq'])
# 
# 		frequencies = ['K','030','Ka','Q','044','V','070','W']
# 		cfreq = [22.8,28.4,33.,41.,44.,61.,70.4,94.]
# 		snside = 'ns' + str(nside).zfill(4)
# 		ddir = dir + snside
# 
# 		info = {}
# 		for ifreq,freq in enumerate(frequencies):
# 			datainfo = {}
# 			datainfo['ifilename'] = ddir + '/map_'+freq+'_'+snside+'_1deg_mK.fits' 
# 			datainfo['infilename'] = ddir + '/rms_adjst_'+freq+'_'+snside+'_1deg_mK.fits' 
# 			datainfo['cfreq'] = cfreq[ifreq]
# 
# 			info[freq] = datainfo
# 
# 		for freq in selection:
# 			fcomp = Frequency()
# 			fcomp.define( info[freq] )
# 			setattr(self, 'f'+freq, fcomp )
# 
# 		if len(maskfile) != 0:
# 			self.mask = Component( filename=ddir+'/'+maskfile )
# 		else:
# 			self.mask = Component( set_ones=True )
# 
# # --------------------------------------
# def fg_model( pars, nside=128l ):
# 	keys = pars.keys()
# 	for key in keys:
# 		if key[0] == 'f':
# 			freq = key.split('_')[0]
# 			break
# 	# print freq
# 	# print keys
# 	# nside = int( np.sqrt(npix/12) )
# 	npix = 12l * nside**2
# 	model = np.zeros( npix )
# 
# 	for key in keys:
# 		if key[0] == 'f':
# 			compname = '_'.join( key.split('_')[1:] )
# #			if ( (compname.lower() != 'monopole') and(compname.lower()[0:6] != 'dipole') ):
# 			if ( (compname.lower() != 'monopole') ):
# 				comp = getattr( gmodel.galaxy, compname )
# 				if len(comp.map) != npix:
# 					info_message('Inconsistent nside: '+str(nside)+', '+len(comp.map)+'. Abort')
# 					sys.exit()
# 				model += comp.map * pars[ key ]
# 			elif compname.lower() == 'monopole':
# 				model += pars[ key ]
# #			elif compname.lower() == 'dipole_x':
# #				model += make_dipole( nside, [pars[key],0.,0.] )
# #			elif compname.lower() == 'dipole_y':
# #				model += make_dipole( nside, [0.,pars[key],0.] )
# #			elif compname.lower() == 'dipole_z':
# #				model += make_dipole( nside, [0.,0.,pars[key]] )
# 			else:
# 				info_message( 'unknown parameter: '+key+'. Abort')
# 				sys.exit()			
# 		elif key[0:6] == 'haze_f':
# 			compname = 'generic'
# 			comp = getattr( gmodel.haze, compname )
# 			model += comp.map * pars[ key ]
# 		elif key == 'global_haze':
# 			comp = getattr( gmodel.haze, freq )
# 			model += comp.map * pars[ key ]
# 
# 	return model
# 
# # --------------------------------------
# def fg_chisq( pars ):
# 	keys = pars.keys()
# 	# print keys
# 	freq = []
# 	for key in keys:
# 		if key[0] == 'f':
# 			freq.append( key.split('_')[0] )
# 	frequencies = set( freq )
# 	# print frequencies
# 	nfreq = len( frequencies )
# 	
# 	chi2 = 0.e0
# 	
# 	for freq in frequencies:
# 		# print freq
# 		channel = getattr( data, freq )
# 		map = channel.imap.map
# 		rms = channel.irms.map
# 		npix = len(map)
# 		nside = int( np.sqrt(npix/12) )
# 		model = np.zeros( npix )
# 		# --- Get the model
# 		for key in keys:
# 			kfreq = key.split('_')[0]
# 			# print kfreq, freq
# 			#print key[0:len(freq)], freq
# 			if kfreq == freq:
# 				compname = '_'.join( key.split('_')[1:] )
# 				# print compname
# #				if ( (compname.lower() != 'monopole') and (compname.lower()[0:6] != 'dipole') ):
# 				if ( (compname.lower() != 'monopole') ):
# 					comp = getattr( gmodel.galaxy, compname )
# 					model += comp.map * pars[ key ]
# 				elif compname.lower() == 'monopole':
# 					model += pars[ key ]
# #				elif compname.lower() == 'dipole_x':
# #					model += make_dipole( nside, [pars[key],0.,0.] )
# #				elif compname.lower() == 'dipole_y':
# #					model += make_dipole( nside, [0.,pars[key],0.] )
# #				elif compname.lower() == 'dipole_z':
# #					model += make_dipole( nside, [0.,0.,pars[key]] )
# 				else:
# 					info_message( 'unknown parameter: '+key+'. Abort')
# 					sys.exit()
# 			elif key[0:6] == 'haze_f':
# 				compname = 'generic'
# 				comp = getattr( gmodel.haze, compname )
# 				model += comp.map * pars[ key ]
# 			elif key == 'global_haze':
# 				comp = getattr( gmodel.haze, freq )
# 				model += comp.map * pars[ key ]
# # When ever there's a frequency mismatch it prints.
# # 			else:
# # 				info_message( key+', '+freq)
# 			
# 		# hp.mollview(map, norm='hist')
# 		# pl.show()
# 		# hp.mollview(rms, norm='hist')
# 		# pl.show()
# 		# hp.mollview(map-model, norm='hist')
# 		# pl.show()
# 		c2 = np.sum( ( (map-model)/rms)**2 * data.mask.map )
# 		# print freq, c2/np.sum( data.mask.map )
# 		chi2 += c2
# 		#print chi2
# 		#sys.exit()
# 
# 	return chi2
# 
# # --------------------------------------
# def get_sed( pars ):
# 	keys = pars.keys()
# 	freq = []
# 	components = []
# 	for key in keys:
# 		tags = key.split('_')
# 		if tags[0][0] == 'f':
# 			if tags[0] not in freq:
# 				freq.append( tags[0] )
# 		if (len( tags )>1) and (tags[0] != 'global'):
# 			if '_'.join(tags[1:]) not in components:
# 				components.append( '_'.join(tags[1:]) )
# 	# print freq
# 	# print components
# 	
# 	nfreq = len(freq)
# 	ncomp = len(components)
# 
# 	sed = {}
# 	for j,c in enumerate(components):
# 		comp = []
# 		for i,f in enumerate(freq):
# 			key = f+'_'+c
# 			# print key
# 			comp.append( pars[key] )
# 		sed[c] = np.array(comp)
# 
# 	return sed
	
# --------------------------------------------------------------------

# frequencies = [ 'K' ,'030','Ka','Q','044','V','070'] #,'W' ]
# fs = ['fK','f030','fKa','fQ','f044','fV','f070'] #, 'fW']
# f = np.array([23., 28., 33, 41., 44., 61., 70.]) #, 93.] )

# frequencies = [ 'K' ]
# fs = ['fK']
# f = np.array([23.])

# cv = conversionfactor( f )
# print cv
# sys.exit()

selection = ['3225','3225h','3235','5114','12232','32232','myhaze']
model_list = [ 'myhaze' ]
# model_list = [ str( raw_input("Select a model among: ['3225','3225h','3235','5114','12232','32232','myhaze']\n-> : ") ) ]

print model_list
# sys.exit()
# model_list = ['3225','3225h','3235','5114','12232','32232']
# model_list = ['3225h','32232']
# model_list = ['5114']
# Templates converted to mK_CMB
# for model in model_list:
#  	make_dmhaze_templates( model=model )
# sys.exit()
dir = '/global/scratch2/sd/dpietrob/DM-haze/data/maps/ns0128/'

# --- Making a more suitable mask ---
# 30 deg radius from [0,-20]
radius = 20
if False:
	m = hp.read_map( dir + 'south_haze_ptsrc_mask_ns0128.fits' )
	# hp.mollview(m)
	print np.sum(m)/len(m)

	x,y,z = hp.ang2vec(np.pi/2+22.5/180.*np.pi ,0.)
	d = make_disc(128, radius, direction=[x,y,z])
	# hp.mollview(d)
	# hp.mollview(d*m)
	print np.sum(m*d)/len(m)

	hp.write_map( dir + 'south'+str(radius)+'_haze_ptsrc_mask_ns0128.fits', d*m*ptsrc*ptsrc1 )
	hp.gnomview( d*m*ptsrc*ptsrc1, rot=[0,-22.5], reso=25, min=-1, max=1)
	hp.graticule(dpar=5,dmer=5)
	# hp.gnomview( ptsrc, rot=[0,-22.5], reso=25, min=-1, max=1)
	# hp.graticule(dpar=5,dmer=5)

	m = hp.read_map( dir + 'fit_mask_p-35_35_t-10-35_ns0128.fits' )
	# hp.mollview(m)
	print np.sum(m)/len(m)

	pl.show()
# sys.exit()

tag = 'southern'+str(radius)	
# mfile = dir + 'south_haze_mask_ns0128.fits'
# mfile = dir + 'fit_mask_p-35_35_t-10-35_ns0128.fits'
m = hp.read_map( dir + 'ebvh_mask_ns0128.fits' )
d = hp.read_map( dir + 'southern_hemisphere_mask_ns0128.fits' )
m *= d
sel = 15
x,y,z = hp.ang2vec(np.pi/2+22.5/180.*np.pi ,0.)
# d = make_disc( 128, sel )
d = make_disc(128, sel, direction=[x,y,z])
m *= d
hp.write_map( dir + 'south_mask_'+str(sel)+'d.fits', m)
# mfile = dir + 'south'+str(radius)+'_haze_mask_ns0128.fits'
mfile = dir + 'south_mask_'+str(sel)+'d.fits'
mask = hp.read_map( mfile )
mask_npix = np.sum( mask )
fsky = np.sum( mask ) / len(mask)
gpix = np.array( (mask == 1.) )
bpix = np.array( (mask == 0.) )


# m = MCMC()

fix_cmb   = True
inc_disc  = True
inc_mono  = False
plot_gnom = True

find_minimum = True
# run_mcmc     = False
# analyze_mcmc = False
# run_2mcmc    = False
# ------
plot_min_sed = True
do_write = False
# ---
# plot_likes   = False
# plot_maps    = False
# plot_sed     = False
# analyze_cov  = False
# make_pics    = False

# data = FreqData( maskfile='south_haze_mask_ns0128.fits')

# ------ Run MCMC ------
for imodel, model in enumerate( model_list) :

	data = Data( DMmodel=model, maskfile='ebvh_mask_ns0128.fits', raw_sync=True )
	gmodel = data.gmodel
	runtag = model+tag
	print runtag

	freq_min = OrderedDict()
	# --- Set parameter values from full sky for some components ---
	if find_minimum:
		info_message( 'Running regression' )

# data = FreqData( maskfile='south_haze_mask_ns0128.fits')

		fx_mnpnt = OrderedDict()
		fx_mnpnt_er = OrderedDict()

		mnpnt = OrderedDict()
		mnpnt_er = OrderedDict()

		chisq = []

		c2 = 0.

		flux = {}
		cfreq = data.data.cfreq
		mKt2cgs = conversionfactor( cfreq, option='mKthermo2cgs' )
			
		haze_flux = []
		haze_flux_er = []
		haze_templ_flux = []
		haze_templ_flux_er = []

		# --- Find minimum for each frequency
		for ifreq, freq in enumerate( data.data.frequencies ):

			h = getattr( gmodel.haze, 'f'+freq )

			# --- Full run for starting point determination
			# Thin disc component removed: chisq gets worse, but global fit does
			# not change drmatically.
			# tfiles = [ gmodel.galaxy.cmb.filename, gmodel.galaxy.dust.filename, gmodel.galaxy.freefree.filename, 
			tfiles = [ gmodel.galaxy.cmb.filename, gmodel.galaxy.dust.filename,
				gmodel.galaxy.freefree.filename, 
				gmodel.galaxy.sync.filename, h.filename,
				dir+'map_monopole_ns0128.fits', dir+'map_xdipole_ns0128.fits',
				dir+'map_ydipole_ns0128.fits', dir+'map_zdipole_ns0128.fits' ]
			if inc_disc:
				tfiles.append( gmodel.galaxy.disc.filename )
				idisc = len(tfiles)-1

			map = getattr( data.data, 'f'+freq )
			r = map_regression( map.imap.map, tfiles, rms=map.irms.map,
				maskfile=data.mask.filename, return_res=True, return_chi2=True )
			chisq.append( str("%.3f" % r[3]) )
			c2 += r[3]

# 			hp.mollview( r[2]*data.mask.map, min=-0.025, max=0.025, title='Full sky residuals: '+freq )
# 			cl = hp.anafast( r[2]*data.mask.map, regression=True )
# 
# 			hp.mollview( gmodel.galaxy.cmb.map, min=-0.3, max=0.3, title='Planck HFI ILC' )
# 			cmb_cl = hp.anafast( gmodel.galaxy.cmb.map, regression=True )
# 
# 			wi = hp.read_map('../data/maps/ns0128/wmap_ilc_ns0128_1deg_mK.fits')
# 			wcmb_cl = hp.anafast( wi, regression=True )
# 			hp.mollview( gmodel.galaxy.cmb.map-wi, min=-0.1, max=0.1, title='Planck HFI ILC - WMAP9 ILC' )
# 
# 			l = np.array( range(3*128) )
# 			bl = hp.gauss_beam(np.pi/180.,3*128-1)
# 			wl = hp.read_cl( '/global/scratch2/sd/dpietrob/Software/Healpix_3.00/data/pixel_window_n0128.fits' )
# 			# print len(wl[0])
# 
# 			fig = plt.figure(30+ifreq, (11,8) )
# 			ax = plt.subplot(111)
# 			y = l[0:375]*(l[0:375]+1)/2./np.pi*cmb_cl[0:375]*1.e6/bl[0:375]**2/wl[0][0:375]**2
# 			plt.plot(l, l*(l+1)/2./np.pi*cl*1.e6/bl**2/wl[0][0:384]**2/fsky, '.-', label=freq+' residuals')
# 			plt.plot(l[0:375], y, 's-', color='r', label='Planck HFI ILC')
# 			plt.plot(l, l*(l+1)/2./np.pi*wcmb_cl*1.e6/bl**2, 'd-', color='m', label='WMAP9 ILC')
# 			plt.xlabel(r'$\ell$')
# 			plt.ylabel(r'$\mathcal{D}_\ell$ [$\mu K^2$]')
# 			ax.set_ylim([-500,6500])
# 			ax.set_xlim([1,3*128])
# 			plt.legend( loc='upper left',prop={'size':10})
# 			# pl.show()
# 			figtag = '../residual_spectra_'+model+'_'+freq
# 			fig.savefig(figtag+'.png')

			fx_mnpnt['f'+freq+'_cmb'] = r[0][0]
			fx_mnpnt['f'+freq+'_dust'] = r[0][1]
 			fx_mnpnt['f'+freq+'_freefree'] = r[0][2]
			fx_mnpnt['f'+freq+'_sync'] = r[0][3]
			fx_mnpnt['f'+freq+'_haze'] = r[0][4]
			if inc_disc:
				fx_mnpnt['f'+freq+'_disc'] = r[0][idisc]

			fx_mnpnt_er['f'+freq+'_cmb'] = r[1][0]
			fx_mnpnt_er['f'+freq+'_dust'] = r[1][1]
  			fx_mnpnt_er['f'+freq+'_freefree'] = r[1][2]
			fx_mnpnt_er['f'+freq+'_sync'] = r[1][3]
			fx_mnpnt_er['f'+freq+'_haze'] = r[1][4]
			if inc_disc:
				fx_mnpnt_er['f'+freq+'_disc'] = r[1][idisc]

			# --- Remove monopole, dipole and CMB from the map
			flat_comp = r[0][5] + \
				gmodel.galaxy.dipole_x.map * r[0][6] + \
				gmodel.galaxy.dipole_y.map * r[0][7] + \
				gmodel.galaxy.dipole_z.map * r[0][8]
			if fix_cmb:
				flat_comp += gmodel.galaxy.cmb.map * r[0][0]
				# gmodel.galaxy.freefree.map * r[0][2]

			# hp.mollview( map.imap.map, min=-0.725, max=0.725, title=freq )
			# hp.mollview( flat_comp, min=-0.25, max=0.25, title=freq )
			# hp.mollview( map.imap.map - flat_comp, min=-0.25, max=0.25, title=freq )

			# --- Run on smal patch after flat component removed
			tfiles = [ gmodel.galaxy.dust.filename,
				gmodel.galaxy.sync.filename,
				h.filename,
				gmodel.galaxy.freefree.filename ]

			if inc_disc:
				tfiles.append( gmodel.galaxy.disc.filename )
				idisc = len( tfiles ) -1
				
			if not fix_cmb:
				tfiles.append( gmodel.galaxy.cmb.filename )
				icmb = len( tfiles ) -1

			if inc_mono:
				tfiles.append( dir+'map_monopole_ns0128.fits' )
				imono = len( tfiles ) -1

			# map = getattr( data, 'f'+freq )
			# --- Setting up the smaller mask in the southern region
			r = map_regression( map.imap.map-flat_comp, tfiles, rms=map.irms.map,
				maskfile=mfile, return_res=True, return_chi2=True )
			chisq.append( str("%.3f" % r[3]) )
			c2 += r[3]

			if plot_gnom:
				if freq == '070' and False:
					ccc = data.gmodel.galaxy.cmb.map
					ccc[ bpix] = -1.6375e30
					hp.gnomview( ccc, min=-0.3, max=0.3, title='CMB', rot=[0,-22.5], reso=25 )
					hp.graticule(dpar=5,dmer=5)

				if freq == 'K' and False:
					ccc = data.gmodel.galaxy.sync.map
					ccc[ bpix] = -1.6375e30
					hp.gnomview( ccc, title='Sync', rot=[0,-22.5], reso=25 )
					hp.graticule(dpar=5,dmer=5)

					ccc = data.gmodel.galaxy.dust.map
					ccc[ bpix] = -1.6375e30
					hp.gnomview( ccc, title='Dust', rot=[0,-22.5], reso=25 )
					hp.graticule(dpar=5,dmer=5)

					ccc = data.gmodel.galaxy.freefree.map
					ccc[ bpix] = -1.6375e30
					hp.gnomview( ccc, title='FF', rot=[0,-22.5], reso=25 )
					hp.graticule(dpar=5,dmer=5)

				if False:
					ccc = map.imap.map
					ccc[ bpix] = -1.6375e30
					hp.gnomview( ccc, min=-0.3, title=freq,	rot=[0,-22.5], reso=25 )
					hp.graticule(dpar=5,dmer=5)

					ccc = r[2]*mask
					ccc[ bpix] = -1.6375e30
					hp.gnomview( ccc, min=-0.05, max=0.05, title='Residuals: '+freq,
						rot=[0,-22.5], reso=25 )
					hp.graticule(dpar=5,dmer=5)

				ccc = (r[2]+r[0][2]*h.map)*mask
				ccc[ bpix] = -1.6375e30
				hp.gnomview( ccc, title='Haze: '+freq, rot=[0,-22.5], reso=25, min=0., max=0.075 )
				hp.graticule(dpar=5,dmer=5)

				ccc = (r[2]+r[0][2]*h.map)*mask
				ccc[ bpix] = -1.6375e30
				hp.mollview( ccc, title='Haze: '+freq, rot=[0,-22.5], min=0., max=0.075 )
				hp.graticule(dpar=20,dmer=20)


				# ccc = (r[0][2]*h.map)*mask
				# ccc[ bpix] = -1.6375e30
				# hp.gnomview( ccc, min=0., title='Haze template: '+freq, rot=[0,-22.5], reso=25 )
				# hp.graticule(dpar=5,dmer=5)
			
			# --- Flux with residuals
			haze_flux.append( np.mean( (r[2] + r[0][2]*h.map )[gpix] ) * mKt2cgs[ifreq] )
			# haze_flux_er.append( np.std( (r[2] + r[0][2]*h.map)[gpix] ) * mKt2cgs[ifreq] )
			# --- Trying to be independent of the haze template
			haze_flux_er.append( np.std( r[2][gpix] ) * mKt2cgs[ifreq] )

			# print '\n'+freq
			# print np.mean( r[0][2]*h.map[gpix] )
			# print np.std( r[0][2]*h.map[gpix] )
			# print np.mean( r[2][gpix] )
			# print np.std( r[2][gpix] )
			# print haze_flux[-1]
			# print haze_flux_er[-1]
			# print np.sqrt( np.std( r[0][2]*h.map[gpix] )**2 + np.std( r[2][gpix] )**2 ) * mKt2cgs[ifreq]
			# print np.sqrt(np.var( h.map[gpix] ))


			# --- Flux from template
			flx = np.mean( r[0][2]*h.map[gpix] ) * mKt2cgs[ifreq]
			haze_templ_flux.append( flx )

			err_rel = r[1] / r[0]
			# haze_templ_flux_er.append( np.std( r[0][2]*h.map[gpix]) * mKt2cgs[ifreq] )
			haze_templ_flux_er.append( err_rel[2] * flx )

			freq_min[ freq ] = r

			mnpnt['f'+freq+'_dust'] = r[0][0]
			mnpnt['f'+freq+'_sync'] = r[0][1]
			mnpnt['f'+freq+'_haze'] = r[0][2]
			mnpnt['f'+freq+'_freefree'] = r[0][3]
			if inc_disc:
				mnpnt['f'+freq+'_disc'] = r[0][idisc]
			if not fix_cmb:
				mnpnt['f'+freq+'_cmb'] = r[0][icmb]
			if inc_mono:
				mnpnt['f'+freq+'_monopole'] = r[0][imono]

			# mnpnt['f'+freq+'_monopole'] = r[0][5]
			# mnpnt['f'+freq+'_dipole_x'] = r[0][6]
			# mnpnt['f'+freq+'_dipole_y'] = r[0][7]
			# mnpnt['f'+freq+'_dipole_z'] = r[0][8]

			mnpnt_er['f'+freq+'_dust'] = r[1][0]
			mnpnt_er['f'+freq+'_sync'] = r[1][1]
			mnpnt_er['f'+freq+'_haze'] = r[1][2]
			mnpnt_er['f'+freq+'_freefree'] = r[1][3]
			if inc_disc:
	 			mnpnt_er['f'+freq+'_disc'] = r[1][idisc]
			if not fix_cmb:
				mnpnt_er['f'+freq+'_cmb'] = r[1][icmb]
			if inc_mono:
	 			mnpnt_er['f'+freq+'_monopole'] = r[1][imono]

			# mnpnt_er['f'+freq+'_monopole'] = r[1][5]
			# mnpnt_er['f'+freq+'_dipole_x'] = r[1][6]
			# mnpnt_er['f'+freq+'_dipole_y'] = r[1][7]
			# mnpnt_er['f'+freq+'_dipole_z'] = r[1][8]

			# pl.show()
		# print haze_flux
		# print haze_templ_flux

		# print haze_flux_er
		# print haze_templ_flux_er

		flux[model+'_haze'] = ( np.array( haze_flux )*1.e20,
			np.array( haze_flux_er )*1.e20, #/np.sqrt(mask_npix), 
			np.array( haze_templ_flux )*1.e20,
			np.array( haze_templ_flux_er )*1.e20 ) #/np.sqrt(mask_npix) )

		if do_write:
			figtag = '../sed_dm_' + model + '_2step_' + tag
			if fix_cmb:
				figtag = figtag + '_fixCMB'
			if not inc_disc:
				figtag = figtag + '_noDISC'
			if inc_mono:
				figtag = figtag + '_MONO'

			file = open(figtag + '.txt', 'w')
			file.write( '# cfreq, res_haze, res_haze err, templ_haze, templ_haze err\n' )
			for ifreq,freq in enumerate(data.data.cfreq):
				# print ifreq, freq
				file.write( '%f \t %f \t %f \t %f \t %f \n' %(freq, flux[model+'_haze'][0][ifreq],
					flux[model+'_haze'][1][ifreq], flux[model+'_haze'][2][ifreq],
					flux[model+'_haze'][3][ifreq] ) )
			file.close()
		# print flux

		info_message(' total chi2 = '+str(c2)+','+str(c2*np.sum(mask)))

		sed = get_sed( mnpnt )
		sed_er = get_sed( mnpnt_er )

		fx_sed = get_sed( fx_mnpnt )
		fx_sed_er = get_sed( fx_mnpnt_er )
# 	
# 		sed_nod = get_sed( mnpnt_nod )
# 		sed_er_nod = get_sed( mnpnt_er_nod )
	
		if plot_min_sed:
			f = data.data.cfreq
			cv = conversionfactor( f, option='antenna2thermo' )
			fig = plt.figure(30+imodel,(18,7.5), dpi=80)

			# --- SED plot based on regression coefficients
			ax = plt.subplot(121)
	
			plt.errorbar( f, sed['dust'], yerr=sed_er['dust'], label='dust', color='b', fmt='s:' )
			plt.errorbar( f, sed['sync'], yerr=sed_er['sync'], label='sync', color='g', fmt='s:' )
			plt.errorbar( f, sed['haze'], yerr=sed_er['haze'], label='haze', color='m', fmt='s:' )
 			plt.errorbar( f, sed['freefree'], yerr=sed_er['freefree'], label='freefree', color='r', fmt='s:' )
			if not fix_cmb:
				plt.errorbar( f, sed['cmb'], yerr=sed_er['cmb'], label='cmb', color='k', fmt='s:' )
			if inc_disc:
				plt.errorbar( f, sed['disc'], yerr=sed_er['disc'], label='disc', color='y', fmt='s:' )

 			plt.errorbar( f, fx_sed['cmb'], yerr=fx_sed_er['cmb'], label='cmb FS', color='k', fmt='.-' )
  			plt.errorbar( f, fx_sed['dust'], yerr=fx_sed_er['dust'], label='dust FS', color='b', fmt='.-' )
  			plt.errorbar( f, fx_sed['sync'], yerr=fx_sed_er['sync'], label='sync FS', color='g', fmt='.-' )
  			plt.errorbar( f, fx_sed['haze'], yerr=fx_sed_er['haze'], label='haze FS', color='m', fmt='.-' )
  			plt.errorbar( f, fx_sed['freefree'], yerr=fx_sed_er['freefree'], label='freefree FS', color='r', fmt='.-' )
			if inc_disc:
	 			plt.errorbar( f, fx_sed['disc'], yerr=fx_sed_er['disc'], label='disc FS', color='y', fmt='s-' )

 
# 			plt.errorbar( f, sed_nod['dust'], yerr=sed_er_nod['dust'], label='dust no-disc', color='b', fmt='s-' )
# 			plt.errorbar( f, sed_nod['cmb'], yerr=sed_er_nod['cmb'], label='cmb no-disc', color='k', fmt='s-' )
# 			plt.errorbar( f, sed_nod['freefree'], yerr=sed_er_nod['freefree'], label='freefree no-disc', color='r', fmt='s-' )
# 			plt.errorbar( f, sed_nod['sync'], yerr=sed_er_nod['sync'], label='sync no-disc', color='g', fmt='s-' )
# 			plt.errorbar( f, sed_nod['haze'], yerr=sed_er_nod['haze'], label='haze no-disc', color='m', fmt='s-' )

			plt.plot( f, (f/f[0])**(-3.1)*cv*fx_sed['sync'][0], ':', label='-3.1', color='g' )
			plt.plot( f, (f/f[0])**(-2.15)*cv*fx_sed['freefree'][0], ':', label='-2.15', color='r' )
#			plt.plot( f, (f/f[0])**(-2.55)*cv*sed['haze'][0], ':', label='-2.55', color='m' )
			plt.plot( f, f*0.+sed['haze'][0], ':', color='m', label='same amplitude' )
			# plt.plot( f, (f/f[0])**(-.1)*cv*sed['haze'][0], ':', label='-2.65', color='m' )

			plt.title('Dark matter model: '+model)
			# plt.title('Dark matter model: '+model+' - Southern region')
			plt.xlabel(r'$\nu$ [GHz]')
			plt.ylabel('Template amplitudes')

			ax.set_yscale('log')
			ax.set_xscale('log')
			ax.set_ylim([5e-3,10])
			ax.set_xlim([8,100])
			plt.legend( loc='upper left',prop={'size':10})

			# --- Monopole dipole plot
# 			ax = plt.subplot(122)
# 	
# 			plt.errorbar( f, sed['monopole'], yerr=sed_er['monopole'], label='monopole', color='k' )
# 			plt.errorbar( f, sed['dipole_x'], yerr=sed_er['dipole_x'], label='dipole_x', color='b' )
# 			plt.errorbar( f, sed['dipole_y'], yerr=sed_er['dipole_y'], label='dipole_y', color='r' )
# 			plt.errorbar( f, sed['dipole_z'], yerr=sed_er['dipole_z'], label='dipole_z', color='m' )
# 			plt.plot( f, f*0., '--', color='g' )
# 			plt.errorbar( f, sed_noh['monopole'], yerr=sed_er_noh['monopole'], label='monopole no-haze', color='k', fmt='--' )
# 			plt.errorbar( f, sed_noh['dipole_x'], yerr=sed_er_noh['dipole_x'], label='dipole_x no-haze', color='b', fmt='--' )
# 			plt.errorbar( f, sed_noh['dipole_y'], yerr=sed_er_noh['dipole_y'], label='dipole_y no-haze', color='r', fmt='--' )
# 			plt.errorbar( f, sed_noh['dipole_z'], yerr=sed_er_noh['dipole_z'], label='dipole_z no-haze', color='m', fmt='--' )
# 
# 			plt.errorbar( f, sed_nod['monopole'], yerr=sed_er_nod['monopole'], label='monopole no-disc', color='k', fmt='s-' )
# 			plt.errorbar( f, sed_nod['dipole_x'], yerr=sed_er_nod['dipole_x'], label='dipole_x no-disc', color='b', fmt='s-' )
# 			plt.errorbar( f, sed_nod['dipole_y'], yerr=sed_er_nod['dipole_y'], label='dipole_y no-disc', color='r', fmt='s-' )
# 			plt.errorbar( f, sed_nod['dipole_z'], yerr=sed_er_nod['dipole_z'], label='dipole_z no-disc', color='m', fmt='s-' )
# 
# 			plt.title(r'$^{\rm red}\chi^2=$'+', '.join(chisq)+r'$\pm0.0028$'+' vs \n'+', '.join(chisq_noh)+' vs \n'+', '.join(chisq_nod), fontsize=12)
# 			# plt.title('Dark matter model: '+model+' - Southern region')
# 			plt.xlabel(r'$\nu$ [GHz]')
# 			plt.ylabel(r'Mono/dipole amplitudes [mK_CMB]')
# 
# 			ax.set_xlim([10,100])
# 			ax.set_xscale('log')
# 			hd = np.reshape([sed['monopole'],sed['dipole_x'],sed['dipole_y'],sed['dipole_z'],sed_noh['monopole'],sed_noh['dipole_x'],sed_noh['dipole_y'],sed_noh['dipole_z']], 8*len(sed['monopole']) )
# 			thrs = 1.1 * np.max( abs( hd ) )
# 			ax.set_ylim([-thrs,thrs])
# 			plt.legend( loc='lower right',prop={'size':10})

			# --- SED plot based on average flux in the region.
			# Mind, this is region dependent!!!
			ax = plt.subplot(122)
			# ff_sed = []
			# sync_sed = []
			# dust_sed = []
			# plt.plot( f, flux[model+'_haze']*10.e17, '.', color='m', label='Haze flux' )
 			plt.errorbar( f, flux[model+'_haze'][0], yerr=flux[model+'_haze'][1], label='Haze flux', color='m', fmt='s-' )

			# print mKt2cgs
			# plt.plot( f, (f/f[0])**(-.1)*flux[model+'_haze'][0][0], ':', label='-2.65', color='m' )
			# plt.plot( f, (f/f[0])**(-.5)*flux[model+'_haze'][0][0], ':', label='-3.05', color='m' )
			# plt.plot( f, (f/f[0])**(-1.)*flux[model+'_haze'][0][0], ':', label='-3.55', color='m' )
			plt.plot( f, f*0., '--', color='k')
			plt.plot( f, (f/f[0])**(-0.55)*flux[model+'_haze'][0][0], '--', color='m', label='-0.55')
			plt.plot( f, (f/f[0])**(-1.1)*flux[model+'_haze'][0][0], ':', color='m', label='-1.1')
 			# plt.errorbar( f, flux[model+'_haze'][2], yerr=flux[model+'_haze'][3], label='Haze template flux', color='m', fmt='.-' )
			# plt.errorbar( f, sed['haze'], yerr=sed_er['haze'], label='haze sed', color='m', fmt='s:' )
			# plt.errorbar( f, fx_sed['haze'], yerr=fx_sed_er['haze'], label='haze FS', color='m', fmt='d:' )
			# ax.set_yscale('log')
			# ax.set_xscale('log')
			ax.set_ylim([-.75,2])
			ax.set_xlim([8,100])
			plt.xlabel(r'$\nu$ [GHz]')
			plt.ylabel(r'Haze flux [erg/(cm$^2$ s Hz sr)] $\times 10^{20}$')
			plt.legend( loc='upper left',prop={'size':10})
			ax.annotate(', '.join(data.data.frequencies), xy=(10, -0.5) )
			ax.annotate(', '.join([ str("%0.3f" %s) for s in flux[model+'_haze'][0] ]), xy=(10, -0.6) )
			ax.annotate(r'$\pm$'+', '.join([ str("%0.3f" %s) for s in flux[model+'_haze'][1] ]), xy=(10, -0.7) )

			figtag = '../sed_dm_'+model+'_2step_'+tag
			if fix_cmb:
				figtag = figtag + '_fixCMB'
			if not inc_disc:
				figtag = figtag + '_noDISC'
			if inc_mono:
				figtag = figtag + '_MONO'

			fig.savefig(figtag + '.png')

		# --- Setting starting point
# 		p = OrderedDict()
# 
# 		limits = {}
# 		limits['dust'] = (0.,11.)
# 		limits['freefree'] = (0.,50)
# 		limits['sync'] = (0.,50)
# 		limits['cmb'] = (0.97,1.03)
# 		# limits['monopole'] = (-0.06,0.06)
# 		# limits['disc'] = (0,100)
# 		limits['monopole'] = (-0.06,0.06)
# 		limits['disc'] = (0.,1)
# 
# 		sig_scaling = 6.
# 		# p['global_haze'] = {'center_value':np.min( sed['haze'] ), 'sigma':np.min(sed_er['haze'])/sig_scaling, 'min':0., 'max':1000}
# 		p['global_haze'] = {'center_value':sed['haze'][0], 'sigma':np.min(sed_er['haze'])/sig_scaling, 'min':0., 'max':1000}
# 
# 		for freq in frequencies:
# 			r = freq_min[ freq ]
# 			p['f'+freq+'_cmb'] = {'center_value':np.max( np.array([0.,r[0][0]])), 'sigma':np.min(sed_er['cmb'])/sig_scaling, 'min':limits['cmb'][0], 'max':limits['cmb'][1]}
# 			p['f'+freq+'_dust'] = {'center_value':np.max( np.array([0.,r[0][1]])), 'sigma':np.min(sed_er['dust'])/sig_scaling, 'min':limits['dust'][0], 'max':limits['dust'][1]}
# 			p['f'+freq+'_freefree'] = {'center_value':np.max( np.array([0.,r[0][2]])), 'sigma':np.min(sed_er['freefree'])/sig_scaling, 'min':limits['freefree'][0], 'max':limits['freefree'][1]}
# 			p['f'+freq+'_sync'] = {'center_value':np.max( np.array([0.,r[0][3]])), 'sigma':np.min(sed_er['sync'])/sig_scaling, 'min':limits['sync'][0], 'max':limits['sync'][1]}
# 			p['f'+freq+'_monopole'] = {'center_value':r[0][5], 'sigma':np.min(sed_er['monopole'])/sig_scaling, 'min':limits['monopole'][0], 'max':limits['monopole'][1]}
# 			p['f'+freq+'_dipole_x'] = {'center_value':r[0][6], 'sigma':np.min(sed_er['dipole_x'])/sig_scaling, 'min':limits['monopole'][0], 'max':limits['monopole'][1]}
# 			p['f'+freq+'_dipole_y'] = {'center_value':r[0][7], 'sigma':np.min(sed_er['dipole_y'])/sig_scaling, 'min':limits['monopole'][0], 'max':limits['monopole'][1]}
# 			p['f'+freq+'_dipole_z'] = {'center_value':r[0][8], 'sigma':np.min(sed_er['dipole_z'])/sig_scaling, 'min':limits['monopole'][0], 'max':limits['monopole'][1]}
# 			p['f'+freq+'_disc'] = {'center_value':np.max( np.array([0.,r[0][9]])), 'sigma':np.min(sed_er['disc'])/sig_scaling, 'min':limits['disc'][0], 'max':limits['disc'][1]}
# 
# 		for key in p.keys():
# 			print key, '\t', p[key]['center_value'], '\t', p[key]['sigma']
# 
# 	# sys.exit()
# 	# ------ Run MCMC ------
# 	if run_mcmc:
# 		bf, ch = m.sampler( p, fg_chisq, nsamples=200000, seed=3, output_tag=runtag,
# 			temperature=1., accept_first=False, update_covmat=False )

if plot_min_sed:
	pl.show()

