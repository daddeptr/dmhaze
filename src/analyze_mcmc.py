import matplotlib.pyplot as plt
import pylab as pl
import pickle
import healpy as hp
import sys

from myclass import *
from mcmc import *
from data_model import *

# --------------------------------------
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
# --------------------------------------
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
# --------------------------------------
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
# --------------------------------------
# def fg_model( gmodel, pars, nside=128l ):
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

# --------------------------------------
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

# --------------------------------------
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

# freq_tags = []
# for f in frequencies:
# 	freq_tags.append( 'f' + f )

# model_list = ['3225','3225h','3235','5114','12232','32232','myhaze']
model_list = ['3225','3225h','3235','5114','12232','32232']

print model_list

tag = '_Lmask_fullFreq_x6sig'	

m = MCMC()

analyze_mcmc = True
run_2mcmc    = False
# ------
plot_min_sed = False
# ---
plot_likes   = True
plot_maps    = True
plot_sed     = True
analyze_cov  = True
make_pics    = True

# data = FreqData( maskfile='south_haze_mask_ns0128.fits')
# data = FreqData()

# ------ Run MCMC ------
# for imodel,model in enumerate( model_list) :
# 	gmodel = FGData( DMmodel=model )
# 	runtag = model+tag
# 	print runtag
# # 	print gmodel.haze
# # 	for f in fs:
# # 		x = getattr(gmodel.haze, f)
# # 		hp.mollview(x.map, title=f)
# # 	pl.show()
# # 	sys.exit()
# 
# 	freq_min = OrderedDict()
# 	if find_minimum:
# 		# --- Set parameter values ---
# 		info_message( 'Running regression' )
# 
# 		dir = '/global/scratch2/sd/dpietrob/DM-haze/data/maps/ns0128/'
# 
# 		mnpnt = OrderedDict()
# 		mnpnt_er = OrderedDict()
# 		mnpnt_noh = OrderedDict()
# 		mnpnt_er_noh = OrderedDict()
# 
# 		chisq = []
# 		chisq_noh = []
# 
# 		c2 = 0.
# 		for freq in frequencies:
# 			h = getattr( gmodel.haze, 'f'+freq )
# 			tfiles = [ gmodel.galaxy.cmb.filename, gmodel.galaxy.dust.filename, gmodel.galaxy.freefree.filename, 
# 				gmodel.galaxy.sync.filename, h.filename, dir+'map_monopole_ns0128.fits', dir+'map_xdipole_ns0128.fits',
# 				dir+'map_ydipole_ns0128.fits', dir+'map_zdipole_ns0128.fits', gmodel.galaxy.disc.filename]
# 
# 			map = getattr( data, 'f'+freq )
# 			r = map_regression( map.imap.map, tfiles, rms=map.irms.map, maskfile=data.mask.filename, return_res=True, return_chi2=True )
# 			chisq.append( str("%.3f" % r[3]) )
# 			c2 += r[3]
# 			# hp.mollview(r[2]*data.mask.map, min=-0.05, max=0.05, title=freq)
# 			# hp.graticule( dpar=20, dmer=20 )
# 
# 			freq_min[ freq ] = r
# 
# 			mnpnt['f'+freq+'_cmb'] = r[0][0]
# 			mnpnt['f'+freq+'_dust'] = r[0][1]
# 			mnpnt['f'+freq+'_freefree'] = r[0][2]
# 			mnpnt['f'+freq+'_sync'] = r[0][3]
# 			mnpnt['f'+freq+'_haze'] = r[0][4]
# 			mnpnt['f'+freq+'_disc'] = r[0][9]
# 			mnpnt['f'+freq+'_monopole'] = r[0][5]
# 			mnpnt['f'+freq+'_dipole_x'] = r[0][6]
# 			mnpnt['f'+freq+'_dipole_y'] = r[0][7]
# 			mnpnt['f'+freq+'_dipole_z'] = r[0][8]
# 
# 			mnpnt_er['f'+freq+'_cmb'] = r[1][0]
# 			mnpnt_er['f'+freq+'_dust'] = r[1][1]
# 			mnpnt_er['f'+freq+'_freefree'] = r[1][2]
# 			mnpnt_er['f'+freq+'_sync'] = r[1][3]
# 			mnpnt_er['f'+freq+'_haze'] = r[1][4]
# 			mnpnt_er['f'+freq+'_monopole'] = r[1][5]
# 			mnpnt_er['f'+freq+'_dipole_x'] = r[1][6]
# 			mnpnt_er['f'+freq+'_dipole_y'] = r[1][7]
# 			mnpnt_er['f'+freq+'_dipole_z'] = r[1][8]
# 			mnpnt_er['f'+freq+'_disc'] = r[1][9]
# 
# 			tfiles = [ gmodel.galaxy.cmb.filename, gmodel.galaxy.dust.filename, gmodel.galaxy.freefree.filename, 
# 				gmodel.galaxy.sync.filename, dir+'map_monopole_ns0128.fits', dir+'map_xdipole_ns0128.fits',
# 				dir+'map_ydipole_ns0128.fits', dir+'map_zdipole_ns0128.fits', gmodel.galaxy.disc.filename]
# 
# 			map = getattr( data, 'f'+freq )
# 			r = map_regression( map.imap.map, tfiles, rms=map.irms.map, maskfile=data.mask.filename, return_res=True, return_chi2=True )
# 			chisq_noh.append( str("%.3f" % r[3]) )
# 
# 			mnpnt_noh['f'+freq+'_cmb'] = r[0][0]
# 			mnpnt_noh['f'+freq+'_dust'] = r[0][1]
# 			mnpnt_noh['f'+freq+'_freefree'] = r[0][2]
# 			mnpnt_noh['f'+freq+'_sync'] = r[0][3]
# 			mnpnt_noh['f'+freq+'_haze'] = r[0][4]
# 			mnpnt_noh['f'+freq+'_monopole'] = r[0][4]
# 			mnpnt_noh['f'+freq+'_dipole_x'] = r[0][5]
# 			mnpnt_noh['f'+freq+'_dipole_y'] = r[0][6]
# 			mnpnt_noh['f'+freq+'_dipole_z'] = r[0][7]
# 			mnpnt_noh['f'+freq+'_disc'] = r[0][8]
# 
# 			mnpnt_er_noh['f'+freq+'_cmb'] = r[1][0]
# 			mnpnt_er_noh['f'+freq+'_dust'] = r[1][1]
# 			mnpnt_er_noh['f'+freq+'_freefree'] = r[1][2]
# 			mnpnt_er_noh['f'+freq+'_sync'] = r[1][3]
# 			mnpnt_er_noh['f'+freq+'_monopole'] = r[1][4]
# 			mnpnt_er_noh['f'+freq+'_dipole_x'] = r[1][5]
# 			mnpnt_er_noh['f'+freq+'_dipole_y'] = r[1][6]
# 			mnpnt_er_noh['f'+freq+'_dipole_z'] = r[1][7]
# 			mnpnt_er_noh['f'+freq+'_disc'] = r[1][8]
# 
# 
# 
# 		# chisq = np.array(chisq)
# 		# chisq_noh = np.array(chisq_noh)
# 
# 		info_message(' total chi2 = '+str(c2)+','+str(c2*np.sum(data.mask.map)))
# 
# 		sed = get_sed( mnpnt )
# 		sed_er = get_sed( mnpnt_er )
# 
# 		sed_noh = get_sed( mnpnt_noh )
# 		sed_er_noh = get_sed( mnpnt_er_noh )
# 	
# 		if plot_min_sed:
# 			cv = conversionfactor( f, option='antenna2thermo' )
# 			fig = plt.figure(10+imodel,(19,7.5))
# 			ax = plt.subplot(121)
# 	
# 			plt.errorbar( f, sed['dust'], yerr=sed_er['dust'], label='dust', color='b' )
# 			plt.errorbar( f, sed['cmb'], yerr=sed_er['cmb'], label='cmb', color='k' )
# 			plt.errorbar( f, sed['freefree'], yerr=sed_er['freefree'], label='freefree', color='r' )
# 			plt.errorbar( f, sed['sync'], yerr=sed_er['sync'], label='sync', color='g' )
# 			plt.errorbar( f, sed['haze'], yerr=sed_er['haze'], label='haze', color='m' )
# 			plt.errorbar( f, sed['disc'], yerr=sed_er['disc'], label='disc', color='y' )
# 
# 			plt.errorbar( f, sed_noh['dust'], yerr=sed_er_noh['dust'], label='dust no-haze', color='b', fmt='--' )
# 			plt.errorbar( f, sed_noh['cmb'], yerr=sed_er_noh['cmb'], label='cmb no-haze', color='k', fmt='--' )
# 			plt.errorbar( f, sed_noh['freefree'], yerr=sed_er_noh['freefree'], label='freefree no-haze', color='r', fmt='--' )
# 			plt.errorbar( f, sed_noh['sync'], yerr=sed_er_noh['sync'], label='sync no-haze', color='g', fmt='--' )
# 			plt.errorbar( f, sed_noh['disc'], yerr=sed_er_noh['disc'], label='disc no-haze', color='y', fmt='--' )
# 
# 			plt.plot( f, (f/f[0])**(-3.1)*cv*sed['sync'][0], ':', label='-3.1', color='g' )
# 			plt.plot( f, (f/f[0])**(-2.15)*cv*sed['freefree'][0], ':', label='-2.15', color='r' )
# 			plt.plot( f, f*0.+sed['haze'][0], '--', color='m', label='same amplitude' )
# 
# 			plt.title('Dark matter model: '+model)
# 			# plt.title('Dark matter model: '+model+' - Southern region')
# 			plt.xlabel(r'$\nu$ [GHz]')
# 			plt.ylabel('Template amplitudes')
# 
# 			ax.set_yscale('log')
# 			ax.set_xscale('log')
# 			ax.set_ylim([1e-2,100])
# 			ax.set_xlim([8,100])
# 			plt.legend( loc='upper left',prop={'size':10})
# 
# 			ax.annotate(', '.join(frequencies), xy=(10, 0.03) )
# 			ax.annotate(', '.join([ str("%0.3f" %s) for s in sed['haze'] ]), xy=(10, 0.02) )
# 			ax.annotate(r'$\pm$'+', '.join([ str("%0.3f" %s) for s in sed_er['haze'] ]), xy=(10, 0.013) )
# 			# --- Monopole dipole plot
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
# 			plt.title(r'$^{\rm red}\chi^2=$'+', '.join(chisq)+r'$\pm0.0028$'+' vs \n'+', '.join(chisq_noh), fontsize=12)
# 			# plt.title('Dark matter model: '+model+' - Southern region')
# 			plt.xlabel(r'$\nu$ [GHz]')
# 			# plt.ylabel('Template amplitudes')
# 
# 			ax.set_xlim([10,100])
# 			ax.set_xscale('log')
# 			hd = np.reshape([sed['monopole'],sed['dipole_x'],sed['dipole_y'],sed['dipole_z'],sed_noh['monopole'],sed_noh['dipole_x'],sed_noh['dipole_y'],sed_noh['dipole_z']], 8*len(sed['monopole']) )
# 			thrs = 1.1 * np.max( abs( hd ) )
# 			ax.set_ylim([-thrs,thrs])
# 			plt.legend( loc='lower right',prop={'size':10})
# 			fig.savefig('../sed_dm_'+model+'_southern.png')
# 
# 		# --- Setting starting point
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
# 
# if plot_min_sed:
# 	pl.show()

# ====================================================================
for imodel,model in enumerate(model_list):
	data = Data( DMmodel=model )
	runtag = model+tag
	print runtag

# ------ Analyze MCMC ------
	if analyze_mcmc:
		# runtag = runtag + '_cmvt'
		file = open( '../mcmc_chain_'+runtag+'.pysav','r' )
		bf, ch = pickle.load( file )
		file.close

		get_cov = True
		if get_cov:
			cov = m.get_covmat( ch, skiplines=0.3 )
			# print cov
			file = open( '../mcmc_chain_covmat_'+runtag+'.pysav','w' )
			pickle.dump( cov, file )
			file.close
			
			w,v = m.get_eig( cov )
			print w

# 		test = False
# 		if test:
# 			c = np.array([[2,1],[1,3]])
# 			print c
# 			print ''
# 			w,v = m.get_eig( c )
# 			print w
# 			print ''
# 			print v
# 			v = -v
# 			print v[:,0]
# 			print v[:,1]
# 		
# 	#		print ''
# 	#		print np.dot(v, np.array([1,0]))
# 	#		print ''
# 	#		print np.dot(v, np.array([0,1]))
# 		
# 			vm1 = np.linalg.inv(v)
# 			print ''
# 			print vm1
# 			print np.dot(v,vm1)
# 			# pt = np.array( [2,3] )
# 			np.random.seed(seed=4)
# 			fig = plt.figure(11, (7,7))
# 			ax = plt.subplot(111)
# 			pt = np.random.randn( 5000,2 )
# 			center = [3,4]
# 			# pt[:,0] += center[0]
# 			# pt[:,1] += center[1]
# 			plt.plot(pt[:,0], pt[:,1], '.', color='k')
# 			plt.plot([0,0], [-20,20],color='k')
# 			plt.plot([-20,20], [0,0],color='k')
# 			plt.plot(center[0],center[1],'.',color='r')
# 			ax.set_xlim(-18,18)
# 			ax.set_ylim(-18,18)
# 
# 			ptr = np.zeros((5000,2),dtype=float)
# 			for ip in range(5000):	
# 				ptr[ip,:] = np.dot( np.dot(np.dot(pt[ip,:],v),np.array([[np.sqrt(w[0]),0],[0,np.sqrt(w[1])]]) ), v.T )
# 
# 			plt.plot( ptr[:,0], ptr[:,1], '.',color='b')
# 			pl.show()
# 
# 
# 			pt = np.array([[1,0],[0,1],[2,2,]])
# 			pt = np.reshape(pt,(3,2))
# 			print pt
# 			print np.dot([1,0],v)
# 			print np.dot([0,1],v)
# 			# vs = np.zeros( (2,2),dtype=float )
# 			# vs[:,0] = -v[:,0] * np.sqrt(w[0])
# 			# vs[:,1] = v[:,1] * np.sqrt(w[1])
# 
# 			# vm1 = np.linalg.inv(v)
# 			# vsm1 = np.linalg.inv(vs)
# 
# 			ptr = np.dot( pt, vm1 )
# 			# print np.dot( pt, v )
# 			print ptr
# 
# 			fig = plt.figure(12)
# 			ax = plt.subplot(111)
# 			plt.plot(pt[:,0], pt[:,1], '.', color='k')
# 			plt.plot(ptr[:,0], ptr[:,1], '.', color='b')
# 			for ip in range(pt.shape[0]):
# 				plt.plot([pt[ip,0],ptr[ip,0]], [pt[ip,1],ptr[ip,1]], '--', color='k')
# 				
# 			# plt.plot(pt[0,:], pt[1,:], '.', color='k')
# 			# ptr += vs[:,0] * np.sqrt(w[0]) * 0.25
# 			# ptr += vs[:,1] * np.sqrt(w[1]) * (-0.37)
# 			# plt.plot(ptr[:,0], ptr[:,1], '.', color='g')
# 
# 			# ptrb = np.dot(ptr, vs)
# 			# plt.plot(ptrb[:,0], ptrb[:,1], '.', color='y')
# 			plt.plot([0,0], [-20,20],color='k')
# 			plt.plot([-20,20], [0,0],color='k')
# 
# 			plt.plot([0.,v[0,0]], [0,v[1,0]],color='r')
# 			plt.plot([0.,v[0,1]], [0,v[1,1]],color='r')
# 			# plt.plot([0.,-vm1[0,0]], [0,-vm1[1,0]],color='b')
# 			# plt.plot([0.,vm1[0,1]], [0,vm1[1,1]],color='b')
# 
# 			ax.set_xlim(-8,8)
# 			ax.set_ylim(-8,8)
# 			pl.show()

		if run_2mcmc:
			bf, ch = m.sampler( p, data.fg_chisq, nsamples=200000, seed=4, output_tag=runtag+'_2nd',
				temperature=1., check_point=runtag, accept_first=False,
				covmat=cov, update_covmat=True )

		lm = m.marginalize( ch, skiplines=0.3, nbins=21 )

		# --- Analyze LIKEs ---
		if plot_likes:
			for i,fi in enumerate(data.data.freq_tags):

				cnt = 330
				fig = plt.figure(i+1)
				plt.title(fi+' likelihoods')
				cnt += 1
				ax = plt.subplot(cnt)
				plt.plot(lm[fi+'_dust'][0],lm[fi+'_dust'][1])
				ax.set_xlabel('dust')
				cnt += 1
				ax = plt.subplot(cnt)
				plt.plot(lm[fi+'_freefree'][0],lm[fi+'_freefree'][1])
				cnt += 1
				ax = plt.subplot(cnt)
				plt.plot(lm[fi+'_cmb'][0],lm[fi+'_cmb'][1])
				cnt += 1
				ax = plt.subplot(cnt)
				plt.plot(lm[fi+'_sync'][0],lm[fi+'_sync'][1])
				cnt += 1
				ax = plt.subplot(cnt)
				plt.plot(lm[fi+'_monopole'][0],lm[fi+'_monopole'][1])
				cnt += 1
				ax = plt.subplot(cnt)
				plt.plot(lm[fi+'_dipole_x'][0],lm[fi+'_dipole_x'][1])
				cnt += 1
				ax = plt.subplot(cnt)
				plt.plot(lm[fi+'_dipole_y'][0],lm[fi+'_dipole_y'][1])
				cnt += 1
				ax = plt.subplot(cnt)
				plt.plot(lm[fi+'_dipole_z'][0],lm[fi+'_dipole_z'][1])
				if make_pics:
					fig.savefig( "../like_"+runtag+"_"+fi+".png" )

			fig = plt.figure(8)
			ax = plt.subplot(111)
			plt.plot(lm['global_haze'][0],lm['global_haze'][1])
			plt.title( 'Haze amplitude for model '+model )
			if make_pics:
				fig.savefig( "../like_hazeAmpl_"+runtag+".png" )
			pl.show()
			
		# --- Analyze MAPs ---
		if plot_maps:
			for i,fi in enumerate(data.data.freq_tags):
				pi = OrderedDict()
				pi['global_haze'] = bf['parameters']['global_haze']
				for key in bf['parameters'].keys():
					if key.split('_')[0] == fi:
						# print key
						pi[key] = bf['parameters'][key]
				# print pi
				bf_map = data.fg_model( pi )
				print data.fg_chisq( pi ) / np.sum(data.mask.map)
				map = getattr(data.data, fi)
				hp.mollview( (map.imap.map - bf_map)*data.mask.map, min=-0.025, max=0.025, title='Residuals @ '+fi+' GHz', unit='mK', fig=20+1 )
				hp.graticule( dpar=20, dmer=20 )
				pl.show()
		
		# --- Analyze SEDs ---
		if plot_sed:
			cfreq = data.data.cfreq
			cv = conversionfactor( cfreq, option='antenna2thermo' )
			sed = get_sed( bf['parameters'] )
			# print sed

			fig = plt.figure(9)
			ax = plt.subplot(111)
			plt.plot( cfreq, sed['dust'], label='dust' )
			plt.plot( cfreq, sed['cmb'], label='cmb' )
			plt.plot( cfreq, sed['freefree'], label='freefree' )
			plt.plot( cfreq, sed['sync'], label='sync' )

			plt.plot( cfreq, (cfreq/cfreq[0])**(-3.1)*cv*sed['sync'][0], '--', label='-3.1', color='k' )
			plt.plot( cfreq, (cfreq/cfreq[0])**(-2.15)*cv*sed['freefree'][0], '--', label='-2.15', color='k' )

			plt.xlabel(r'$\nu$')
			plt.ylabel('Template amplitudes')

			ax.set_yscale('log')
			ax.set_xscale('log')
			ax.set_ylim([5e-2,20])
			ax.set_xlim([10,100])
			plt.legend( loc='upper left')
			if make_pics:
				fig.savefig( "../seds_"+runtag+".png" )

			pl.show()
