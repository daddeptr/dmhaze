import numpy as np
import sys
import healpy as hp
from collections import OrderedDict

from myclass import *

# --------------------------------------
class Component( MyClass ):
	def __init__(self, filename='', set_ones=False, scaling=1.):
		self.filename = filename
		if len(filename) == 0:
			if not set_ones:
				self.map = 0.
			else:
				self.map = 1. * scaling
		else:
			self.map = hp.read_map( filename ) * scaling

# --------------------------------------
class Frequency( MyClass ):
	def __init__( self ):
		self.imap = Component()
		self.irms = Component()
		self.cfreq = 0.

	def define( self, datainfo ):
		# intensity
		imap = getattr( self, 'imap' )
		file = datainfo['ifilename']
		setattr( imap, 'filename', file)
		map = hp.read_map( file )
		setattr( imap, 'map', map)
		setattr( self, 'imap', imap)

		# rms
		irms = getattr( self, 'irms' )
		file = datainfo['infilename']
		setattr( irms, 'filename', file)
		map = hp.read_map( file )
		setattr( irms, 'map', map)
		setattr( self, 'irms', irms)

		# cfreq
		setattr( self, 'cfreq', datainfo['cfreq'])

# --------------------------------------
class Galaxy( MyClass ):
	def __init__( self, galactic_templates, selection=[] ):
		if len(selection) == 0:
			selection = haze_templates.keys()
		for key in selection:
			file = galactic_templates[key]
			comp = Component( file )
			setattr( self, key, comp)

# --------------------------------------
class Haze( MyClass ):
	def __init__( self, haze_templates, selection=[] ):
		if len(selection) == 0:
			selection = haze_templates.keys()
		for key in selection:
			file = haze_templates[key]
			comp = Component( file )
			setattr( self, key, comp)
					
# --------------------------------------
class Data( MyClass ): 
# 
# 	# --------------------------------------
# 	class Component( MyClass ):
# 		def __init__(self, filename='', set_ones=False, scaling=1.):
# 			self.filename = filename
# 			if len(filename) == 0:
# 				if not set_ones:
# 					self.map = 0.
# 				else:
# 					self.map = 1. * scaling
# 			else:
# 				self.map = hp.read_map( filename ) * scaling
# 
# 	# ----------------------------------
# 	class Frequency( MyClass ):
# 		def __init__( self ):
# 			self.imap = Component()
# 			self.irms = Component()
# 			self.cfreq = 0.
# 
# 		def define( self, datainfo ):
# 			# intensity
# 			imap = getattr( self, 'imap' )
# 			file = datainfo['ifilename']
# 			setattr( imap, 'filename', file)
# 			map = hp.read_map( file )
# 			setattr( imap, 'map', map)
# 			setattr( self, 'imap', imap)
# 
# 			# rms
# 			irms = getattr( self, 'irms' )
# 			file = datainfo['infilename']
# 			setattr( irms, 'filename', file)
# 			map = hp.read_map( file )
# 			setattr( irms, 'map', map)
# 			setattr( self, 'irms', irms)
# 
# 			# cfreq
# 			setattr( self, 'cfreq', datainfo['cfreq'])
# 
# 	# ------------------------------
# 	class Galaxy( MyClass ):
# 		def __init__( self, galactic_templates, selection=[] ):
# 			if len(selection) == 0:
# 				selection = haze_templates.keys()
# 			for key in selection:
# 				file = galactic_templates[key]
# 				comp = Component( file )
# 				setattr( self, key, comp)
# 		
# 	# ------------------------------
# 	class Haze( MyClass ):
# 		def __init__( self, haze_templates, selection=[] ):
# 			if len(selection) == 0:
# 				selection = haze_templates.keys()
# 			for key in selection:
# 				file = haze_templates[key]
# 				comp = Component( file )
# 				setattr( self, key, comp)
					
	# --------------------------------------
	class FGData( MyClass ):

		# ----------------------------------
		def __init__(self, selection=[], dir='/global/scratch2/sd/dpietrob/DM-haze/data/maps/',
			nside=128, DMmodel='32232', raw_sync=False):
# 
# 			# ------------------------------
# 			class Galaxy( MyClass ):
# 				def __init__( self, galactic_templates, selection=[] ):
# 					if len(selection) == 0:
# 						selection = haze_templates.keys()
# 					for key in selection:
# 						file = galactic_templates[key]
# 						comp = Component( file )
# 						setattr( self, key, comp)
# 					
# 			# ------------------------------
# 			class Haze( MyClass ):
# 				def __init__( self, haze_templates, selection=[] ):
# 					if len(selection) == 0:
# 						selection = haze_templates.keys()
# 					for key in selection:
# 						file = haze_templates[key]
# 						comp = Component( file )
# 						setattr( self, key, comp)
# 					
			snside = str(nside).zfill(4)
			fgdir = dir+'ns'+snside+'/'

			# --- Define templates
			galactic_templates = {}
			galactic_templates['cmb']      = fgdir + 'planck_hfi_ilc_ns'+snside+'_1deg_mK.fits'
			galactic_templates['dust']     = fgdir + 'map_fds_dust_94_ns'+snside+'_1deg.fits'
			galactic_templates['freefree'] = fgdir + 'map_halpha_ns'+snside+'_1deg_norm.fits'
			if not raw_sync:
				galactic_templates['sync']     = fgdir + 'map_haslam_dsds_ns'+snside+'_1deg_norm.fits'
			else:
				galactic_templates['sync']     = fgdir + 'map_haslam_nofilt_ns'+snside+'_1deg_norm.fits'
			galactic_templates['dipole_x'] = fgdir + 'map_xdipole_ns'+snside+'.fits'
			galactic_templates['dipole_y'] = fgdir + 'map_ydipole_ns'+snside+'.fits'
			galactic_templates['dipole_z'] = fgdir + 'map_zdipole_ns'+snside+'.fits'
			galactic_templates['disc']     = fgdir + 'map_mydisc_ns'+snside+'_new.fits'
			# galactic_templates['bubbles']     = fgdir + 'fermi_bubbles_ns'+snside+'.fits'
			galactic_templates['bubbles']     = fgdir + 'map_myhaze_ns'+snside+'_new.fits'

			if ( (len(selection) == 0) or (selection[0].lower() == 'all') ):
				selection = ['cmb', 'dust', 'freefree', 'sync', 'dipole_x', 'dipole_y', 'dipole_z', 'disc', 'bubbles']

	# 		galaxy = Galaxy( selection )
	# 		galaxy.define( galactic_templates, selection )
			galaxy = Galaxy( galactic_templates, selection=selection )
		
			#print galaxy


			haze_templates = {}
	# 		haze_templates['single'] = fgdir + 'map_myhaze_ns'+snside+'.fits'
			#flist = ['fK', 'f030']
	# --- Generic template for all frequencies
			flist = ['generic']
			for item in flist:
				haze_templates[ item ]   = fgdir + 'map_myhaze_ns'+snside+'_new.fits'

	# --- Tailored template per frequency
			# DMmodel = '3225' 
			flist = ['K','030','Ka','Q','044','V','070','W']
			for item in flist:
				if DMmodel != 'myhaze':
					haze_templates[ 'f'+item ]   = fgdir + 'haze_model_54_'+DMmodel+'_ns'+snside+'_'+item+'_norm.fits'
				else:
					haze_templates[ 'f'+item ]   = fgdir + 'map_'+DMmodel+'_ns'+snside+'_new_'+item+'.fits'

	# 		haze = Haze( frequency_list='haze_f'+flist)
	# 		haze.define( haze_templates )
			haze = Haze( haze_templates )

			#print haze
		
			self.galaxy = galaxy
			self.gal_list = selection
			self.haze = haze
			self.haz_list = haze_templates.keys()

	# --------------------------------------
	class FreqData( MyClass ):
		
		# ----------------------------------
		def __init__(self, selection=['K','030','Ka','Q','044','V','070','W'],
			dir='/global/scratch2/sd/dpietrob/DM-haze/data/maps/',
			nside=128):

			frequencies = ['K','030','Ka','Q','044','V','070','W']
			# cfreq = [22.8, 28.4, 33., 41., 44., 61., 70.4, 94.]
			# cfreq = np.array( [23., 28., 33, 41., 44., 61., 70., 93.] )
			cfreq = {'K':23., '030':28., 'Ka':33, 'Q':41., '044':44., 'V':61., '070':70., 'W':93.} # Andrey's values
			snside = 'ns' + str(nside).zfill(4)
			ddir = dir + snside

			info = {}
			for freq in frequencies:
				datainfo = {}
				datainfo['ifilename'] = ddir + '/map_'+freq+'_'+snside+'_1deg_mK.fits' 
				datainfo['infilename'] = ddir + '/rms_adjst_'+freq+'_'+snside+'_1deg_mK.fits' 
				datainfo['cfreq'] = cfreq[freq]
				info[freq] = datainfo

			self.frequencies = []
			self.freq_tags = []
			self.cfreq = []

			for freq in selection:
				fcomp = Frequency()
				fcomp.define( info[freq] )
				setattr(self, 'f'+freq, fcomp )
				self.frequencies.append( freq )
				self.freq_tags.append( 'f'+freq )
				self.cfreq.append( cfreq[freq] )

			self.cfreq = np.array( self.cfreq )
			
	# --------------------------------------
	def fg_chisq(self, pars ):
		keys = pars.keys()
		# print keys
		freq = []
		for key in keys:
			if key[0] == 'f':
				freq.append( key.split('_')[0] )
		frequencies = set( freq )
		# print frequencies
		nfreq = len( frequencies )
	
		chi2 = 0.e0
	
		for freq in frequencies:
			# print freq
			channel = getattr( self.data, freq )
			map = channel.imap.map
			rms = channel.irms.map
			npix = len(map)
			nside = int( np.sqrt(npix/12) )
			model = np.zeros( npix )
			# --- Get the model
			for key in keys:
				kfreq = key.split('_')[0]
				# print kfreq, freq
				#print key[0:len(freq)], freq
				if kfreq == freq:
					compname = '_'.join( key.split('_')[1:] )
					# print compname
	#				if ( (compname.lower() != 'monopole') and (compname.lower()[0:6] != 'dipole') ):
					if ( (compname.lower() != 'monopole') ):
						comp = getattr( self.gmodel.galaxy, compname )
						model += comp.map * pars[ key ]
					elif compname.lower() == 'monopole':
						model += pars[ key ]
	#				elif compname.lower() == 'dipole_x':
	#					model += make_dipole( nside, [pars[key],0.,0.] )
	#				elif compname.lower() == 'dipole_y':
	#					model += make_dipole( nside, [0.,pars[key],0.] )
	#				elif compname.lower() == 'dipole_z':
	#					model += make_dipole( nside, [0.,0.,pars[key]] )
					else:
						info_message( 'unknown parameter: '+key+'. Abort')
						sys.exit()
				elif key[0:6] == 'haze_f':
					compname = 'generic'
					comp = getattr( self.gmodel.haze, compname )
					model += comp.map * pars[ key ]
				elif key == 'global_haze':
					comp = getattr( self.gmodel.haze, freq )
					model += comp.map * pars[ key ]
	# When ever there's a frequency mismatch it prints.
	# 			else:
	# 				info_message( key+', '+freq)
			
			# hp.mollview(map, norm='hist')
			# pl.show()
			# hp.mollview(rms, norm='hist')
			# pl.show()
			# hp.mollview(map-model, norm='hist')
			# pl.show()
			c2 = np.sum( ( (map-model)/rms)**2 * self.mask.map )
			# print freq, c2/np.sum( data.mask.map )
			chi2 += c2
			#print chi2
			#sys.exit()

		return chi2

	# --------------------------------------
	def fg_model(self, pars, nside=128l ):
		keys = pars.keys()
		for key in keys:
			if key[0] == 'f':
				freq = key.split('_')[0]
				break
		# print freq
		# print keys
		# nside = int( np.sqrt(npix/12) )
		npix = 12l * nside**2
		model = np.zeros( npix )

		for key in keys:
			if key[0] == 'f':
				compname = '_'.join( key.split('_')[1:] )
	#			if ( (compname.lower() != 'monopole') and(compname.lower()[0:6] != 'dipole') ):
				if ( (compname.lower() != 'monopole') ):
					comp = getattr( self.gmodel.galaxy, compname )
					if len(comp.map) != npix:
						info_message('Inconsistent nside: '+str(nside)+', '+len(comp.map)+'. Abort')
						sys.exit()
					model += comp.map * pars[ key ]
				elif compname.lower() == 'monopole':
					model += pars[ key ]
	#			elif compname.lower() == 'dipole_x':
	#				model += make_dipole( nside, [pars[key],0.,0.] )
	#			elif compname.lower() == 'dipole_y':
	#				model += make_dipole( nside, [0.,pars[key],0.] )
	#			elif compname.lower() == 'dipole_z':
	#				model += make_dipole( nside, [0.,0.,pars[key]] )
				else:
					info_message( 'unknown parameter: '+key+'. Abort')
					sys.exit()			
			elif key[0:6] == 'haze_f':
				compname = 'generic'
				comp = getattr( self.gmodel.haze, compname )
				model += comp.map * pars[ key ]
			elif key == 'global_haze':
				comp = getattr( self.gmodel.haze, freq )
				model += comp.map * pars[ key ]

		return model

	# --------------------------------------
	def __init__(self, frequencies=[ 'K','030','Ka','Q','044','V','070'],
		maskfile='ebvh_mask_top-half_disc_40_ns0128.fits', nside=128l,
		dir='/global/scratch2/sd/dpietrob/DM-haze/data/maps/',
		DMmodel='32232', raw_sync=False, regression=0 ):

		self.data = self.FreqData( selection=frequencies, nside=nside, dir=dir )
		self.gmodel = self.FGData( dir=dir, DMmodel=DMmodel, raw_sync=raw_sync)

		if len(maskfile) != 0:
			snside = 'ns' + str(nside).zfill(4)
			if dir[-1] == '/':
				ddir = dir + snside
			else:
				ddir = dir + '/' + snside
			
			self.mask = Component( filename=ddir+'/'+maskfile )
		else:
			self.mask = Component( set_ones=True )
		
		if regression == 1:
			for tag in self.data.freq_tags:
				freq = getattr( self.data, tag )
				comp = getattr( freq, 'imap' )				
				tmp = remove_dipole( comp.map, mask=self.mask.filename, onlymonopole=True )[0]
				setattr( comp, 'map', tmp )
				setattr( freq, 'imap', comp )
				setattr( self.data, tag, freq )
			
			for tag in self.gmodel.gal_list:
				comp = getattr( self.gmodel.galaxy, tag )
				tmp = remove_dipole( comp.map, mask=self.mask.filename, onlymonopole=True )[0]
				setattr( comp, 'map', tmp )
				setattr( self.gmodel.galaxy, tag, comp )
				
			
		if regression == 2:
			for tag in self.data.freq_tags:
				freq = getattr( self.data, tag )
				comp = getattr( freq, 'imap' )
				
				tmp = remove_dipole( comp.map, mask=self.mask.filename )[0]
				setattr( comp, 'map', tmp )
				setattr( freq, 'imap', comp )
				setattr( self.data, tag, freq )
			
			for tag in self.gmodel.gal_list:
				comp = getattr( self.gmodel.galaxy, tag )
				tmp = remove_dipole( comp.map, mask=self.mask.filename )[0]
				setattr( comp, 'map', tmp )
				setattr( self.gmodel.galaxy, tag, comp )
				
# --------------------------------------
def get_sed( pars ):
	keys = pars.keys()
	freq = []
	components = []
	for key in keys:
		tags = key.split('_')
		if tags[0][0] == 'f':
			if tags[0] not in freq:
				freq.append( tags[0] )
		if (len( tags )>1) and (tags[0] != 'global'):
			if '_'.join(tags[1:]) not in components:
				components.append( '_'.join(tags[1:]) )
	# print freq
	# print components
	
	nfreq = len(freq)
	ncomp = len(components)

	sed = {}
	for j,c in enumerate(components):
		comp = []
		for i,f in enumerate(freq):
			key = f+'_'+c
			# print key
			comp.append( pars[key] )
		sed[c] = np.array(comp)

	return sed

# --------------------------------------------------------------------
# frequencies = [ 'K' ,'030','Ka','Q','044','V','070'] #,'W' ]

# cfreq = np.array( [23., 28., 33, 41., 44., 61., 70.] ) #, 93.] )

# data = FreqData()

