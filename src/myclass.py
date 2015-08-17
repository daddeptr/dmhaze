import inspect
# import numpy as np
# import healpy as hp
# import sys

# --------------------------------------
def isclass( class_obj ):
	return inspect.isclass( class_obj )
	
# --------------------------------------
def myrepr( class_obj ):
	attrs = dir( class_obj )
	out = []
	for i,attr in enumerate(attrs):
		if attr[0:2] != "__":
			out.append(attr)
	return out

# --------------------------------------
def myattrs( class_obj ):
	tags = dir( class_obj )
	attrs = []
	for attr in tags:
		if (attr[0] != "_") and (attr != 'keys') and (attr != 'vals') \
			and (attr != 'default') and (attr != 'set_pars') and (attr != 'allowed_models'):# \
#			 and (attr != 'info'):
			attrs.append(attr)
	return attrs
	
# --------------------------------------
def printme( class_obj ):
	name = str( class_obj.__class__)
	values = name+": \n"
	#exclude = ['keys','vals','default','set_pars','allowed_models','info']
	exclude = ['keys','vals','default','set_pars','allowed_models']
	for item in myrepr( class_obj ):
		if item not in exclude:
			# assuming that class names begins w/ capital letters
			if item[0] != item[0].upper():
				try:
					attr = str( getattr(class_obj,item) )
					if (len(attr) >0) and (attr[0] != '<'):
						values += "   "+item+" = "+attr+"\n"
					else:
						values += "   "+item+" = "+attr.replace("\n","\n   ")+"\n"
				except ValueError, NameError:
					value += "   "+item+" = NOT DEFINED\n"
	return values

# --------------------------------------
class MyClass( object ):
	def __repr__(self):
		structure = "Attributes: \n"
		for item in myrepr( self ):
			structure += "   "+item+"\n"
		return structure

#	def __str__(self):
#		return str(self.__dict__)
	def __str__(self):
		return printme( self )

	def __eq__(self, other):
		if self.__class__ != other.__class__:
			return False
		else:
			return self.__dict__ == other.__dict__

# ----------------------------------
def info_message( message, code='' ):
	print code + ' >>> ' + message

# ----------------------------------
def conversionfactor( frequencies, option='cgs2mKthermo', T_cmb=2.726e0 ):
	from numpy import exp, array, ndarray, float64
	from sys import exit

	options = [ 
		'brightness2mKthermo',
		'brightness2mKantenna',
		'mKthermo2brightness',
		'mKantenna2brightness',
		'brightness2muKthermo',
		'brightness2muKantenna',
		'muKthermo2brightness',
		'muKantenna2brightness',
		'antenna2thermo',
		'thermo2antenna',
		'cgs2muKthermo',
		'cgs2mKthermo',
		'muKthermo2cgs',
		'mKthermo2cgs']

	if option not in options:
		info_message( ' conversion UNKNOWN. Abort' )
		info_message( ' Allowed conversions are:' )
		print options
		exit()
             
	if not isinstance(frequencies, ndarray):
		frequencies = array( frequencies, dtype=float64 )
# +
#  NAME:
# 	conversionfactor
# 
#  PURPOSE:
#        Delivers conversion factor between milli/mu K(antenna),
#        milli/mu K(thermodynamic) and brightness MJy/sr. Grand-unified subroutine.      
# 
#  CALLING SEQUENCE:
#  	A=conversionfactor(frequencies,/antenna2thermo)
# 
#  INPUTS:
#       FREQUENCIES REAL ARRAY    List of frequencies (GHz).
# 
#  OPTIONAL INPUTS:
# 
#  OPTIONAL INPUT KEYWORDS:
#       /brightness2mKthermo           Perform conversion MJy/sr -> mKthermo
#       /brightness2mKantenna          MJy/sr -> mKantenna
#       /mKthermo2brightness           mKthermo -> MJy/sr
#       /mKantenna2brightness          mKantenna -> MJy/sr
#       /brightness2muKthermo          Perform conversion MJy/sr -> muKthermo
#       /brightness2muKantenna         MJy/sr -> muKantenna
#       /muKthermo2brightness          muKthermo -> MJy/sr
#       /muKantenna2brightness         muKantenna -> MJy/sr
#       /antenna2thermo                Antenna -> thermodynamic
#       /thermo2antenna                Thermodynamic -> antenna
# 
#  NOTES
# 
#  SIDE EFFECTS
# 
#  EXAMPLES
#        Convert from antenna to thermodynamic units.
# 
#        correctionfactor=conversionfactor(70,/antenna2thermo)
# 
#  COMMONS USED : 
# 
#  PROCEDURES USED: 
# 
#  MODIFICATION HISTORY:
#  	May 2006, Samuel Leach, SISSA
# 


# ############Fundamental constants##############
	c = 2.99792458e10                 #  cm / s
	h = 6.626176e-27                  #  erg * s
	k = 1.380662e-16                  #  erg / K

	x = h * frequencies * 1.e9 / k / T_cmb

# Just compute ALL conversion factors and be done with it.
	conversion_thermo2antenna = x**2 * exp(x)/(exp(x)-1.)**2
	conversion_antenna2thermo = 1./conversion_thermo2antenna

	conversion_mKantenna2brightness = 1.e-3 * 2. * k * ( frequencies * 1.e9/c)**2 * 1.e17
	conversion_brightness2mKantenna = 1./conversion_mKantenna2brightness

	conversion_muKantenna2brightness = 1.e-6 * 2. * k * ( frequencies * 1.e9/c)**2 * 1.e17
	conversion_brightness2muKantenna = 1./conversion_muKantenna2brightness

	conversion_brightness2mKthermo = conversion_brightness2mKantenna * conversion_antenna2thermo
	conversion_mKthermo2brightness = 1./conversion_brightness2mKthermo

	conversion_brightness2muKthermo = conversion_brightness2muKantenna * conversion_antenna2thermo
	conversion_muKthermo2brightness = 1./conversion_brightness2muKthermo

	conversion_brightness2cgs = 1.e6 * 1.e-23
	conversion_cgs2brightness = 1. / conversion_brightness2cgs

	conversion_cgs2muKthermo = conversion_cgs2brightness * conversion_brightness2muKthermo
	conversion_muKthermo2cgs = 1. / conversion_cgs2muKthermo

	conversion_cgs2mKthermo = conversion_cgs2brightness * conversion_brightness2mKthermo
	conversion_mKthermo2cgs = 1. / conversion_cgs2mKthermo

# Default output
	if option == 'antenna2thermo':
		conversion = conversion_antenna2thermo
	elif option == 'thermo2antenna':
		conversion = conversion_thermo2antenna
	elif option == 'brightness2mKthermo':
		conversion = conversion_brightness2mKthermo
	elif option == 'brightness2mKantenna':
		conversion = conversion_brightness2mKantenna
	elif option == 'brightness2muKthermo':
		conversion = conversion_brightness2muKthermo
	elif option == 'brightness2muKantenna':
		conversion = conversion_brightness2muKantenna
	elif option == 'mKthermo2brightness':
		conversion = conversion_mKthermo2brightness
	elif option == 'mKantenna2brightness':
		conversion = conversion_mKantenna2brightness
	elif option == 'muKthermo2brightness':
		conversion = conversion_muKthermo2brightness
	elif option == 'muKantenna2brightness':
		conversion = conversion_muKantenna2brightness
	elif option == 'cgs2muKthermo':
		conversion = conversion_cgs2muKthermo
	elif option == 'cgs2mKthermo':
		conversion = conversion_cgs2mKthermo
	elif option == 'muKthermo2cgs':
		conversion = conversion_muKthermo2cgs
	elif option == 'mKthermo2cgs':
		conversion = conversion_mKthermo2cgs

	return conversion

# ----------------------------------
def xyGaussian(nside,sigx,sigy,moll=False):
	import healpy as hp
	import numpy as np
	import pylab as pl
	
	npix = 12l*nside**2
	ipix = np.array( range(npix) )

	sigx = sigx * np.pi/180.	
	sigy = sigy * np.pi/180.	

	theta, phi = hp.pixelfunc.pix2ang( nside, ipix )
	theta = np.pi/2. - theta

	snd = np.array( phi > np.pi )
	
	phi[snd] = phi[snd] - 2.*np.pi
	
	g = np.zeros( npix, dtype='float' )
	
	g = np.exp( -0.5 * ( (theta/sigy)**2 + (phi/sigx)**2 ) )

	if moll:
		hp.mollview(g, fig=1)
		pl.show()
	
	if sigx == sigy:
		vec0 = np.reshape( hp.pixelfunc.ang2vec( np.pi/2, 0. ), (3) )
		vec = np.reshape( hp.pixelfunc.pix2vec( nside, ipix ), (npix,3) )
	
		angdist = np.acos( np.dot(vec, vec0) )
	
		ga = np.exp( -0.5*(angdist/sigx)**2 )
		hp.mollview(ga, fig=2)
		pl.show()

	return g

# ----------------------------------
def defaultHaze(frequencies=[22.8,28.4,33.,41.,44.,61.,70.4, 94.], nside=128, indx=-2.55,
	sigx=15., sigy=25., labels=['K','030','Ka','Q','044','V','070','W']):
	import numpy as np
	import healpy as hp
	
	npix = 12l*nside**2
	nfreq= len(frequencies)
	if not isinstance(frequencies, np.ndarray):
		frequencies = np.array(frequencies, dtype=float)

	g = xyGaussian( nside, sigx, sigy )
	nuref = 22.8
	
	cv = conversionfactor( frequencies, option='antenna2thermo' )
	for ifreq, freq in enumerate(frequencies):
		gf = g * (freq/nuref)**indx * cv[ifreq]
		hp.write_map( '/global/scratch2/sd/dpietrob/DM-haze/data/maps/ns'+str(nside).zfill(4)+'/map_myhaze_ns'+str(nside).zfill(4)+'_new_'+labels[ifreq]+'.fits', gf)

	return True

# --------------------------------------
def make_dmhaze_templates( model='3225'):
	import pyfits as pf
	import healpy as hp
	
	dir = '/global/homes/d/dpietrob/myscratch/DM-haze/data/dmhaze/'
	frequencies = ['K', '030', 'Ka', 'Q', '044', 'V', '070', 'W']
	f = [23., 28., 33., 41., 44., 61., 70., 93.]
	cv = conversionfactor( f )
	print cv
	nfreq = len( frequencies )
	# haze = hp.read_map( dir+'synchrotron_healpix_54_'+model+'.fits', field=(range(nfreq)) )
	haze = hp.read_map( dir+'synchrotron_healpix_54_'+model, field=(range(nfreq)) )
	
	for ifreq, freq in enumerate(frequencies):
		# filename = dir + '../maps/ns0128/haze_model_54_'+model+'_ns0128_'+freq+'.fits'
		# hp.write_map( filename, haze[ifreq] )
		filename = dir + '../maps/ns0128/haze_model_54_'+model+'_ns0128_'+freq+'_norm.fits'
		hp.write_map( filename, haze[ifreq] * cv[ifreq] )
		# hp.mollview( haze[ifreq] )
		# pl.show()

# --------------------------------------
def make_dipole( nside, direction ):
	import numpy as np
	import healpy as hp

	npix = 12l * nside**2
	dipole = np.zeros( npix, dtype=float )
	ipix = np.array( range(npix) )
#	vecs = np.reshape( hp.pixelfunc.pix2vec( nside, ipix, nest=True ), (npix,3) )
	vecs = hp.pixelfunc.pix2vec( nside, ipix )
	direction = np.array( direction )
#	dipole = np.dot( vecs, direction )
	dipole = vecs[0] * direction[0] + vecs[1] * direction[1] + vecs[2] * direction[2]
#	print dipole.shape
#	hp.mollview( dipole )
#	pl.show()
#	sys.exit()
	return dipole
	
# --------------------------------------
def make_sky_cut( nside, angle_deg ):

	import numpy as np
	import healpy as hp
	
	info_message( 'nside     = '+str(nside) )
	info_message( 'angle_deg = '+str(angle_deg) )

	npix = 12l*nside**2

	mask = np.ones(npix, dtype=float)
	ipix = np.array(range(npix))

	theta, phi = hp.pixelfunc.pix2ang( nside, ipix )

	r2d = 180. / np.pi
	theta_deg = 90. - theta * r2d

	bp = np.array( abs(theta_deg) < angle_deg )

	mask[bp] = 0.	

	return mask

# --------------------------------------
def make_disc( nside, radius, direction=[1,0,0], radians=False, minus=False):

	import numpy as np
	import healpy as hp
	
	info_message( 'nside  = ' + str(nside) )
	info_message( 'radius = ' + str(radius) )
	npix = 12l*nside**2

	mask = np.zeros( npix )
	d2r = np.pi / 180.
	if not radians:
		radius *= d2r

	listpix = hp.query_disc( nside, direction, radius )

	mask[listpix] = 1.

	if minus:
		mask = 1. - mask

	return mask

# --------------------------------------
def map_regression( inmap, template_files, label='Residuals', return_fit=False, maskfile=None,
	do_plot=False, regression=0, rms=1., cut=0., scale=1.,
	debug=False, show=False, return_res=False, return_chi2=False, verbose=True ):

	import numpy as np
	import healpy as hp

	if maskfile != None:
		if verbose:
			info_message( '- Using mask:'+maskfile )

	if debug:   
		print type( inmap )

	if isinstance( inmap, str ):
		if verbose:
			info_message('reading map from "'+inmap+'"')
		inmap = hp.read_map( infile )
	elif not isinstance( inmap, np.ndarray ):
		info_message( ' - ERROR: input failure. Abort.')
		sys.exit()

	
	templates_are_files = False
	if not isinstance(template_files, str):
		if not isinstance(template_files, np.ndarray):
			template_files = list( template_files )
		if isinstance( template_files[0], str):
			templates_are_files = True

	npix = len( inmap )
	nside = int( np.sqrt( float(npix) / 12 ) )

	if maskfile != None:
		mask = hp.read_map( maskfile )
		if len(mask) != npix:
			hp.pixel_func.ud_grade( mask, nside )
	else:
		mask = np.ones( npix, dtype=float )

	if show and verbose:
		hp.mollview( mask, norm='hist' )
		pl.show()

	if cut != 0:
		c = make_sky_cut(nside, cut)
		mask *= c

	bp = np.array( mask == 0.)
	gp = np.array( mask == 1.)
	ngp = np.sum( gp )
	nbp = np.sum( bp )

	if  regression == 1:
		if verbose:
			info_message( 'removing monopole' )
		inmap = remove_dipole( inmap, onlymonopole=True, cut=cut, mask=mask )[0]
	if  regression == 2:
		inmap = remove_dipole( inmap, cut=cut, mask=mask )[0]
		if verbose:
			info_message( 'removing monopole and dipole' )

	if not isinstance(rms, np.ndarray):
		if isinstance(rms, str):
			rms = hp.read_map( rms )
		else:
			rms = np.ones(npix, dtype=np.float64)

	if show and verbose:
		hp.mollview( rms, norm='hist' )
		pl.show()

	if show and verbose:
		cc = inmap
		cc[bp] = hp.UNSEEN
		hp.mollview( cc, min=-300, max=300, title='!6Input map')
		pl.show()

	nfg = len( template_files )
	info_message( 'nfg = ' + str(nfg) )
	# A = np.zeros( nfg, dtype=np.float64 )
	A_err = np.zeros( nfg, dtype=np.float64 )
	oneOsigma = np.zeros( (nfg, nfg), dtype=np.float64 )
	# sigma = np.zeros( (nfg, nfg), dtype=np.float64 )

	temp = np.zeros( (npix,nfg), dtype=np.float64 )
	for ifg, template in enumerate(template_files):
		if templates_are_files:
			t = hp.read_map( template )
		else:
			t = template
			template_files[ifg] = 'template_'+str(ifg+1).zfill(2)
		if regression == 1:
			t = remove_dipole( t, onlymonopole=True, cut=cut, mask=mask )[0]
		if regression == 2:
			t = remove_dipole( t, cut=cut, mask=mask )[0]
		temp[:,ifg] = t
		if (show or debug) and verbose:
			hp.mollview( temp[:,ifg], norm='asinh' ) 
			pl.show()

	T = np.zeros( (nfg, nfg), dtype=np.float64 )
	B = np.zeros( nfg, dtype=np.float64 )

	warn = True
	for ifg in range(nfg):
		"""
       if keyword_set(do_corr) then begin
           if warn then begin
               print, ' - Correlate used!'
               warn = False
           endif
###       B[ifg] = total( (map[gp]-mean(map[gp]))*(temp[gp,ifg]-mean(temp[gp,ifg])) / rms[gp]^2 ) / ngp
           B[ifg] = correlate( map[gp]/rms[gp], temp[gp,ifg]/rms[gp], /covariance, /double )
           for jfg=ifg,nfg-1 do begin
###           T[ifg,jfg] = total( (temp[gp,ifg]-mean(temp[gp,ifg]))*(temp[gp,jfg]-mean(temp[gp,jfg])) / rms[gp]^2 ) / ngp
               T[ifg,jfg] = correlate( temp[gp,ifg]/rms[gp], temp[gp,jfg]/rms[gp], /covariance, /double )
               T[jfg,ifg] = T[ifg,jfg]
               oneOsigma[ifg,jfg] = 2. * T[ifg,jfg]
               oneOsigma[jfg,ifg] = oneOsigma[ifg,jfg]
           endfor
       endif else begin
		"""
# --- Changed
		if warn:
			info_message( ' - Code changed: correlate is not used anymore!' )
			warn = False

		B[ifg] = np.sum( inmap[gp]/rms[gp] * temp[gp,ifg]/rms[gp] )
		for jfg in range(ifg,nfg):
			T[ifg,jfg] = np.sum( temp[gp,ifg] * temp[gp,jfg] / rms[gp]**2 )
			T[jfg,ifg] = T[ifg,jfg]

			oneOsigma[ifg,jfg] = 2. * T[ifg,jfg]
			oneOsigma[jfg,ifg] = oneOsigma[ifg,jfg]

	if debug:
		print T, B

	Tm1 = np.linalg.inv( T )
       
	A = np.dot( Tm1, B )
       
	sigma = np.linalg.inv( oneOsigma )

	for kfg in range( nfg ):
		if (sigma[kfg,kfg] > 0.):
			A_err[kfg] = np.sqrt( sigma[kfg,kfg] )
		else:
			info_message( 'Sigma Error' )

	if verbose:
		info_message( ' ============================================================' )
		info_message( 'Coefficients:' )
		for ifg in range( nfg ):
			print template_files[ifg], A[ifg],' +/-',A_err[ifg]
		info_message( ' ============================================================' )
       
	fit = np.zeros( npix, dtype=np.float64 )
	for ifg in range( nfg ):
		fit += temp[:,ifg] * A[ifg]

	if (show or debug) and verbose:
		cc = fit
		cc[bp] = hp.UNSEEN
		thrs = np.min(abs(np.array([np.min(fit),np.max(fit)])))
		hp.mollview( cc, title='Template Linear Combination', min=-thrs, max=thrs )
		pl.show()

	res = inmap - fit
	if (show or debug) and verbose:
		cc = res
		cc[bp] = -1.6375e30
		thrs = np.min(abs(np.array([np.min(fit),np.max(fit)])))
		hp.mollview( cc, title='Residual Foregrounds', min=-thrs, max=thrs )

	# info_message( ' - total chi2 = ' + str( np.sum( (res/rms)**2 ) / npix ) + ', ' + str( 1./np.sqrt(npix) ) )
	chi2 = np.sum( (res[gp]/rms[gp])**2 ) / ngp
	if verbose:
		info_message( ' - total chi2 = ' + str( chi2 ) + ', ' + str( 1./np.sqrt(ngp) ) )
	# info_message( ' - total chi2 = ' + str( np.sum( (res[bp]/rms[bp])**2 ) / nbp ) + ', ' + str( 1./np.sqrt(nbp) ) )

	if verbose:
		info_message( ' --- End of Procedure ---' )
	
	outcome = [A, A_err]
	if return_fit:
		outcome.append( fit )
	if return_res:
		outcome.append( res )
	if return_chi2:
		outcome.append( chi2 )


	return tuple(outcome)

# --------------------------------------
def remove_dipole( inmap, return_fit=False, mask=None, rms=1., cut=0.,
	debug=False, show=False, return_A=False, verbose=True, onlymonopole=False,
	onlydipole=False ):

	import numpy as np
	import healpy as hp

	if debug:
		print type( inmap )

	if isinstance( inmap, str ):
		if verbose:
			info_message('reading map from "'+inmap+'"')
		inmap = hp.read_map( infile )
	elif not isinstance( inmap, np.ndarray ):
		info_message( ' - ERROR: input failure. Abort.')
		sys.exit()

	npix = len( inmap )
	nside = int( np.sqrt( float(npix) / 12 ) )

	# --- setting up mask
	if mask != None:
		if isinstance(mask,str):
			if verbose:
				info_message( '- Using mask:'+mask )
   
			mask = hp.read_map( mask )
		elif not isinstance( mask, np.ndarray):
			info_message( 'mask type not understood. Abort.' )
			sys.exit()
			
		if len(mask) != npix:
			hp.pixel_func.ud_grade( mask, nside )
	else:
		mask = np.ones( npix, dtype=float )

	if show and verbose:
		hp.mollview( mask, norm='hist' )
		pl.show()

	if cut != 0:
		c = make_sky_cut(nside, cut)
		mask *= c

	bp = np.array( mask == 0.)
	gp = np.array( mask == 1.)
	ngp = np.sum( gp )
	nbp = np.sum( bp )

	if not isinstance(rms, np.ndarray):
		if isinstance(rms, str):
			rms = hp.read_map( rms )
		else:
			rms = np.ones(npix, dtype=np.float64)

	if show and verbose:
		hp.mollview( rms, norm='hist' )
		pl.show()

	if show and verbose:
		cc = inmap
		cc[bp] = hp.UNSEEN
		hp.mollview( cc, min=-300, max=300, title='!6Input map')
		pl.show()

	if onlymonopole:
		nfg = 1
	elif onlydipole:
		nfg = 3
	else:
		nfg = 4

	if nfg == 1:
		fit = np.mean( inmap[gp] )
		info_message( 'monopole = '+str('%f' %fit) )
	if True:
		info_message( 'nfg = ' + str(nfg) )
		A_err = np.zeros( nfg, dtype=np.float64 )
		oneOsigma = np.zeros( (nfg, nfg), dtype=np.float64 )

		temp = np.zeros( (npix,nfg), dtype=np.float64 )
		if onlymonopole:
			template_files = ['monopole']
			temp[:,0] = np.ones( npix )
		elif onlydipole:
			template_files = ['dipole_x','dipole_y','dipole_z']
			temp[:,0] = make_dipole( nside, [1.,0.,0.] )
			temp[:,1] = make_dipole( nside, [0.,1.,0.] )
			temp[:,2] = make_dipole( nside, [0.,0.,1.] )
		else:
			template_files = ['monopole','dipole_x','dipole_y','dipole_z']
			temp[:,0] = 1.
			temp[:,1] = make_dipole( nside, [1.,0.,0.] )
			temp[:,2] = make_dipole( nside, [0.,1.,0.] )
			temp[:,3] = make_dipole( nside, [0.,0.,1.] )

		if (show or debug) and verbose:
			for ifg in range(nfg):
				hp.mollview( temp[:,ifg], norm='asinh' ) 
			pl.show()

		T = np.zeros( (nfg, nfg), dtype=np.float64 )
		B = np.zeros( nfg, dtype=np.float64 )

		warn = True
		for ifg in range(nfg):
# --- Changed
			if warn:
				info_message( ' - Code changed: correlate is not used anymore!' )
				warn = False

			B[ifg] = np.sum( inmap[gp]/rms[gp] * temp[gp,ifg]/rms[gp] )
			for jfg in range(ifg,nfg):
				T[ifg,jfg] = np.sum( temp[gp,ifg] * temp[gp,jfg] / rms[gp]**2 )
				T[jfg,ifg] = T[ifg,jfg]

				oneOsigma[ifg,jfg] = 2. * T[ifg,jfg]
				oneOsigma[jfg,ifg] = oneOsigma[ifg,jfg]

		if debug:
			print T, B

		Tm1 = np.linalg.inv( T )
       
		A = np.dot( Tm1, B )
       
		sigma = np.linalg.inv( oneOsigma )

		for kfg in range( nfg ):
			if (sigma[kfg,kfg] > 0.):
				A_err[kfg] = np.sqrt( sigma[kfg,kfg] )
			else:
				info_message( 'Sigma Error' )

		if verbose:
			info_message( ' ============================================================' )
			info_message( 'Coefficients:' )
			for ifg in range( nfg ):
				print template_files[ifg], A[ifg],' +/-',A_err[ifg]
			if onlydipole:
				vec = np.array( A )
			elif not onlymonopole:
				vec = np.array( A[1:4] )
				theta, phi = hp.vec2ang( vec )
				theta_deg = 90.-theta*180./np.pi
				phi_deg = phi*180./np.pi
				dip_ampl = np.sqrt( np.dot(vec,vec) )
				info_message( '(theta,phi) = ('+str('%f' %theta_deg)+','+str('%f' %phi_deg)+')' )
				info_message( 'A_dip       = '+str('%f' %dip_ampl) )
			info_message( ' ============================================================' )
       
		fit = np.zeros( npix, dtype=np.float64 )
		for ifg in range( nfg ):
			fit += temp[:,ifg] * A[ifg]

		if (show or debug) and verbose:
			cc = fit
			cc[bp] = hp.UNSEEN
			thrs = np.min(abs(np.array([np.min(fit),np.max(fit)])))
			hp.mollview( cc, title='Template Linear Combination', min=-thrs, max=thrs )
			pl.show()

	res = np.array(inmap - fit)
	if (show or debug) and verbose:
		cc = res
		cc[bp] = -1.6375e30
		thrs = np.min(abs(np.array([np.min(fit),np.max(fit)])))
		hp.mollview( cc, title='Residual Foregrounds', min=-thrs, max=thrs )

	# info_message( ' - total chi2 = ' + str( np.sum( (res/rms)**2 ) / npix ) + ', ' + str( 1./np.sqrt(npix) ) )
	# chi2 = np.sum( (res[gp]/rms[gp])**2 ) / ngp
	# if verbose:
	# 	info_message( ' - total chi2 = ' + str( chi2 ) + ', ' + str( 1./np.sqrt(ngp) ) )
	# info_message( ' - total chi2 = ' + str( np.sum( (res[bp]/rms[bp])**2 ) / nbp ) + ', ' + str( 1./np.sqrt(nbp) ) )

	if verbose:
		info_message( ' --- End of Procedure ---' )
	
	outcome = [ res ]
	if return_fit:
		outcome.append( fit )
	if return_A:
		outcome.append( (A, A_err) )


	return tuple(outcome)

# --------------------------------------------------------------------
def asinh10( x ):
	import numpy as np
	return np.log10( 0.5e0*( x+np.sqrt(x**2+4.e0)) )

# # --------------------------------------------------------------------
# def read_map( filename, field=(0,1,2) ):
# 	#import healpy as hp
# 	try:
# 		m = hp.read_map( filename, field=field )
# 	except IndexError:
# 		print " psm_fitsio.read_map >>> IndexError on '%s', try reading 1 field" %filename
# 		m = hp.read_map( filename )
# 
# 	if isinstance(m,np.ndarray):
# 		return m
# 	elif isinstance(m,tuple):
# 		c = len(m)
# 		r = len(m[0])
# 		a = np.zeros( (r,c), dtype=np.float64 )
# 		for i in range(c):
# 			a[:,i] = m[i]
# 		return a
# 		
# # --------------------------------------------------------------------
# def write_map( filename, map, ordering='RING' ):
# 	#import healpy as hp
# 	code = 'psm_hpx.write_map'
# 	nd = map.ndim
# 	sh = map.shape
# 	print code+' >>> writing: '+filename
# 	#print code+' >>> ', nd, sh
# 	if nd == 1:
# 		if ordering[0:4].upper() == 'RING':
# 			hp.write_map( filename, map )
# 		elif ordering[0:4].upper() == 'NEST':
# 			hp.write_map( filename, map, nest=True )
# 		else:
# 			info_message( 'psm_hpx.write_map >>> "ordering" not understood: '+ordering+'. Abort' )
# 			sys.exit()	
# 	elif (nd == 2) and (sh[1] == 3):
# 		if ordering[0:4].upper() == 'RING':
# 			hp.write_map( filename, (map[:,0], map[:,1], map[:,2]) )
# 		elif ordering[0:4].upper() == 'NEST':
# 			hp.write_map( filename, (map[:,0], map[:,1], map[:,2]), nest=True )
# 		else:
# 			info_message( 'psm_hpx.write_map >>> "ordering" not understood: '+ordering+'. Abort' )
# 			sys.exit()	
# 	else:
# 		print " psm_hpx.write_map >>> Arbitrary number of maps not implemented (1/3). Abort: '%s', %2d" %(filename,sh)
# 		sys.exit()
# 
