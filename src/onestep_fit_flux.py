import matplotlib.pyplot as plt
import pylab as pl
import pickle
import healpy as hp
import sys

from myclass import *
from mcmc import *
from data_model import *

def onestep_fit_flux(model_list=['myhaze'], mfile='ebvh_mask_ns0128.fits', tag='galmask',
	fix_cmb=False, inc_disc=True, inc_mono=True, inc_dipo=True, inc_haze=True, inc_bubbles=True,
	plot_gnom=False,
	plot_moll=True, find_minimum=True, regression=0, generic_haze=False, plot_mondipo=True, 
	plot_min_sed=True, write_png=True, raw_sync=True, frequencies=[ 'K','030','Ka','Q','044','V','070']):

# 	selection = ['3225','3225h','3235','5114','12232','32232','myhaze']
# 	print model_list
# 	for model in model_list:
# 		if model not in selection:
# 			info_message('ERROR: unknown model. Abort')
# 			print model_list
# 			sys.exit()


	dir = '/global/scratch2/sd/dpietrob/DM-haze/data/maps/ns0128/'

	# Mask bright object n the field
	if False:
		m = hp.read_map( dir+'south_haze_mask_ns0128.fits' )

		x,y,z = hp.ang2vec(np.pi/2+21.5/180.*np.pi ,0.75/180.*np.pi)
		ptsrc = make_disc(128, 2, direction=[x,y,z], minus=True)
		print np.sum(m*ptsrc)/len(m)

		x,y,z = hp.ang2vec(np.pi/2+19.5/180.*np.pi ,359.5/180.*np.pi)
		ptsrc1 = make_disc(128, 2, direction=[x,y,z], minus=True)
		print np.sum(m*ptsrc*ptsrc1)/len(m)

		hp.write_map( dir + 'south_haze_ptsrc_mask_ns0128.fits', m*ptsrc*ptsrc1 )
		hp.gnomview( m*ptsrc*ptsrc1, rot=[0,-22.5], reso=25, min=-1, max=1)
		hp.graticule(dpar=5,dmer=5)

	# --- Making a more suitable mask ---
	# 30 deg radius from [0,-20]
	radius = -1 #30
	if False:
		m = hp.read_map( dir + 'south_haze_ptsrc_mask_ns0128.fits' )
		# hp.mollview(m)
		print np.sum(m)/len(m)

		x,y,z = hp.ang2vec(np.pi/2+22.5/180.*np.pi ,0.)
		d = make_disc(128, radius, direction=[x,y,z])
		# hp.mollview(d)
		# hp.mollview(d*m)
		print np.sum(m*d)/len(m)

		hp.write_map( dir + 'south'+str(radius)+'_haze_ptsrc_mask_ns0128.fits', d*m )
		m = hp.read_map( dir + 'south'+str(radius)+'_haze_ptsrc_mask_ns0128.fits' )
		hp.gnomview( m, rot=[0,-22.5], reso=25, min=-1, max=1)
		hp.graticule(dpar=5,dmer=5)
		# hp.gnomview( ptsrc, rot=[0,-22.5], reso=25, min=-1, max=1)
		# hp.graticule(dpar=5,dmer=5)

		m = hp.read_map( dir + 'fit_mask_p-35_35_t-10-35_ns0128.fits' )
		# hp.mollview(m)
		print np.sum(m)/len(m)

		pl.show()
	# sys.exit()

	if False:
		if radius>0:
			tag = 'southern'+str(radius)	
			mfile = 'south'+str(radius)+'_haze_ptsrc_mask_ns0128.fits'
		else:
			tag = 'southern'
			mfile = 'south_haze_ptsrc_mask_ns0128.fits'

	mask = hp.read_map( dir + mfile )
	mask_npix = np.sum( mask )
	fsky = np.sum( mask ) / len(mask)
	gpix = np.array( (mask == 1.) )
	bpix = np.array( (mask == 0.) )

	m = hp.read_map( dir + mfile )

	if False:
		x,y,z = hp.ang2vec(np.pi/2+22.5/180.*np.pi ,0.)
		d = make_disc(128, 20, direction=[x,y,z])
	else:
		# haze_fit_region = 'fit_mask_p-35_35_t-10-35_ns0128.fits'
		haze_fit_region = 'fit_mask_p-20_20_t-10-35_ns0128.fits'
		d = hp.read_map( dir + haze_fit_region )
	haze_region = (m*d == 1)
	not_haze_region = (m*d == 0)

	if not inc_haze:
		tag = tag + '_nohaze'

	if regression >= 1:
		inc_mono = False
		tag = tag + '_nomono'
		if regression == 2:
			inc_dipo = False
			tag = tag + 'dipo'

	# ------ Run MCMC ------
	for imodel, model in enumerate( model_list) :
		print model
		# data = Data( DMmodel=model, maskfile='ebvh_mask_ns0128.fits', raw_sync=True )
		data = Data( DMmodel=model, maskfile=mfile, raw_sync=raw_sync, frequencies=frequencies )
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
			haze_flux_er_scaled = []
			# haze_templ_flux = []
			# haze_templ_flux_er = []
			bubbles_flux = []
			bubbles_flux_er = []

			# --- Find minimum for each frequency
			for ifreq, freq in enumerate( data.data.frequencies ):

				if generic_haze:
					h = getattr( gmodel.haze, 'generic' )
				else:
					h = getattr( gmodel.haze, 'f'+freq )

				# --- Full run for starting point determination
				# Thin disc component removed: chisq gets worse, but global fit does
				# not change drmatically.
				# tfiles = [ gmodel.galaxy.cmb.filename, gmodel.galaxy.dust.filename, gmodel.galaxy.freefree.filename, 
				tfiles = [ gmodel.galaxy.dust.filename,
					gmodel.galaxy.freefree.filename, 
					gmodel.galaxy.sync.filename]
				idust = 0
				iff   = 1
				isync = 2
				# ihaze = 3
				if inc_haze:
					tfiles.append( h.filename )
					ihaze = len(tfiles)-1

				if inc_mono:
					tfiles.append( dir+'map_monopole_ns0128.fits' )
					imono = len(tfiles)-1
				if inc_dipo:
					tfiles.append( dir+'map_xdipole_ns0128.fits' )
					idipox = len(tfiles)-1

					tfiles.append( dir+'map_ydipole_ns0128.fits' )
					idipoy = len(tfiles)-1

					tfiles.append( dir+'map_zdipole_ns0128.fits' )
					idipoz = len(tfiles)-1
				if inc_disc:
					tfiles.append( gmodel.galaxy.disc.filename )
					idisc = len(tfiles)-1
				if inc_bubbles:
					tfiles.append( gmodel.galaxy.bubbles.filename )
					ibubbles = len(tfiles)-1
				if not fix_cmb:
					tfiles.append( gmodel.galaxy.cmb.filename )
					icmb = len(tfiles)-1


				map = getattr( data.data, 'f'+freq )
				regmap = map.imap.map
				if fix_cmb:
					regmap -= gmodel.galaxy.cmb.map
					
				r = map_regression( regmap, tfiles, rms=map.irms.map,
					maskfile=data.mask.filename, regression=regression, 
					return_res=True, return_chi2=True )
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

				if not fix_cmb:
					fx_mnpnt['f'+freq+'_cmb'] = r[0][icmb]
					fx_mnpnt_er['f'+freq+'_cmb'] = r[1][icmb]
				fx_mnpnt['f'+freq+'_dust'] = r[0][idust]
				fx_mnpnt_er['f'+freq+'_dust'] = r[1][idust]

				fx_mnpnt['f'+freq+'_freefree'] = r[0][iff]
				fx_mnpnt_er['f'+freq+'_freefree'] = r[1][iff]

				fx_mnpnt['f'+freq+'_sync'] = r[0][isync]
				fx_mnpnt_er['f'+freq+'_sync'] = r[1][isync]

				if inc_haze:
					fx_mnpnt['f'+freq+'_haze'] = r[0][ihaze]
					fx_mnpnt_er['f'+freq+'_haze'] = r[1][ihaze]
				if inc_disc:
					fx_mnpnt['f'+freq+'_disc'] = r[0][idisc]
					fx_mnpnt_er['f'+freq+'_disc'] = r[1][idisc]
				if inc_bubbles:
					fx_mnpnt['f'+freq+'_bubbles'] = r[0][ibubbles]
					fx_mnpnt_er['f'+freq+'_bubbles'] = r[1][ibubbles]
				if inc_mono:
					fx_mnpnt['f'+freq+'_monopole'] = r[0][imono]
					fx_mnpnt_er['f'+freq+'_monopole'] = r[1][imono]
				if inc_dipo:
					fx_mnpnt['f'+freq+'_dipole_x'] = r[0][idipox]
					fx_mnpnt_er['f'+freq+'_dipole_x'] = r[1][idipox]

					fx_mnpnt['f'+freq+'_dipole_y'] = r[0][idipoy]
					fx_mnpnt_er['f'+freq+'_dipole_y'] = r[1][idipoy]

					fx_mnpnt['f'+freq+'_dipole_z'] = r[0][idipoz]
					fx_mnpnt_er['f'+freq+'_dipole_z'] = r[1][idipoz]

				# gnomview projection
				if plot_gnom:
					ccc = map.imap.map
					ccc[ not_haze_region] = -1.6375e30
					hp.gnomview( ccc, min=-0.3, title=freq,	rot=[0,-22.5], reso=25 )
					hp.graticule(dpar=5,dmer=5)

					if ifreq == 0 and False:
						if fix_cmb:
							ccc = data.gmodel.galaxy.cmb.map
						else:
							ccc = data.gmodel.galaxy.cmb.map * r[0][icmb]
						ccc[ not_haze_region ] = -1.6375e30
						hp.gnomview( ccc, min=-0.3, max=0.3, title='CMB', rot=[0,-22.5], reso=25 )
						hp.graticule(dpar=5,dmer=5)

						ccc = data.gmodel.galaxy.sync.map
						ccc[ not_haze_region ] = -1.6375e30
						hp.gnomview( ccc, title='Sync', rot=[0,-22.5], reso=25 )
						hp.graticule(dpar=5,dmer=5)

						ccc = data.gmodel.galaxy.dust.map
						ccc[ not_haze_region ] = -1.6375e30
						hp.gnomview( ccc, title='Dust', rot=[0,-22.5], reso=25 )
						hp.graticule(dpar=5,dmer=5)

						ccc = data.gmodel.galaxy.freefree.map
						ccc[ not_haze_region ] = -1.6375e30
						hp.gnomview( ccc, title='FF', rot=[0,-22.5], reso=25 )
						hp.graticule(dpar=5,dmer=5)

					ccc = r[2]*mask
					ccc[ not_haze_region ] = -1.6375e30
					hp.gnomview( ccc, min=-0.05, max=0.05, title='Residuals: '+freq,
						rot=[0,-22.5], reso=25 )
					hp.graticule(dpar=5,dmer=5)

					if inc_haze:
						ccc = (r[2]+r[0][ihaze]*h.map)*mask
						ccc[ not_haze_region ] = -1.6375e30
						hp.gnomview( ccc, title='Haze: '+freq, rot=[0,-22.5], reso=25 )
						hp.graticule(dpar=5,dmer=5)

					# ccc = (r[0][2]*h.map)*mask
					# ccc[ not_haze_region ] = -1.6375e30
					# hp.gnomview( ccc, min=0., title='Haze template: '+freq, rot=[0,-22.5], reso=25 )
					# hp.graticule(dpar=5,dmer=5)
			
				# mollview projection
				if plot_moll:
# 					ccc = regmap
# 					ccc[ bpix] = -1.6375e30
# 					hp.mollview( ccc, min=-0.3, title=freq,	rot=[0,-22.5], unit='$\mu K$')
# 					hp.graticule(dpar=20,dmer=20)

					if ifreq == 0 and False:
						if fix_cmb:
							ccc = data.gmodel.galaxy.cmb.map
						else:
							ccc = data.gmodel.galaxy.cmb.map * r[0][icmb]
						ccc[ bpix] = -1.6375e30
						hp.mollview( ccc, min=-0.3, max=0.3, title='CMB', rot=[0,-22.5] )
						hp.graticule(dpar=20,dmer=20)

						ccc = data.gmodel.galaxy.sync.map
						ccc[ bpix] = -1.6375e30
						hp.mollview( ccc, title='Sync', rot=[0,-22.5])
						hp.graticule(dpar=20,dmer=20)

						ccc = data.gmodel.galaxy.dust.map
						ccc[ bpix] = -1.6375e30
						hp.mollview( ccc, title='Dust', rot=[0,-22.5])
						hp.graticule(dpar=20,dmer=20)

						ccc = data.gmodel.galaxy.freefree.map
						ccc[ bpix] = -1.6375e30
						hp.mollview( ccc, title='FF', rot=[0,-22.5])
						hp.graticule(dpar=20,dmer=20)

					ccc = r[2]*mask
					ccc[ bpix] = -1.6375e30
					hp.mollview( ccc, min=-0.05, max=0.05, title='Residuals: '+freq,
						rot=[0,-22.5], unit='$\mu K$' )
					hp.graticule(dpar=20,dmer=20)

					if inc_haze:
						ccc = (r[2]+r[0][ihaze]*h.map)*mask
						ccc[ bpix] = -1.6375e30
						hp.mollview( ccc, title='Haze: '+freq, rot=[0,-22.5], unit='$\mu K$', min=-0.05, max=np.max(ccc))
						hp.graticule(dpar=20,dmer=20)

					# ccc = (r[0][2]*h.map)*mask
					# ccc[ bpix] = -1.6375e30
					# hp.gnomview( ccc, min=0., title='Haze template: '+freq, rot=[0,-22.5], reso=25 )
					# hp.graticule(dpar=5,dmer=5)
			
				# --- Flux with residuals
	# 			haze_flux.append( np.mean( (r[2] + r[0][ihaze]*h.map )[gpix] ) * mKt2cgs[ifreq] )
	# 			# --- Trying to be independent of the haze template: residuals only
	# 			haze_flux_er.append( np.std( r[2][gpix] ) * mKt2cgs[ifreq] )
	# 			haze_flux_er_scaled.append( np.sqrt( np.var( r[2][gpix] ) + (haze_flux[-1]*r[1][ihaze])**2)* mKt2cgs[ifreq] )
				# --- Trying haze selection
				if inc_haze:
					haze_flux.append( np.mean( (r[2] + r[0][ihaze]*h.map )[haze_region] ) * mKt2cgs[ifreq] )
					# --- Trying to be independent of the haze template: residuals only
					haze_flux_er.append( np.std( r[2][haze_region] ) * mKt2cgs[ifreq] )
					haze_flux_er_scaled.append( np.sqrt( np.var( r[2][haze_region] ) + (haze_flux[-1]*r[1][ihaze])**2)* mKt2cgs[ifreq] )

				freq_min[ freq ] = r

			if inc_haze:
				flux[model+'_haze'] = ( np.array( haze_flux )*1.e20,
					np.array( haze_flux_er )*1.e20,
					np.array( haze_flux_er_scaled )*1.e20 )

				if write_png:
					figtag = '../sed_dm_' + model + '_1step_' + tag

					file = open(figtag + '.txt', 'w')
					file.write( '# cfreq, res_haze, res_haze err, templ_haze, templ_haze err\n' )
					for ifreq,freq in enumerate(data.data.cfreq):
					# print ifreq, freq
						file.write( '%f \t %f \t %f \t %f \n' %(freq, flux[model+'_haze'][0][ifreq],
							flux[model+'_haze'][1][ifreq], flux[model+'_haze'][2][ifreq] ) )
					file.close()
			# print flux

			info_message(' total chi2 = '+str(c2)+','+str(c2*np.sum(mask)))

			fx_sed = get_sed( fx_mnpnt )
			fx_sed_er = get_sed( fx_mnpnt_er )
	
			if plot_min_sed:
				f = data.data.cfreq
				cv = conversionfactor( f, option='antenna2thermo' )
				fig = plt.figure(30+imodel,(16,4.5), dpi=80)

				# --- SED plot based on regression coefficients
				ax = plt.subplot(131)
	
				if not fix_cmb:
					plt.errorbar( f, fx_sed['cmb'], yerr=fx_sed_er['cmb'], label='cmb', color='k', fmt='.-' )
				plt.errorbar( f, fx_sed['dust'], yerr=fx_sed_er['dust'], label='dust', color='b', fmt='.-' )
				plt.errorbar( f, fx_sed['sync'], yerr=fx_sed_er['sync'], label='sync', color='g', fmt='.-' )
				if inc_haze:
					plt.errorbar( f, fx_sed['haze'], yerr=fx_sed_er['haze'], label='haze', color='m', fmt='.-' )
				plt.errorbar( f, fx_sed['freefree'], yerr=fx_sed_er['freefree'], label='freefree', color='r', fmt='.-' )
				if inc_disc:
					plt.errorbar( f, fx_sed['disc'], yerr=fx_sed_er['disc'], label='disc', color='y', fmt='.-' )
					plt.plot( f, (f/f[0])**(-3.1)*cv*fx_sed['disc'][0], ':', label='-3.1', color='y' )
				if inc_bubbles:
					plt.errorbar( f, fx_sed['bubbles'], yerr=fx_sed_er['bubbles'], label='bubbles', color='c', fmt='.-' )
					plt.plot( f, (f/f[0])**(-2.56)*cv*fx_sed['bubbles'][0], ':', label='-2.56', color='c' )
					plt.plot( f, (f/f[0])**(-3.1)*cv*fx_sed['bubbles'][0], ':', label='-3.1', color='c' )

				plt.plot( f, (f/f[0])**(-3.1)*cv*fx_sed['sync'][0], ':', label='-3.1', color='g' )
				plt.plot( f, (f/f[0])**(-2.15)*cv*fx_sed['freefree'][0], ':', label='-2.15', color='r' )
	#			plt.plot( f, (f/f[0])**(-2.55)*cv*sed['haze'][0], ':', label='-2.55', color='m' )
				# plt.plot( f, f*0.+fx_sed['haze'][0], ':', color='m', label='same amplitude' )
				# plt.plot( f, (f/f[0])**(-.1)*cv*sed['haze'][0], ':', label='-2.65', color='m' )

				# plt.title('Dark matter model: '+model)
				# plt.title('Dark matter model: '+model+' - Southern region')
				plt.xlabel(r'$\nu$ [GHz]')
				plt.ylabel('Template amplitudes')

				ax.set_yscale('log')
				ax.set_xscale('log')
				ax.set_ylim([1e-3,1000])
				ax.set_xlim([8,100])
				plt.legend( loc='upper left',prop={'size':10})

				# --- Monopole dipole plot
				if regression == 0:
					ax = plt.subplot(132)
		
					plt.errorbar( f, fx_sed['monopole'], yerr=fx_sed_er['monopole'], label='monopole', color='k', fmt='.-' )
					plt.errorbar( f, fx_sed['dipole_x'], yerr=fx_sed_er['dipole_x'], label='dipole_x', color='b', fmt='.-' )
					plt.errorbar( f, fx_sed['dipole_y'], yerr=fx_sed_er['dipole_y'], label='dipole_y', color='r', fmt='.-' )
					plt.errorbar( f, fx_sed['dipole_z'], yerr=fx_sed_er['dipole_z'], label='dipole_z', color='m', fmt='.-' )
					if inc_disc:
						plt.errorbar( f, fx_sed['disc'], yerr=fx_sed_er['disc'], label='disc', color='y', fmt='.-' )
					plt.plot( f, f*0., '--', color='g' )
	# 	 			plt.errorbar( f, sed_noh['monopole'], yerr=sed_er_noh['monopole'], label='monopole no-haze', color='k', fmt='--' )
	# 	 			plt.errorbar( f, sed_noh['dipole_x'], yerr=sed_er_noh['dipole_x'], label='dipole_x no-haze', color='b', fmt='--' )
	# 	 			plt.errorbar( f, sed_noh['dipole_y'], yerr=sed_er_noh['dipole_y'], label='dipole_y no-haze', color='r', fmt='--' )
	# 	 			plt.errorbar( f, sed_noh['dipole_z'], yerr=sed_er_noh['dipole_z'], label='dipole_z no-haze', color='m', fmt='--' )
	# 	 
	# 	 			plt.errorbar( f, sed_nod['monopole'], yerr=sed_er_nod['monopole'], label='monopole no-disc', color='k', fmt='s-' )
	# 	 			plt.errorbar( f, sed_nod['dipole_x'], yerr=sed_er_nod['dipole_x'], label='dipole_x no-disc', color='b', fmt='s-' )
	# 	 			plt.errorbar( f, sed_nod['dipole_y'], yerr=sed_er_nod['dipole_y'], label='dipole_y no-disc', color='r', fmt='s-' )
	# 	 			plt.errorbar( f, sed_nod['dipole_z'], yerr=sed_er_nod['dipole_z'], label='dipole_z no-disc', color='m', fmt='s-' )
	 
	#	 			plt.title(r'$^{\rm red}\chi^2=$'+', '.join(chisq)+r'$\pm0.0028$'+' vs \n'+', '.join(chisq_noh)+' vs \n'+', '.join(chisq_nod), fontsize=12)
					# plt.title('Dark matter model: '+model+' - Southern region')
					plt.xlabel(r'$\nu$ [GHz]')
					plt.ylabel(r'Mono/dipole amplitudes [mK_CMB]')
	 
					ax.set_xlim([10,100])
					ax.set_ylim([-0.1,0.1])
					ax.set_xscale('log')
					hd = np.reshape([fx_sed['monopole'],fx_sed['dipole_x'],fx_sed['dipole_y'],fx_sed['dipole_z']], 4*len(fx_sed['monopole']) )
					thrs = 1.1 * np.max( abs( hd ) )
					ax.set_ylim([-thrs,thrs])
					plt.legend( loc='lower right',prop={'size':10})

				if inc_haze:
					# --- SED plot based on average flux in the region.
					# Mind, this is region dependent!!!
					ax = plt.subplot(133)
					plt.errorbar( f, flux[model+'_haze'][0], yerr=flux[model+'_haze'][1], label='Haze flux', color='m', fmt='s-' )
					plt.plot( f, (f/f[0])**(-0.55)*flux[model+'_haze'][0][0], ':', color='m', label=r'$\beta$=-0.55')
					plt.plot( f, (f/f[0])**(-1.1)*flux[model+'_haze'][0][0], '--', color='m', label=r'$\beta$=-1.1')
					plt.plot( f, f*0., '--', color='k')
					# plt.errorbar( f, flux[model+'_haze'][2], yerr=flux[model+'_haze'][3], label='Haze template flux', color='m', fmt='.-' )
					# plt.errorbar( f, sed['haze'], yerr=sed_er['haze'], label='haze sed', color='m', fmt='s:' )
					# plt.errorbar( f, fx_sed['haze'], yerr=fx_sed_er['haze'], label='haze FS', color='m', fmt='d:' )
					# ax.set_yscale('log')
					# ax.set_xscale('log')
					ax.set_ylim([-.75,2.])
					ax.set_xlim([8,100])
					plt.xlabel(r'$\nu$ [GHz]')
					plt.ylabel(r'Haze flux [erg/(cm$^2$ s Hz sr)] $\times 10^{20}$')
					plt.legend( loc='upper left',prop={'size':10})
					ax.annotate(', '.join(data.data.frequencies), xy=(10, -0.5) )
					ax.annotate(', '.join([ str("%0.3f" %s) for s in flux[model+'_haze'][0] ]), xy=(10, -0.6) )
					ax.annotate(r'$\pm$'+', '.join([ str("%0.3f" %s) for s in flux[model+'_haze'][1] ]), xy=(10, -0.7) )

				figtag = 'pics/sed_dm_'+model+'_1step_'+tag
				if write_png:
					fig.savefig(figtag + '.png')

				print r' \hline'
				print r' {\bf Component} & {' + r'} & {\bf '.join(data.data.frequencies) + r'} \\'
				print r' \hline'
				print r' \hline'
				for key in fx_sed.keys():
					print r'' + key.replace('_',r'\_') + ' & '+ ' & '.join( [str('%f') %s for s in fx_sed[key]] ) + r'\\'
					print r' error & '+ ' & '.join( [str('%f') %s for s in fx_sed_er[key]] ) + r'\\'
					print r' \hline'
				print r' \hline'
					
				for ifreq, freq in enumerate( data.data.frequencies ):
					print freq + r' & ' + chisq[ifreq] + r'\\'
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

	return fx_sed