import matplotlib.pyplot as plt
import pylab as pl
import pickle
import healpy as hp
import sys
import getopt
import os

from myclass import *
from mcmc import *
from data_model import *

# --------------------------------------------------------------------
def twostep_fit_flux( model=None,
	dir = '/global/scratch2/sd/dpietrob/DM-haze/data/maps/ns0128/',
	mfile = 'south_haze_ptsrc_mask_ns0128.fits', # mask used for the region fit
	sig_scaling = 8,
	haze_fit_mask = 'fit_mask_p-20_20_t-10-35_ns0128.fits',
	fix_cmb   = True,
	fix_ff    = False, # Does/does not assign freefree parameter
	inc_disc  = False,
	## inc_bubbles  = True,
	## :DP May 2015 to get Col.4 of Tab.3 running with no bubbles
	inc_bubbles  = False,
	inc_haze  = False,
	inc_mono  = True,
	inc_dipo  = True,
	plot_gnom = False,
	seed=0,
	regression = 0,
	tag = '',
	fsmaskf = 'ebvh_mask_ns0128.fits',
	raw_sync = False
	 ):
##	fsmaskf = 'ebvh_mask_ns0128.fits' # the one used for the MCMC; the one above matches the Tab.3 onestep
##	fsmaskf = 'ebvh_mask_top-half_disc_40_ns0128.fits'

	if model == None:
		info_message( 'ERROR: dm_model not provided. Stop' )
		sys.exit()

	tag = tag + '_' + model + '_2stp_'+str(sig_scaling).zfill(2)+'x'
	tag = tag + '_' + str(seed)
	print tag

	print model

# --- Templates converted to mK_CMB
	fname = dir + 'haze_model_54_'+model+'_ns0128_K_norm.fits'
	if not os.path.isfile(fname):
		info_message( ' >>> making DM model templates...' )
	  	make_dmhaze_templates( model=model )
#
#
#
# --- Making a more suitable mask ---
# 30 deg radius from [0,-22.5]
# radius = -1
# if False:
# 	m = hp.read_map( dir + 'south_haze_ptsrc_mask_ns0128.fits' )
# 	# hp.mollview(m)
# 	print np.sum(m)/len(m)
# 
# 	x,y,z = hp.ang2vec(np.pi/2+22.5/180.*np.pi ,0.)
# 	d = make_disc(128, radius, direction=[x,y,z])
# 	# hp.mollview(d)
# 	# hp.mollview(d*m)
# 	print np.sum(m*d)/len(m)
# 
# 	hp.write_map( dir + 'south'+str(radius)+'_haze_ptsrc_mask_ns0128.fits', d*m*ptsrc*ptsrc1 )
# 	hp.gnomview( d*m*ptsrc*ptsrc1, rot=[0,-22.5], reso=25, min=-1, max=1)
# 	hp.graticule(dpar=5,dmer=5)
# 	# hp.gnomview( ptsrc, rot=[0,-22.5], reso=25, min=-1, max=1)
# 	# hp.graticule(dpar=5,dmer=5)
# 
# 	m = hp.read_map( dir + 'fit_mask_p-35_35_t-10-35_ns0128.fits' )
# 	# hp.mollview(m)
# 	print np.sum(m)/len(m)
# 
# 	pl.show()
# # sys.exit()
# 
# if radius>0:
# 	tag = 'southern'+str(radius)	
# 	mfile = 'south'+str(radius)+'_haze_ptsrc_mask_ns0128.fits'
# else:
# 	tag = 'southern'
# 	mfile = 'south_haze_ptsrc_mask_ns0128.fits'
# 
# # tag = tag + '_2stp'
# # tag = tag + '_2stp_10x'
# sig_scaling = 10
# 
# tag = tag + '_2stp_6freq_fxFF_'+str(sig_scaling).zfill(2)+'x'

	mask = hp.read_map( dir + mfile )
	mask_npix = np.sum( mask )
	fsky = np.sum( mask ) / len(mask)
	gpix = np.array( (mask == 1.) )
	bpix = np.array( (mask == 0.) )

	# --- Define haze flux computation region
	d = hp.read_map( dir + haze_fit_mask )
	haze_region = (mask*d == 1)
	not_haze_region = (mask*d == 0)

	# --- What to do 
	find_minimum = True
	run_mcmc     = True
	# ------
	plot_min_sed = False
	do_write     = False

	# ------ Run MCMC ------
	# setup  = Data( DMmodel=model, maskfile='ebvh_mask_ns0128.fits', raw_sync=True )
	# setup  = Data( DMmodel=model, maskfile=mfile, raw_sync=True )
	setup  = Data( DMmodel=model, maskfile=mfile, raw_sync=raw_sync )
	data   = setup.data
	gmodel = setup.gmodel

	freq_min = OrderedDict()
	# --- Set parameter values from full sky for some components ---
	if find_minimum:
		info_message( 'Running regression' )

		fx_mnpnt = OrderedDict()
		fx_mnpnt_er = OrderedDict()

		mnpnt = OrderedDict()
		mnpnt_er = OrderedDict()

		chisq = []

		fx_c2 = 0.
		c2 = 0.

		flux = {}
		cfreq = data.cfreq
		mKt2cgs = conversionfactor( cfreq, option='mKthermo2cgs' )
			
		haze_flux = []
		haze_flux_er = []
		haze_templ_flux = []
		haze_templ_flux_er = []

		# --- Find minimum for each frequency
		for ifreq, freq in enumerate( data.frequencies ):

			h = getattr( gmodel.haze, 'f'+freq )

			# --- Full run for starting point determination
			# Thin disc component removed: chisq gets worse, but global fit does
			# not change drmatically.
			# tfiles = [ gmodel.galaxy.cmb.filename, gmodel.galaxy.dust.filename, gmodel.galaxy.freefree.filename, 
			tfiles = [ gmodel.galaxy.cmb.filename, gmodel.galaxy.dust.filename,
				gmodel.galaxy.freefree.filename, 
				gmodel.galaxy.sync.filename ]
			icmb = 0
			idust = 1
			iff = 2
			isync = 3
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

			map = getattr( data, 'f'+freq )
			r = map_regression( map.imap.map, tfiles, rms=map.irms.map,
				maskfile=dir+fsmaskf, return_res=True, return_chi2=True )
				# maskfile=dir+'map_monopole_ns0128.fits', return_res=True, return_chi2=True )
			chisq.append( str("%.3f" % r[3]) )
			fx_c2 += r[3]
			
			# hp.mollview( r[2], title='FS '+freq, min=-0.3, max=0.3)
			
			fx_mnpnt['f'+freq+'_cmb'] = r[0][icmb]
			fx_mnpnt['f'+freq+'_dust'] = r[0][idust]
			fx_mnpnt['f'+freq+'_freefree'] = r[0][iff]
			fx_mnpnt['f'+freq+'_sync'] = r[0][isync]
			if inc_haze:
				fx_mnpnt['f'+freq+'_haze'] = r[0][ihaze]
			if inc_mono:
				fx_mnpnt['f'+freq+'_monopole'] = r[0][imono]
			if inc_dipo:
				fx_mnpnt['f'+freq+'_dipole_x'] = r[0][idipox]
				fx_mnpnt['f'+freq+'_dipole_y'] = r[0][idipoy]
				fx_mnpnt['f'+freq+'_dipole_z'] = r[0][idipoz]
			if inc_disc:
				fx_mnpnt['f'+freq+'_disc'] = r[0][idisc]
			if inc_bubbles:
				fx_mnpnt['f'+freq+'_bubbles'] = r[0][ibubbles]

			fx_mnpnt_er['f'+freq+'_cmb'] = r[1][icmb]
			fx_mnpnt_er['f'+freq+'_dust'] = r[1][idust]
			fx_mnpnt_er['f'+freq+'_freefree'] = r[1][iff]
			fx_mnpnt_er['f'+freq+'_sync'] = r[1][isync]
			if inc_haze:
				fx_mnpnt_er['f'+freq+'_haze'] = r[1][ihaze]
			if inc_mono:
				fx_mnpnt_er['f'+freq+'_monopole'] = r[1][imono]
			if inc_dipo:
				fx_mnpnt_er['f'+freq+'_dipole_x'] = r[1][idipox]
				fx_mnpnt_er['f'+freq+'_dipole_z'] = r[1][idipoz]
				fx_mnpnt_er['f'+freq+'_dipole_y'] = r[1][idipoy]
			if inc_disc:
				fx_mnpnt_er['f'+freq+'_disc'] = r[1][idisc]
			if inc_bubbles:
				fx_mnpnt_er['f'+freq+'_bubbles'] = r[1][ibubbles]

			# --- Remove monopole, dipole and CMB from the map
			if inc_mono:
				flat_comp = r[0][imono]
			if inc_dipo:
				flat_comp += gmodel.galaxy.dipole_x.map * r[0][idipox] + \
				gmodel.galaxy.dipole_y.map * r[0][idipoy] + \
				gmodel.galaxy.dipole_z.map * r[0][idipoz]
			if fix_cmb:
				flat_comp += gmodel.galaxy.cmb.map * r[0][icmb]
				# gmodel.galaxy.freefree.map * r[0][2]

			# --- Run on small patch after flat component removed
			tfiles = [ gmodel.galaxy.dust.filename,
				gmodel.galaxy.sync.filename,
				gmodel.galaxy.freefree.filename ]

			if inc_haze:
				tfiles.append( h.filename )
				ihaze = len( tfiles ) -1
			if inc_disc:
				tfiles.append( gmodel.galaxy.disc.filename )
				idisc = len( tfiles ) -1
				
			if inc_bubbles:
				tfiles.append( gmodel.galaxy.bubbles.filename )
				ibubbles = len( tfiles ) -1
				
			if not fix_cmb:
				tfiles.append( gmodel.galaxy.cmb.filename )
				icmb = len( tfiles ) -1

			#if inc_mono:
			#	tfiles.append( dir+'map_monopole_ns0128.fits' )
			#	imono = len( tfiles ) -1

			# map = getattr( data, 'f'+freq )
			# --- Setting up the smaller mask in the southern region
			imap = getattr( map, 'imap' )
			setattr(imap, 'map', map.imap.map-flat_comp )
			setattr(map, 'imap', imap )
			
			r = map_regression( map.imap.map, tfiles, rms=map.irms.map,
				maskfile=dir+mfile, return_res=True, return_chi2=True )
			chisq.append( str("%.3f" % r[3]) )
			c2 += r[3]

			# --- Plot gnomview projections
			if plot_gnom:
				ccc = map.imap.map
				ccc[ bpix] = -1.6375e30
				hp.gnomview( ccc, title=freq, rot=[0,-22.5], reso=25 )
				hp.graticule(dpar=5,dmer=5)

				if freq == '070':
					ccc = gmodel.galaxy.cmb.map
					ccc[ bpix] = -1.6375e30
					hp.gnomview( ccc, min=-0.3, max=0.3, title='CMB', rot=[0,-22.5], reso=25 )
					hp.graticule(dpar=5,dmer=5)

				if freq == 'K':
					ccc = gmodel.galaxy.sync.map
					ccc[ bpix] = -1.6375e30
					hp.gnomview( ccc, title='Sync', rot=[0,-22.5], reso=25 )
					hp.graticule(dpar=5,dmer=5)

					ccc = gmodel.galaxy.dust.map
					ccc[ bpix] = -1.6375e30
					hp.gnomview( ccc, title='Dust', rot=[0,-22.5], reso=25 )
					hp.graticule(dpar=5,dmer=5)

					ccc = gmodel.galaxy.freefree.map
					ccc[ bpix] = -1.6375e30
					hp.gnomview( ccc, title='FF', rot=[0,-22.5], reso=25 )
					hp.graticule(dpar=5,dmer=5)

				ccc = r[2]*mask
				ccc[ bpix] = -1.6375e30
				hp.gnomview( ccc, min=-0.05, max=0.05, title='Residuals: '+freq,
					rot=[0,-22.5], reso=25 )
				hp.graticule(dpar=5,dmer=5)

				ccc = (r[2]+r[0][2]*h.map)*mask
				ccc[ bpix] = -1.6375e30
				hp.gnomview( ccc, title='Haze: '+freq, rot=[0,-22.5], reso=25 )
				hp.graticule(dpar=5,dmer=5)

				# ccc = (r[0][2]*h.map)*mask
				# ccc[ bpix] = -1.6375e30
				# hp.gnomview( ccc, min=0., title='Haze template: '+freq, rot=[0,-22.5], reso=25 )
				# hp.graticule(dpar=5,dmer=5)
				pl.show()
			
			# --- Flux with residuals
			# haze_flux.append( np.mean( (r[2] + r[0][2]*h.map )[gpix] ) * mKt2cgs[ifreq] )
			haze_flux.append( np.mean( (r[2] + r[0][2]*h.map )[haze_region] ) * mKt2cgs[ifreq] )
			# haze_flux_er.append( np.std( (r[2] + r[0][2]*h.map)[gpix] ) * mKt2cgs[ifreq] )
			# --- Trying to be independent of the haze template
			# haze_flux_er.append( np.std( r[2][gpix] ) * mKt2cgs[ifreq] )
			haze_flux_er.append( np.sqrt( np.var( r[2][haze_region] ) + (haze_flux[-1]*r[1][2])**2)* mKt2cgs[ifreq] )

			freq_min[ freq ] = r

			mnpnt['f'+freq+'_dust'] = r[0][0]
			mnpnt['f'+freq+'_sync'] = r[0][1]
			mnpnt['f'+freq+'_freefree'] = r[0][2]
			if inc_haze:
				mnpnt['f'+freq+'_haze'] = r[0][ihaze]
			if inc_disc:
				mnpnt['f'+freq+'_disc'] = r[0][idisc]
			if inc_bubbles:
				mnpnt['f'+freq+'_bubbles'] = r[0][ibubbles]
			if not fix_cmb:
				mnpnt['f'+freq+'_cmb'] = r[0][icmb]
			#if inc_mono:
			#	mnpnt['f'+freq+'_monopole'] = r[0][imono]

			mnpnt_er['f'+freq+'_dust'] = r[1][0]
			mnpnt_er['f'+freq+'_sync'] = r[1][1]
			mnpnt_er['f'+freq+'_freefree'] = r[1][2]
			if inc_haze:
				mnpnt_er['f'+freq+'_haze'] = r[1][ihaze]
			if inc_disc:
				mnpnt_er['f'+freq+'_disc'] = r[1][idisc]
			if inc_bubbles:
				mnpnt_er['f'+freq+'_bubbles'] = r[1][ibubbles]
			if not fix_cmb:
				mnpnt_er['f'+freq+'_cmb'] = r[1][icmb]
			#if inc_mono:
			#	mnpnt_er['f'+freq+'_monopole'] = r[1][imono]

		flux[model+'_haze'] = ( np.array( haze_flux )*1.e20,
			np.array( haze_flux_er )*1.e20 )

		if do_write:
			figtag = '../sed_dm_' + model + '_2step_' + tag
			if fix_cmb:
				figtag = figtag + '_fixCMB'
			if not inc_disc:
				figtag = figtag + '_noDISC'
			if inc_mono:
				figtag = figtag + '_MONO'

			file = open(figtag + '.txt', 'w')
			file.write( '# cfreq, res_haze, res_haze err\n' )
			for ifreq,freq in enumerate(data.cfreq):
				# print ifreq, freq
				file.write( '%f \t %f \t %f \n' %(freq, flux[model+'_haze'][0][ifreq],
					flux[model+'_haze'][1][ifreq] ) )
			file.close()
		# print flux

		info_message(' total chi2 = '+str(c2)+','+str(c2*np.sum(mask)))

		sed = get_sed( mnpnt )
		sed_er = get_sed( mnpnt_er )

		fx_sed = get_sed( fx_mnpnt )
		fx_sed_er = get_sed( fx_mnpnt_er )
	
		if plot_min_sed:
			f = data.cfreq
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

			plt.plot( f, (f/f[0])**(-3.1)*cv*fx_sed['sync'][0], ':', label='-3.1', color='g' )
			plt.plot( f, (f/f[0])**(-2.15)*cv*fx_sed['freefree'][0], ':', label='-2.15', color='r' )
			plt.plot( f, f*0.+sed['haze'][0], ':', color='m', label='same amplitude' )

			plt.title('Dark matter model: '+model)
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
			plt.errorbar( f, flux[model+'_haze'][0], yerr=flux[model+'_haze'][1], label='Haze flux', color='m', fmt='s-' )
			plt.plot( f, f*0., '--', color='k')

			ax.set_ylim([-.75,2])
			ax.set_xlim([8,100])
			plt.xlabel(r'$\nu$ [GHz]')
			plt.ylabel(r'Haze flux [erg/(cm$^2$ s Hz sr)] $\times 10^{20}$')
			plt.legend( loc='upper left',prop={'size':10})
			ax.annotate(', '.join(data.frequencies), xy=(10, -0.5) )
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

			pl.show()

		print r'$^{\rm red}\chi^2=$'+', '.join(chisq)
		print "Total"
		print fx_c2, c2
		print "Reduced"
		print fx_c2/len(data.frequencies), c2/len(data.frequencies)

		return


# ------------

#hp.mollview(hp.read_map('/global/scratch2/sd/dpietrob/DM-haze/data/maps/ns0128/south_haze_ptsrc_mask_ns0128.fits'))
#pl.show()


