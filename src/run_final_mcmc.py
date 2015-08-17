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

# !!! Possible bug  ihaze=4... No, it's correct to set monopole and dipole
# --------------------------------------------------------------------
def run_final_mcmc( model=None,
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
	inc_mono  = True,
	inc_dipo  = True,
	plot_gnom = False,
	seed=0,
	regression = 0,
	tag = '',
	nsamples = 750000
	 ):

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
	setup  = Data( DMmodel=model, maskfile=mfile, raw_sync=True )
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
				gmodel.galaxy.sync.filename, h.filename ]
			icmb = 0
			idust = 1
			iff = 2
			isync = 3
			ihaze = 4

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
				maskfile=dir+'ebvh_mask_ns0128.fits', return_res=True, return_chi2=True )
				# maskfile=dir+'map_monopole_ns0128.fits', return_res=True, return_chi2=True )
			chisq.append( str("%.3f" % r[3]) )
			fx_c2 += r[3]
			
			# hp.mollview( r[2], title='FS '+freq, min=-0.3, max=0.3)
			
			fx_mnpnt['f'+freq+'_cmb'] = r[0][icmb]
			fx_mnpnt['f'+freq+'_dust'] = r[0][idust]
			fx_mnpnt['f'+freq+'_freefree'] = r[0][iff]
			fx_mnpnt['f'+freq+'_sync'] = r[0][isync]
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
				h.filename,
				gmodel.galaxy.freefree.filename ]

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
			mnpnt['f'+freq+'_haze'] = r[0][2]
			mnpnt['f'+freq+'_freefree'] = r[0][3]
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
			mnpnt_er['f'+freq+'_haze'] = r[1][2]
			mnpnt_er['f'+freq+'_freefree'] = r[1][3]
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

		# --- Setting starting point
		p = OrderedDict()

		limits = {}
		limits['dust'] = (0.,11.)
		limits['freefree'] = (0.,50)
		limits['sync'] = (0.,50)
		limits['cmb'] = (0.97,1.03)
		# limits['monopole'] = (-0.06,0.06)
		# limits['disc'] = (0,100)
		# limits['monopole'] = (-0.06,0.06)
		limits['disc'] = (0.,1)
		limits['bubbles'] = (0.,10.) #This is actually the Haze, not DM

		# sig_scaling = 6.
		# p['global_haze'] = {'center_value':np.min( sed['haze'] ), 'sigma':np.min(sed_er['haze'])/sig_scaling, 'min':0., 'max':1000}
		p['global_haze'] = {'center_value':sed['haze'][0], 'sigma':np.min(sed_er['haze'])/sig_scaling, 'min':0., 'max':1000.}

		for ifreq, freq in enumerate( data.frequencies ):
			m = getattr( data, 'f'+freq )
			# hp.mollview( m.imap.map )
		
			r = freq_min[ freq ]
			if not fix_cmb:
				p['f'+freq+'_cmb'] = {'center_value':np.max( np.array([0.,sed['cmb'][ifreq]])), 'sigma':np.min(sed_er['cmb'])/sig_scaling, 'min':limits['cmb'][0], 'max':limits['cmb'][1]}
			p['f'+freq+'_dust'] = {'center_value':np.max( np.array([0.,sed['dust'][ifreq]])), 'sigma':np.min(sed_er['dust'])/sig_scaling, 'min':limits['dust'][0], 'max':limits['dust'][1]}
			if not fix_ff:
				p['f'+freq+'_freefree'] = {'center_value':np.max( np.array([0.,sed['freefree'][ifreq]])), 'sigma':np.min(sed_er['freefree'])/sig_scaling, 'min':limits['freefree'][0], 'max':limits['freefree'][1]}
			p['f'+freq+'_sync'] = {'center_value':np.max( np.array([0.,sed['sync'][ifreq]])), 'sigma':np.min(sed_er['sync'])/sig_scaling, 'min':limits['sync'][0], 'max':limits['sync'][1]}
			# p['f'+freq+'_monopole'] = {'center_value':r[0][5], 'sigma':np.min(sed_er['monopole'])/sig_scaling, 'min':limits['monopole'][0], 'max':limits['monopole'][1]}
			# p['f'+freq+'_dipole_x'] = {'center_value':r[0][6], 'sigma':np.min(sed_er['dipole_x'])/sig_scaling, 'min':limits['monopole'][0], 'max':limits['monopole'][1]}
			# p['f'+freq+'_dipole_y'] = {'center_value':r[0][7], 'sigma':np.min(sed_er['dipole_y'])/sig_scaling, 'min':limits['monopole'][0], 'max':limits['monopole'][1]}
			# p['f'+freq+'_dipole_z'] = {'center_value':r[0][8], 'sigma':np.min(sed_er['dipole_z'])/sig_scaling, 'min':limits['monopole'][0], 'max':limits['monopole'][1]}
			if inc_disc:
				p['f'+freq+'_disc'] = {'center_value':np.max( np.array([0.,r[0][9]])), 'sigma':np.min(sed_er['disc'])/sig_scaling, 'min':limits['disc'][0], 'max':limits['disc'][1]}
			if inc_bubbles:
				p['f'+freq+'_bubbles'] = {'center_value':np.max( np.array([0.,sed['bubbles'][ifreq]])), 'sigma':np.min(sed_er['bubbles'])/sig_scaling, 'min':limits['bubbles'][0], 'max':limits['bubbles'][1]}

		ppp = {}
		for key in p.keys():
			print key, '\t', p[key]['center_value'], '\t', p[key]['sigma']
			ppp[key] = p[key]['center_value']

		info_message( 'Starting point chi2: '+str(setup.fg_chisq(ppp) ) )

	# sys.exit()
	# ------ Run MCMC ------
	mcmc = MCMC()
	if run_mcmc:
		bf, ch = mcmc.sampler( p, setup.fg_chisq, nsamples=nsamples, seed=seed, output_tag=tag,
			temperature=1., accept_first=False, update_covmat=False )

# --------------------------------------------------------------------
def main(argv):
	code = "DM-run_final_mcmc"

	paramfile = {}
	try:
		opts, args = getopt.getopt(argv,"m:s:o:hv",["dm_model=","seed=","outtag=","help","verbose"])
	except getopt.GetoptError:
		info_message( "unknown argument.\nCalling sequence\nrun_final_mcmc.py -m <dark_matter_model> -s <seed> -o <outtag> [-vh]", code=code )
		sys.exit(2)

	for opt, arg in opts:
		if opt in ('-h', '--help'):
			info_message( "unknown argument.\nCalling sequence\nrun_final_mcmc.py -m <dark_matter_model> -s <seed> -o <outtag> [-vh]", code=code )
			sys.exit()
		elif opt in ('-v', '--verbose'):
			Verbose = True
			if Verbose:
				info_message( "Verbose turned on.", code=code )
		elif opt in ("-o", "--outtag"):
			paramfile['o'] = arg
		elif opt in ("-m", "--dm_model"):
			paramfile['m'] = arg
		elif opt in ("-s", "--seed"):
			paramfile['s'] = int(arg)

	if ('m' not in paramfile.keys()) or ('s' not in paramfile.keys()):
		info_message( "Input parameter error. Stop.", code=code)
		if ('m' not in paramfile.keys()):
			info_message( "darm matter model missing", code=code)
		if ('s' not in paramfile.keys()):
			info_message( "mcmc seed missing", code=code)
		sys.exit(2)
	if ('o' not in paramfile.keys()):
		info_message( "output tag missing", code=code)
		paramfile['o'] = 'mcmc_chain'

	print paramfile

	run_final_mcmc( model=paramfile['m'], seed=paramfile['s'], tag=paramfile['o'] )

# --------------------------------------------------------------------
if __name__ == "__main__":
	main(sys.argv[1:])
	
	# ------ Analyze MCMC ------
# 	if analyze_mcmc:
# 		# runtag = runtag + '_cmvt'
# 		file = open( '../mcmc_chain_'+runtag+'.pysav','r' )
# 		bf, ch = pickle.load( file )
# 		file.close
# 
# 		info_message( 'Best fit chi2: '+str(bf['chisq'] ) )
# 
# 		if get_cov:
# 			cov = mcmc.get_covmat( ch, skiplines=0.3 )
# 			# print cov
# 			file = open( '../mcmc_chain_covmat_'+runtag+'.pysav','w' )
# 			pickle.dump( cov, file )
# 			file.close
# 
# 		if run_2mcmc:
# 			bf, ch = mcmc.sampler( p, fg_chisq, nsamples=200000, seed=4, output_tag=runtag+'_2nd',
# 				temperature=1., check_point=runtag, accept_first=False,
# 				covmat=cov, update_covmat=True )
# 
# 		lm = mcmc.marginalize( ch, skiplines=0.3, nbins=21 )
# 
# 		if plot_chisq:
# 			c2 = []
# 			for l in ch:
# 				# print l
# 				c2.append(l['chi2'])
# 			c2 = np.array(c2, dtype=float)
# 			fig = plt.figure(40, (13,7), dpi=60)
# 			ax = plt.subplot(121)
# 			plt.plot(c2)
# 			# ax.set_yscale('log')
# 			ax.set_xscale('log')	
# 			ax = plt.subplot(122)
# 			plt.plot(c2[2500:])
# 			# ax.set_yscale('log')
# 			# ax.set_xscale('log')	
# 			pl.show()
# 
# 		if plot_likes:
# 			for i,fi in enumerate(data.freq_tags):
# 
# 				cnt = 130
# 				fig = plt.figure(i+1, (13,5), dpi=60)
# 				plt.title(fi+' likelihoods')
# 				cnt += 1
# 				ax = plt.subplot(cnt)
# 				plt.plot(lm[fi+'_dust'][0],lm[fi+'_dust'][1])
# 				ax.set_xlabel('dust')
# 				cnt += 1
# 				if not fix_ff:
# 					ax = plt.subplot(cnt)
# 					plt.plot(lm[fi+'_freefree'][0],lm[fi+'_freefree'][1])
# 					cnt += 1
# 				ax = plt.subplot(cnt)
# 				plt.plot(lm[fi+'_sync'][0],lm[fi+'_sync'][1])
# 				cnt += 1
# 				if make_pics:
# 					fig.savefig( "../like_"+runtag+"_"+fi+".png" )
# 
# 			fig = plt.figure(8)
# 			ax = plt.subplot(111)
# 			plt.plot(lm['global_haze'][0],lm['global_haze'][1])
# 			plt.title( 'Haze amplitude for model '+model )
# 			if make_pics:
# 				fig.savefig( "../like_hazeAmpl_"+runtag+".png" )
# 
# 		if plot_maps:
# 			for i,fi in enumerate(data.freq_tags):
# 				pi = OrderedDict()
# 				pi['global_haze'] = bf['parameters']['global_haze']
# 				for key in bf['parameters'].keys():
# 					if key.split('_')[0] == fi:
# 						# print key
# 						pi[key] = bf['parameters'][key]
# 						info_message( key + "\t" + str(pi[key]) )
# 						
# 				bf_map = setup.fg_model( pi )
# 				print setup.fg_chisq( pi ) / np.sum(setup.mask.map)
# 				map = getattr(data, fi)
# 
# 				# hp.mollview(bf_map, min=0., max=0.3, title=fi, unit='mK')
# 				# hp.graticule( dpar=20, dmer=20 )
# 
# 				ccc = (map.imap.map - bf_map)
# 				ccc[bpix] = hp.UNSEEN
# 				hp.gnomview( ccc, title='Residuals @ '+fi+' GHz', unit='mK', rot=[0,-22.5], reso=25 )
# 				hp.graticule( dpar=5, dmer=5 )
# 
# 				ccc = bf_map
# 				ccc[bpix] = hp.UNSEEN
# 				hp.gnomview( ccc, title='Solution @ '+fi+' GHz', unit='mK', rot=[0,-22.5], reso=25 )
# 				hp.graticule( dpar=5, dmer=5 )
# 		
# 			ccc = gmodel.galaxy.cmb.map
# 			ccc[bpix] = hp.UNSEEN
# 			hp.gnomview( ccc, title='CMB', unit='mK', rot=[0,-22.5], reso=25 )
# 			hp.graticule( dpar=5, dmer=5 )
# 		
# 		# pl.show()
# 
# 		# --- Analyze SEDs ---
# 		if plot_sed:
# 			f = data.cfreq
# 			cv = conversionfactor( f, option='antenna2thermo' )
# 			bf_sed = get_sed( bf['parameters'] )
# 			# print sed
# 
# 			fig = plt.figure(30, (13,8))
# 			ax = plt.subplot(111)
# 
# 			plt.errorbar( f, sed['dust'], yerr=sed_er['dust'], label='dust', color='b', fmt='s:' )
# 			plt.errorbar( f, sed['sync'], yerr=sed_er['sync'], label='sync', color='g', fmt='s:' )
# 			plt.errorbar( f, sed['haze'], yerr=sed_er['haze'], label='haze', color='m', fmt='s:' )
# 			if not fix_ff:
# 				plt.errorbar( f, sed['freefree'], yerr=sed_er['freefree'], label='freefree', color='r', fmt='s:' )
# 			if not fix_cmb:
# 				plt.errorbar( f, sed['cmb'], yerr=sed_er['cmb'], label='cmb', color='k', fmt='s:' )
# 			if inc_disc:
# 				plt.errorbar( f, sed['disc'], yerr=sed_er['disc'], label='disc', color='y', fmt='s:' )
# 
# 			plt.errorbar( f, fx_sed['cmb'], yerr=fx_sed_er['cmb'], label='cmb FS', color='k', fmt='.--' )
# 			plt.errorbar( f, fx_sed['dust'], yerr=fx_sed_er['dust'], label='dust FS', color='b', fmt='.--' )
# 			plt.errorbar( f, fx_sed['sync'], yerr=fx_sed_er['sync'], label='sync FS', color='g', fmt='.--' )
# 			plt.errorbar( f, fx_sed['haze'], yerr=fx_sed_er['haze'], label='haze FS', color='m', fmt='.--' )
# 			if not fix_ff:
# 				plt.errorbar( f, fx_sed['freefree'], yerr=fx_sed_er['freefree'], label='freefree FS', color='r', fmt='.--' )
# 			if inc_disc:
# 				plt.errorbar( f, fx_sed['disc'], yerr=fx_sed_er['disc'], label='disc FS', color='y', fmt='s--' )
# 
# 			plt.plot( f, (f/f[0])**(-3.1)*cv*fx_sed['sync'][0], ':', label='-3.1', color='g' )
# 			if not fix_ff:
# 				plt.plot( f, (f/f[0])**(-2.15)*cv*fx_sed['freefree'][0], ':', label='-2.15', color='r' )
# 			plt.plot( f, f*0.+sed['haze'][0], ':', color='m', label='same amplitude' )
# 
# 			plt.title('Dark matter model: '+model)
# 			plt.xlabel(r'$\nu$ [GHz]')
# 			plt.ylabel('Template amplitudes')
# 
# 			ax.set_yscale('log')
# 			ax.set_xscale('log')
# 			ax.set_ylim([1e-1,10])
# 			ax.set_xlim([10,100])
# 
# 
# 			plt.plot( f, bf_sed['dust'], 'd-', label='BF dust', color='b' )
# 			# plt.plot( f, sed['cmb'], label='cmb' )
# 			if not fix_ff:
# 				plt.plot( f, bf_sed['freefree'], 'd-', label='BF freefree', color='r' )
# 			plt.plot( f, bf_sed['sync'], 'd-', label='BF sync', color='g' )
# 			plt.plot( f[0], bf['parameters']['global_haze'], 'd-', label='BF haze', color='m' )
# 
# 			# plt.plot( f, (f/f[0])**(-3.1)*cv*bf_sed['sync'][0], '--', label='-3.1', color='g' )
# 			# plt.plot( f, (f/f[0])**(-2.15)*cv*bf_sed['freefree'][0], '--', label='-2.15', color='r' )
# 
# 			plt.legend( loc='upper left',prop={'size':10})
# 
# 			if make_pics:
# 				fig.savefig( "../seds_"+runtag+".png" )
	#
	# pl.show()


