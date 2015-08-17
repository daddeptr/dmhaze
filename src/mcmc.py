import numpy as np
from collections import OrderedDict
import pickle
from myclass import info_message, MyClass
	
# --------------------------------------
class MCMC( MyClass ):

	# ----------------------------------
	def sampler( self, parameters, chisq_func, nsamples=10000l,
		Verbosity=False, burnin=-1, seed=0, temperature=1., output_tag='', output_dir='../chains/',
		check_point=None, accept_first=False, covmat=None, update_covmat=False ):
		# parameters should be a dictionary of dictionaries (one for each parameter)
		# containing center_value, min, max,
		# and the standard deviation (sigma) of the parameter. Samples are drawn from a Gaussian
		# distribution with the provided standard deviation
		
		# ------
		def write_chain( chain ):
			info_message( ' Saving chains...' )
			ofile = output_dir+'mcmc_chain'+output_tag+'.txt'
			#print ofile
			ch = []
			for c in chains:
				a = [ c['weight'], c['chi2'] ]
				a.extend( c['parameters'].values() )
				#print a
				ch.append( a )
			ch = np.array( ch )
			#print ch[0:2,:]
			#print ch.shape
			np.savetxt( ofile, ch, fmt='%-16.7E' )
			info_message( ' Done.\n' )				

		# ------
		def write_parnames( par_keys, latexpars=None ):
			basename = output_dir+'mcmc_chain'+output_tag+'.txt'
			basename = basename[0:-6]
			info_message( 'writing paramnames file...' )
			parnames = open( basename+'.paramnames', 'w')
			if latexpars == None:
				latexpars = {}
				for key in par_keys:
					latexpars[key] = key
			for key in par_keys:
				parnames.write( key + '\t' + latexpars[key] + '\n' )
			parnames.close()
		# ------

		code = 'sampler'

		info_message( 'Starting MCMC...' )
		info_message( '\t output_tag: '+output_tag )

		if len(output_tag) != 0 and output_tag[0] != '_':
			output_tag = '_'+output_tag
		if check_point != None:
			if check_point[0] != '_':
				check_point = '_'+check_point

		loginf = 1.e2
		zero_like = np.exp( -0.5*loginf )
		
		npars = len( parameters )
		if Verbosity:
			info_message( 'INFO - Number of parameters: '+str(npars), code=code )

		np.random.seed(seed=seed)
		step_jumps = np.random.randn( nsamples+1, npars )
		like_eval = np.random.rand( nsamples+1 )
		
		par_keys = parameters.keys()
			
		# old_parameters = {}
		old_parameters = OrderedDict()
		for ikey,key in enumerate(par_keys):
			par = parameters[key]
			if accept_first:
				old_parameters[key] = par['center_value']
			else:
				old_parameters[key] = par['center_value'] + par['sigma'] * step_jumps[0,ikey]
			if ( old_parameters[ key ] < par['min']):
				old_parameters[ key ] = par['center_value']
				if ( old_parameters[ key ] < par['min']):
					old_parameters[ key ] = par['min']
			if ( old_parameters[ key ] > par['max']):
				old_parameters[ key ] = par['center_value']
				if ( old_parameters[ key ] > par['max']):
					old_parameters[ key ] = par['max']
					
		# --- Reading covmat if provided
		use_cov = False # Is set to True later on in the code if update_covmat is True
		if isinstance(covmat,str):
			try:
				# print covmat
				f = open( covmat )
				covmat = pickle.load( f )
				f.close()
			except:
				info_message( 'ERROR - Error loading the covariance matrix file:'+covmat+'. Abort.' )
				sys.exit()

		if covmat != None:
			info_message( 'INFO - Covariance matrix provided as array.' )
			if (covmat.ndim==2) and (covmat.shape[0]==covmat.shape[1]) and (covmat.shape[0]==npars):
				use_cov = True
				for ikey,key in enumerate(par_keys):
					# old_pars_matrix[ikey,ikey] = old_parameters[ key ]
					use_cov = ( use_cov and parameters[key]['sigma'] != 0. )
				if use_cov:
					w, v = self.get_eig( covmat )
					# vm1 = np.linalg.inv( v )
					wsig = np.sqrt( w )
					S = np.zeros( (npars,npars),dtype=np.float64)
					for ip in range(npars):
						S[ip,ip] = np.sqrt( w[ip] )
					Sm1 = np.linalg.inv( S )
				else:
					info_message( 'WARNING - Something went wrong with thte covariance matrix provided.\nMaybe some parameters were fixed. Ignored.' )
					
		counter = 0
		accepted = 0
		sample_check = 0
		skipped = 0
		nskp = 0
		nout_bounds = 0

		chains = []		
		old_chi2 = chisq_func(old_parameters)
		old_like = np.exp( -0.5 * old_chi2 / temperature )
		#print old_parameters
		best_fit = {'chisq': old_chi2, 'parameters': old_parameters}
		if Verbosity:
			print ' --> ', best_fit

		if check_point != None:
			info_message( 'INFO - Starting from last point of the chain...' )
			try:
				filename = '../mcmc_chain'+check_point+'.pysav'
				info_message( filename )
				f = open( filename, 'r' )
				best_fit, ch = pickle.load( f )
				f.close()
				old_chi2 = best_fit['chisq']
				old_parameters = best_fit['parameters']
				burnin = 0
			except:
 				info_message( 'ERROR - Error loading check_point: '+filename+'. Skip.' )
				
		if burnin >= 0:
			info_message( 'INFO - burnin set to '+str(burnin) )
			
		# --- Loop over samples
		for isample in range(1,nsamples+1):
			# --- Print some statistics
			if ( (isample % (nsamples/20)) == 0 ):
				info_message( ' steps completed ----- '+str( int(float(isample)/nsamples*100))+'%' )
				info_message( ' sample_check = ' + str(sample_check) )
				info_message( ' nsamples     = ' + str(nsamples) )
				info_message( ' skipped      = ' + str(skipped) )
				info_message( ' nout_bounds  = ' + str(nout_bounds) )
				info_message( ' chain length = ' + str(len(chains)) )
				info_message( ' acceptance_rate = ' + str(float(len(chains))/(isample)*100) )
				info_message( ' burnin       = ' + str(burnin) )
				print best_fit
# 				info_message( ' Saving chains...' )
# 				ch = (best_fit, chains)
# 				f = open( output_dir+'/mcmc_chain'+output_tag+'.pysav', 'w' )
# 				pickle.dump( ch, f )
# 				f.close()
# 				self.py2getdist( output_dir+'/mcmc_chain'+output_tag+'.pysav' )
# 				info_message( ' Done.\n' )				
				## :DP May 2015 - OoM fix
				write_chain( chains )
				write_parnames( par_keys )
							
			counter += 1
			
			# --- draw new parameters
			new_parameters = OrderedDict()
			out_of_bounds = False

			if not use_cov:
				for ikey,key in enumerate(par_keys):
					if parameters[key]['sigma'] != 0.:
						new_parameters[key] = old_parameters[key] + parameters[key]['sigma'] * step_jumps[isample,ikey]
			else:
				# info_message( 'Here' )
				old_pars = np.zeros( npars, dtype=np.float64 )

				# step_diag = wsig * step_jumps[isample,:]
				step_diag = np.dot( np.dot( np.dot( step_jumps[isample,:], v ), S), v.T )

				for ikey,key in enumerate(par_keys):
					old_pars[ikey] = old_parameters[ key ]

				# old_pars = np.dot( np.dot( np.dot( old_pars, v ), S), v.T )
				new_pars = old_pars + step_diag
				# new_pars = np.dot( np.dot( np.dot( new_pars_diag, v.T ), Sm1), v ) 

				for ikey,key in enumerate(par_keys):
					new_parameters[key] = new_pars[ikey]
					
				# print new_pars
				if False:
					print ' old <-- ', old_parameters
					print ' -   ', old_chi2
					print ' new --> ', new_parameters
					print  ' -   ', chisq_func( new_parameters )
				# sys.exit()

			for ikey,key in enumerate(par_keys):
				if new_parameters[key] < parameters[key]['min']:
					# new_parameters[key] = parameters[key]['min']
					# out_of_bounds = True
					nout_bounds += 1
					# print key, new_parameters[key], parameters[key]['min']
					new_parameters[key] = old_parameters[key]
					# break
				if new_parameters[key] > parameters[key]['max']:
					# new_parameters[key] = parameters[key]['max']
					# out_of_bounds = True
					nout_bounds += 1
					# print key, new_parameters[key], parameters[key]['max']
					new_parameters[key] = old_parameters[key]
					# break
			
			if not out_of_bounds:
				# --- compute new likelihood
				new_chi2 = chisq_func( new_parameters )
				new_like = np.exp( -0.5 * new_chi2 / temperature )
				#print new_like
			
				# --- compare likelihoods
	# 			like_ratio = new_like/old_like
				delta_chi2 = (new_chi2-old_chi2)
				like_ratio = np.exp( -0.5 * delta_chi2 / temperature )
				if delta_chi2 > 0 and False:
					info_message( 'like_ratio = '+str(delta_chi2)+', '+str(like_ratio)+', '+str(like_eval[isample]) )
					info_message( str(delta_chi2 / temperature ) )
	# 			print old_chi2, new_chi2, old_chi2-new_chi2
	# 			print old_like, new_like, np.exp( -0.5*(new_chi2-old_chi2) )

				if (like_ratio >= like_eval[isample]):
	# 			if (like_ratio >= like_eval[isample]) or (new_chi2 <= old_chi2):
					accepted += 1
					#print old_parameters
					#print new_parameters
					#print ' --> ', like_ratio, like_eval[isample]
					old_parameters = new_parameters
					old_like = new_like
					old_chi2 = new_chi2
					sample_check += counter
					if accepted > burnin:
						chains.append( {'weight':counter, 'chi2':new_chi2, 'parameters':new_parameters} )
						#print chains[-1]
						#print ''
						if new_chi2 < best_fit['chisq']:
							best_fit = {'chisq': new_chi2, 'parameters': old_parameters}
							# print '\n --> ', isample, best_fit['chisq']
							# info_message( ' steps completed ----- '+str( int(float(isample)/nsamples*100))+'%' )
							# info_message( ' sample_check = ' + str(sample_check) )
							# info_message( ' nsamples     = ' + str(nsamples) )
							# info_message( ' skipped      = ' + str(skipped) )
							# info_message( ' nout_bounds  = ' + str(nout_bounds) )
							# info_message( ' chain length = ' + str(len(chains)) )
							# info_message( ' acceptance_rate = ' + str(float(len(chains))/(isample)*100) )
							# info_message( ' burnin       = ' + str(burnin) )
							
					# --- Update covariance matrix
					if use_cov and update_covmat and ( (accepted % 1000 ) == 0 ):
						info_message( 'INFO - Accepted '+str(accepted)+': updating covariance matrix...' )
						covmat = self.get_covmat( chains, skiplines=0 )
						w, v = self.get_eig( covmat )
						# vm1 = np.linalg.inv( v )
						wsig = np.sqrt( w )
						S = np.zeros( (npars,npars),dtype=np.float64)
						for ip in range(npars):
							S[ip,ip] = np.sqrt( w[ip] )
						Sm1 = np.linalg.inv( S )	
						info_message( 'Done.' )
							
					# --- Use covariance matrix if not provided
					if (covmat == None) and update_covmat and (accepted == 100 ):
						use_cov = True
						info_message( 'INFO - Accepted '+str(accepted)+': updating covariance matrix...' )
						covmat = self.get_covmat( chains, skiplines=0 )
						w, v = self.get_eig( covmat )
						# vm1 = np.linalg.inv( v )
						wsig = np.sqrt( w )
						S = np.zeros( (npars,npars),dtype=np.float64)
						for ip in range(npars):
							S[ip,ip] = np.sqrt( w[ip] )
						Sm1 = np.linalg.inv( S )	
						info_message( 'Done.' )
							
					counter = 0
					nskp = 0
				else:
					nskp += 1
					skipped += 1
			else:
				nskp += 1
				skipped += 1
				
		sample_check += counter
		info_message( ' - sample_check = ' + str(sample_check) )
		info_message( ' - nsamples     = ' + str(nsamples) )
		info_message( ' - skipped      = ' + str(skipped) )
		info_message( ' - nout_bounds  = ' + str(nout_bounds) )
		info_message( ' - chain length = ' + str(len(chains)) )
		info_message( ' - acceptance_rate = ' + str(float(len(chains))/nsamples*100) )
		info_message( ' - burnin       = ' + str(burnin) )

# 		ch = (best_fit, chains)
# 		f = open( output_dir+'/mcmc_chain'+output_tag+'.pysav', 'w' )
# 		pickle.dump( ch, f )
# 		f.close()
# 		self.py2getdist( output_dir+'/mcmc_chain'+output_tag+'.pysav' )
		## :DP - May 2015 OoM fix
		write_chain( chains )

		print best_fit
		#print chains[-2:]
		f = open( output_dir+'/mcmc_chain_bestfit'+output_tag+'.txt', 'w' )
		f.write( 'chisq = '+str(best_fit['chisq'])+'\n' )
		for key in np.sort(best_fit['parameters'].keys()):
			f.write( key + ' = '+str(best_fit['parameters'][key])+'\n' )
		f.close()

		return best_fit, chains

	# ----------------------------------
	def marginalize( self, chains, nbins=21, skiplines=0.3, thin=1):
		# ------------------------------
		def extract_parameter( chains, par_key):
			w = []
			v = []
			clen = len( chains )
			for i in range(clen):
				w.append( chains[i]['weight'])
				v.append( chains[i]['parameters'][par_key])
			v = np.array(v)
			w = np.array(w)
			return v,w
			
		# ------------------------------
		def compute_histogram( value, weight, nbins=21):
			minv = np.min( value )
			maxv = np.max( value )
			deltav = (maxv-minv)/(nbins)
			intervals = np.array(range(nbins+1)) * deltav + minv
			#print intervals
			v = []
			m = []
			for ibin in range(nbins):
				v.append( (intervals[ibin]+intervals[ibin+1])/2. )
				select = np.array( (value >= intervals[ibin]) ) * np.array( (value < intervals[ibin+1]) )
				m.append( np.sum(weight[select]) )
			v = np.array(v)
			m = np.array(m, dtype=np.float) / np.sum(m)
			return v, m

		chlen = len( chains )
		if isinstance(skiplines,float):
			chlen = int( chlen * skiplines )		
		elif isinstance(skiplines, int):
			chlen = skiplines
		chains = chains[chlen::thin]
		info_message( 'marginalize: Processing '+str(len(chains))+' samples' )
		npars = len( chains[0]['parameters'] )
		par_keys = chains[0]['parameters'].keys()
		
		marginals = {}
		for key in par_keys:
			v,w = extract_parameter( chains, key )
			v,m = compute_histogram( v, w, nbins=nbins )
			marginals[key] = (v,m)
		return marginals

	# ----------------------------------
	def get_covmat( self, chains, skiplines=0.3, thin=1, print_key=False ):
		chlen = len( chains )
		if isinstance(skiplines,float):
			chlen = int( chlen * skiplines )		
		elif isinstance(skiplines, int):
			chlen = skiplines
		chains = chains[chlen::thin]
		chlen = len( chains )
		info_message( 'get_covmat: Processing '+str(chlen)+' samples' )
		npars = len( chains[0]['parameters'] )
		par_keys = chains[0]['parameters'].keys()

		covmat = np.zeros( (npars,npars), dtype=np.float64 )
		meanpars = np.zeros( npars, dtype=np.float64 )
		
		# Compute mean parameter values and reshaping parameter chains
		pchains = []
		# print_key = False
		for ch in chains:
			prow = []
			for ikey,key in enumerate(par_keys):
				if print_key:
					print key
				prow.append( ch['parameters'][key] )
				meanpars[ikey] += ch['parameters'][key]
			pchains.append( prow )
			if print_key:
				print_key = False

		meanpars /= chlen
		pchains = np.array( pchains )
		pchains.shape
		for ich in range(chlen):
			for ipar in range(npars):
				covmat[ipar,:] += (pchains[ich,ipar]-meanpars[ipar]) * (pchains[ich,:]-meanpars[:])		

		covmat /= chlen
		
		return covmat
	
	# ----------------------------------
	def get_eig( self, covmat ):
		import numpy as np
		w, vr = np.linalg.eigh( covmat )

		return w, vr
		
	# ----------------------------------
# 	def py2getdist( self, chains_file, latexpars=None ):
# 
# 		info_message( ' >>> py2getdist: processing "'+chains_file+'"...' )
# 
# 		info_message( 'importing chains...' )
# 		chains = open( chains_file, 'r' )
# 		bf, ch = pickle.load( chains )
# 		chains.close
# 
# 		info_message( 'chain length = '+str(len(ch)) )
# 
# 		if latexpars == None:
# 			latexpars = {}
# 			for key in ch[0]['parameters'].keys():
# 				latexpars[key] = key
# 
# 		info_message( 'saving chains in getdist format...' )
# 		basename = chains_file[0:-6]
# 		gd_file = open( basename+'.txt', 'w')
# 		pars = bf['parameters'].keys()
# 
# 		gd_chains = []
# 		first_sample = True
# 		for sample in ch:
# 			smp = []
# 			smp.append( str( '%13.7e' %sample['weight'] ) )
# 			smp.append( str( '%13.7e' %sample['chi2']/2. ) )
# 			for key in pars:
# 				smp.append( str( '%13.7e' %sample['parameters'][key] ) )
# 
# 			gd_chains.append(smp)
# 			gd_file.write( '\t'.join( smp )+'\n' )
# 
# 			if first_sample:
# 				info_message( 'writing paramnames file...' )
# 				parnames = open( basename[0:-2]+'.paramnames', 'w')
# 				for key in pars:
# 					parnames.write( key + '\t' + latexpars[key] + '\n' )
# 				parnames.close()
# 				first_sample = False
# 
# 		gd_file.close()

	
