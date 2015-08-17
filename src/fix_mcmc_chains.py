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
def fix_mcmc_chains( tag = '' ):

	ldir = "../chains"
	files = [ f for f in os.listdir( ldir ) if f.endswith(".txt") and f.startswith( tag ) ]
	print " - fixing %i files..." %len( files )
#	print files[0:10]
	for f in files:
		#if not (iif % 50):
## problem: 341, mcmc_chain_bubbles_FINAL_long_22113_2stp_08x_1.txt
		print f
		a = np.loadtxt( ldir+"/"+f )
		print a.shape
		print a[0,:]
		## Commented for DM
		#a[:,1] /= 2.
# 		np.savetxt( ldir+"/fixed_"+f, a, fmt='%-16.7E', delimiter='' )
		np.savetxt( ldir+"/fixed_"+f, a, fmt='%-16.7E' )

	return True
	
# --------------------------------------------------------------------
def main(argv):
	code = "fix_mcmc_chains"

	paramfile = {}
	try:
		opts, args = getopt.getopt(argv,"o:hv",["dm_model=","seed=","outtag=","help","verbose"])
	except getopt.GetoptError:
		info_message( "unknown argument.\nCalling sequence\nfix_mcmc_chains.py -o <outtag> [-vh]", code=code )
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

	if ('o' not in paramfile.keys()):
		info_message( "output tag missing", code=code)
		paramfile['o'] = 'mcmc_chain'

	print paramfile

	fix_mcmc_chains( tag=paramfile['o'] )

# --------------------------------------------------------------------
if __name__ == "__main__":

	## python run_final_mcmc.py -o $tag

	main(sys.argv[1:])
	
'''
model=1115
tag="bubbles_FINAL_long"

tmpfile="/global/scratch2/sd/dpietrob/DM-haze/distparams_${tag}_${model}.ini"

cp /global/scratch2/sd/dpietrob/DM-haze/distparams.ini ${tmpfile}

sed "6 c\ file_root = chains/mcmc_chain_${tag}_${model}_2stp_08x" -i ${tmpfile}
sed "7 c\ out_root = mcmc_chain_${tag}_${model}_2stp_08x" -i ${tmpfile}

#!/bin/bash

mkdir -v scripts

for mass in {1..5}
do
	for channel in 1 2 #2 #1
	do
		for slope in {1..3}
		do
			for mag in 13 21 32 4 5 #13 21 32 4 5
			do
				model=${mass}${channel}${slope}${mag}"h"
				echo $model
				# ls "data/maps/ns0128/haze_model_54*${model}*"
				cp test_run.pbs scripts/run_$model.pbs
				tmpfile="scripts/run_${model}.pbs"
				sed "9 c\#PBS -N mcmc_bb_${model}" -i ${tmpfile}
				sed "10 c\#PBS -e mcmc_bb_${model}.err" -i ${tmpfile}
				sed "11 c\#PBS -o mcmc_bb_${model}.out" -i ${tmpfile}
				sed "19 c\pbsdsh /global/scratch2/sd/dpietrob/DM-haze/scripts/runtask_${model}.sh" -i ${tmpfile}

				cp runtask.sh scripts/runtask_$model.sh
				tmpfile="scripts/runtask_$model.sh"
				sed "17 c\model=${model}" -i ${tmpfile}
				
				qsub scripts/run_$model.pbs
			done
		done
	done
done
'''