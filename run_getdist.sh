#!/bin/bash

#module swap pgi intel ; module swap openmpi openmpi-intel ; module load mkl ; module load python ; module load cython ; module load idl ; export CFITSIO_LIB=/global/scratch2/sd/dpietrob/Software/cfitsio_ifort; echo Intel loaded; export HEALPIX=$GSCRATCH2"/Tools/Healpix_3.00"; source /global/homes/d/dpietrob/.healpix/3_00_Linux/idl.sh; echo myHealpix loaded
export getdist='/global/scratch2/sd/dpietrob/Software/xfpipe/cosmomc_cosmoslik/getdist'
which $getdist

#tag="bubbles_FINAL_long"
tag="DM_FINAL_long"
echo $tag

#for mass in {1..5}
#do
#	for channel in 1 2 #2
#	do
#		for slope in {1..3}
#		do
#			for mag in 13 21 32 4 5 #13 21 32 4 5
#			do
#				for pw in "" "h" #13 21 32 4 5
#				do
#					model=${mass}${channel}${slope}${mag}${pw}
				for model in 12132h #11113 
				do
					echo $model
					# ls "data/maps/ns0128/haze_model_54*${model}*"
					tmpfile="/global/scratch2/sd/dpietrob/DM-haze/distparams_${tag}_${model}.ini"
					echo $tmpfile
					## To fix -loglike error "fixed_" added
#					cp default.paramnames chains/fixed_mcmc_chain_${tag}_${model}_2stp_08x.paramnames
#					sed "6 c\file_root = chains/fixed_mcmc_chain_${tag}_${model}_2stp_08x" -i ${tmpfile}
#					sed "7 c\out_root = fixed_mcmc_chain_${tag}_${model}_2stp_08x" -i ${tmpfile}
					cp default.paramnames chains/fixed_mcmc_chain_${tag}_${model}_2stp_08x.paramnames
					sed "/bubbles/d" -i chains/fixed_mcmc_chain_${tag}_${model}_2stp_08x.paramnames
					## head chains/fixed_mcmc_chain_${tag}_${model}_2stp_08x.paramnames
#					sed -i "5 s/^/#/" chains/fixed_mcmc_chain_${tag}_${model}_2stp_08x.paramnames
					sed "6 c\file_root = chains/fixed_mcmc_chain_${tag}_${model}_2stp_08x" -i ${tmpfile}
					sed "7 c\out_root = fixed_mcmc_chain_${tag}_${model}_2stp_08x" -i ${tmpfile}
					$getdist $tmpfile
				done
#				done
#			done
#		done
#	done
#done