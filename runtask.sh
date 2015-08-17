#!/bin/bash
cd /global/scratch2/sd/dpietrob/DM-haze/src/
pwd

source /project/projectdirs/cmb/modules/hpcports_NERSC.sh
module unload openmpi/1.4.5 pgi/12.9
hpcports gnu
module load astromatic-hpcp
## module load matplotlib-hpcp
module load scipy-hpcp
module load cmb
echo "carver environment loaded"
module list

nch=$(printf "%1.1d" $(($PBS_VNODENUM+1 )) )     # 1, 2, 3,...
## echo "My chain number: "$nch
model=1115
## tag="bubbles_FINAL_long"
## :DP May 2015 to get Col.4 of Tab.3
tag="DM_FINAL_long"

tmpfile="/global/scratch2/sd/dpietrob/DM-haze/distparams_${tag}_${model}.ini"

cp /global/scratch2/sd/dpietrob/DM-haze/distparams.ini ${tmpfile}

sed "6 c\ file_root = chains/mcmc_chain_${tag}_${model}_2stp_08x" -i ${tmpfile}
sed "7 c\ out_root = mcmc_chain_${tag}_${model}_2stp_08x" -i ${tmpfile}
sed "13 c\ chain_num = ${PBS_NUM_PPN}" -i ${tmpfile}

echo "python run_final_mcmc.py -m ${model} -s ${nch} -o ${tag}"
## which python
python run_final_mcmc.py -m $model -s $nch -o $tag
