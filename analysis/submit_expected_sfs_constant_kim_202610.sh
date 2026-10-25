#!/bin/bash
#$ -cwd
#$ -N constkim_sfs
#$ -l h_data=4G,h_rt=8:00:00
#$ -pe shared 8
#$ -o /u/scratch/e/eewade/logs/constkim_sfs.txt
#$ -j y
# Expected nonsyn SFS under constant size + Kim et al. 2017 DFE (Fig 2B of the 202610 manuscript figures).
. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate dadi_env
python /u/project/awbigham/eewade/daily/202610/sweeps/expected_sfs_constant_kim_202610.py 100 1000
