#!/bin/bash
#$ -N recurrent202610.$TASK_ID
#$ -o /u/scratch/e/eewade/logs/recurrent202610.$TASK_ID.txt
#$ -j y
#$ -cwd
#$ -l h_rt=23:00:00,h_vmem=5G
#$ -pe share 4
#$ -M $USER@mail
#$ -m a
# 3520 tasks = 4 (s,p) pairs x 2 r x 1 dfe x 20 replicates x 22 segments (generate_recurrent_task_args_202610.sh)
#$ -t 1-3520
. /u/local/Modules/default/init/modules.sh
module load slim/5.1

# Fresh output dir: SFS files are appended per segment, never reuse an old dir
PATH_OUT=/u/scratch/e/eewade/recurrentsweep_output_202610
mkdir -p "$PATH_OUT"

SLIM_SCRIPT=/u/project/awbigham/eewade/daily/202610/sweeps/RecurrentSweep_v2_genic.slim
ARGS_FILE=/u/project/awbigham/eewade/daily/202610/sweeps/recurrent_task_args.txt
LINE=$(awk -v id="$SGE_TASK_ID" '$1==id {print}' "$ARGS_FILE")
echo "$LINE"

REPLICATE=$(echo $LINE | awk '{print $2}')
SEGMENT=$(echo $LINE | awk '{print $3}')
RUN_NAME_OUTPUT=$(echo $LINE | awk '{print $4}')
MU=$(echo $LINE | awk '{print $5}')
BEN_S=$(echo $LINE | awk '{print $6}')
BEN_FRAC=$(echo $LINE | awk '{print $7}')
RECOM_RATE=$(echo $LINE | awk '{print $8}')
SEQLENGTH=$(echo $LINE | awk '{print $9}')
FIXED_GENS=$(echo $LINE | awk '{print $10}')
DFE_MEAN=$(echo $LINE | awk '{print $11}')
DFE_SHAPE=$(echo $LINE | awk '{print $12}')

slim -d "run_name_output='${RUN_NAME_OUTPUT}'" -d "segment=${SEGMENT}" -d MU=${MU} -d BEN_S=${BEN_S} -d BEN_FRAC=${BEN_FRAC} -d RECOM_RATE=${RECOM_RATE} -d SEQLENGTH=${SEQLENGTH} -d FIXED_GENS=${FIXED_GENS} -d DFE_MEAN=${DFE_MEAN} -d DFE_SHAPE=${DFE_SHAPE} -d "PATH='${PATH_OUT}'" "$SLIM_SCRIPT"

# gzip this segment's n=100 VCF
v="${PATH_OUT}/vcf_${RUN_NAME_OUTPUT}_seg${SEGMENT}.vcf"
if [ -e "$v" ]; then gzip -f "$v"; fi
