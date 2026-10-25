#!/bin/bash
#$ -N tracksweep202610.$TASK_ID
#$ -o /u/scratch/e/eewade/logs/tracksweep202610.$TASK_ID.txt
#$ -j y
#$ -cwd
# #$ -l highp
#$ -l h_rt=23:00:00,h_vmem=5G
#$ -pe share 4
#$ -M $USER@mail
# Notify when
#$ -m bea
# 5280 tasks = 12 param sets (8 sweep + 4 neutral) x 20 replicates x 22 segments
#$ -t 1-5280
. /u/local/Modules/default/init/modules.sh
module load slim/5.1

# Output and temp-state directories (same for every task -- not varied per replicate/segment)
PATH_OUT=/u/scratch/e/eewade/tracksweep_output_202610
TEMPDIR_OUT=/u/scratch/e/eewade/tracksweep_temp_202610
mkdir -p "$PATH_OUT" "$TEMPDIR_OUT"

SLIM_SCRIPT=/u/project/awbigham/eewade/daily/202610/sweeps/ConstantPopSizeOneSweep_v4_genic.slim
ARGS_FILE=/u/project/awbigham/eewade/daily/202610/sweeps/single_task_args.txt
LINE=$(awk -v id="$SGE_TASK_ID" '$1==id {print}' $ARGS_FILE)
echo $LINE

REPLICATE=$(echo $LINE | awk '{print $2}')
SEGMENT=$(echo $LINE | awk '{print $3}')
RUN_NAME_OUTPUT=$(echo $LINE | awk '{print $4}')
MU=$(echo $LINE | awk '{print $5}')
BEN_S=$(echo $LINE | awk '{print $6}')
RECOM_RATE=$(echo $LINE | awk '{print $7}')
SEQLENGTH=$(echo $LINE | awk '{print $8}')
FIXED_GENS=$(echo $LINE | awk '{print $9}')
DFE_MEAN=$(echo $LINE | awk '{print $10}')
DFE_SHAPE=$(echo $LINE | awk '{print $11}')

echo $REPLICATE
echo $SEGMENT
echo $RUN_NAME_OUTPUT
echo $MU
echo $BEN_S
echo $RECOM_RATE
echo $SEQLENGTH
echo $FIXED_GENS
echo $DFE_MEAN
echo $DFE_SHAPE

slim -d "run_name_output='${RUN_NAME_OUTPUT}'" -d "segment=${SEGMENT}" -d MU=${MU} -d BEN_S=${BEN_S} -d RECOM_RATE=${RECOM_RATE} -d SEQLENGTH=${SEQLENGTH} -d FIXED_GENS=${FIXED_GENS} -d DFE_MEAN=${DFE_MEAN} -d DFE_SHAPE=${DFE_SHAPE} -d "PATH='${PATH_OUT}'" -d "TEMPDIR='${TEMPDIR_OUT}'" ${SLIM_SCRIPT}

# gzip this segment's pre-/post-sweep n=100 VCFs (vcf_<run>_seg<seg>_{pre,post}.vcf)
for v in "${PATH_OUT}/vcf_${RUN_NAME_OUTPUT}_seg${SEGMENT}_pre.vcf" "${PATH_OUT}/vcf_${RUN_NAME_OUTPUT}_seg${SEGMENT}_post.vcf"; do
  if [ -e "$v" ]; then gzip -f "$v"; fi
done
