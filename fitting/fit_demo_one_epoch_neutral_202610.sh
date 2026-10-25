#!/bin/bash
#$ -N dadi202610_1ep.$TASK_ID
#$ -o /u/scratch/e/eewade/logs/dadi202610_1ep.$TASK_ID.txt
#$ -j y
#$ -cwd
#$ -l h_rt=4:00:00,h_vmem=8G
# 240 tasks = 40 neutral human-DFE run_names (2 recombination rates x 20 replicates)
#             x 2 sample sizes x 3 distances (1, 3, 5 Mb)
#$ -t 1-240

# Constant-size (one_epoch, dadi snm) demography fit on the syn SFS of the
# single-sweep neutral controls.
# Same SFS conversion, L_syn and dadi settings as fit_demo_dfe_human_202610.sh.
# Folded runs: qsub -v FOLD_MODE=fold -N dadi202610_1ep_fold \
#   -o '/u/scratch/e/eewade/logs/dadi202610_1ep_fold.$TASK_ID.txt' fit_demo_one_epoch_neutral_202610.sh
#
# Output: data/working-data/202610/dadi_analysis[_folded]/<run_name>/n<n>_d<d>mb/one_epoch/demography/

. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate dadi_env

FIT_DEMOG_SCRIPT=/u/project/awbigham/eewade/daily/202610/sweeps/fit_demo_fixedgrid_test.py
ARGS_FILE=/u/project/awbigham/eewade/daily/202610/sweeps/single_task_args.txt
SLIM_OUTDIR=/u/home/e/eewade/project-awbigham/data/working-data/202610/tracksweep_output
ANALYSIS_OUTDIR=/u/home/e/eewade/project-awbigham/data/working-data/202610/dadi_analysis

FOLD_MODE=${FOLD_MODE:-nofold}
case "$FOLD_MODE" in
  nofold) FOLD_FLAG=--no-fold ;;
  fold)   FOLD_FLAG=--fold; ANALYSIS_OUTDIR=${ANALYSIS_OUTDIR}_folded ;;
  *) echo "ERROR: FOLD_MODE must be fold or nofold, got '$FOLD_MODE'"; exit 1 ;;
esac

EXON_SYN_FRACTION=0.302
PTS="1200 1600 2000"
MODEL=one_epoch

SAMPLE_SIZES=(100 1000)
DISTANCES_MB=(1 3 5)

mapfile -t RUN_NAMES < <(awk '!seen[$4]++ {print $4}' "$ARGS_FILE" | grep '_neutral_dfehuman_' | sort -V)

N_RUNS=${#RUN_NAMES[@]}
N_SIZES=${#SAMPLE_SIZES[@]}
N_DISTS=${#DISTANCES_MB[@]}
N_PER_RUN=$((N_SIZES * N_DISTS))

IDX=$((SGE_TASK_ID - 1))
RUN_IDX=$((IDX / N_PER_RUN))
REM=$((IDX % N_PER_RUN))
SIZE_IDX=$((REM / N_DISTS))
DIST_IDX=$((REM % N_DISTS))

if [ "$RUN_IDX" -ge "$N_RUNS" ]; then
  echo "ERROR: task $SGE_TASK_ID out of range ($N_RUNS run_names x $N_PER_RUN combos = $((N_RUNS * N_PER_RUN)))"
  exit 1
fi

RUN_NAME="${RUN_NAMES[$RUN_IDX]}"
SAMPLE_SIZE="${SAMPLE_SIZES[$SIZE_IDX]}"
DIST_MB="${DISTANCES_MB[$DIST_IDX]}"

OUTDIR="$ANALYSIS_OUTDIR/$RUN_NAME/n${SAMPLE_SIZE}_d${DIST_MB}mb/${MODEL}"
DEMOG_OUTDIR="$OUTDIR/demography"
DEMOG_FILE="$DEMOG_OUTDIR/${MODEL}_demography.txt"

echo "Task $SGE_TASK_ID -> $RUN_NAME n=$SAMPLE_SIZE d=${DIST_MB}mb model=$MODEL fold_mode=$FOLD_MODE"

if grep -q "Best fit parameters:" "$DEMOG_FILE" 2>/dev/null; then
  echo "Already complete ($DEMOG_FILE) -- skipping."
  exit 0
fi

syn_file="$SLIM_OUTDIR/syn_sfs_${RUN_NAME}_n${SAMPLE_SIZE}_d${DIST_MB}mb.tsv"
bp_file="$SLIM_OUTDIR/genic_bp_${RUN_NAME}.csv"
for f in "$syn_file" "$bp_file"; do
  if [ ! -s "$f" ]; then echo "ERROR: missing $f"; exit 1; fi
done

# Sum segments and convert to dadi's unfolded Spectrum format.
mkdir -p "$DEMOG_OUTDIR"
SYN_DADI="$OUTDIR/${RUN_NAME}_syn_n${SAMPLE_SIZE}_d${DIST_MB}mb.sfs"
awk -v n="$SAMPLE_SIZE" '
  NF == 0 { next }
  { for (i = 1; i <= NF; i++) sum[i] += $i; ncols = NF }
  END {
    if (ncols != n) { print "ERROR: expected " n " columns, got " ncols > "/dev/stderr"; exit 1 }
    printf "%d unfolded\n0", n + 1
    for (i = 1; i <= ncols; i++) printf " %s", sum[i]
    printf "\n"
  }' "$syn_file" > "$SYN_DADI" || exit 1

# L_syn = exon_bp * 0.302 + intron_bp, summed over the 22 segments.
N_SEG_BP=$(awk -F',' -v d="$DIST_MB" '$1 != "segment" && $2 == d' "$bp_file" | wc -l)
if [ "$N_SEG_BP" -ne 22 ]; then
  echo "ERROR: expected 22 segments at d=${DIST_MB}mb in $bp_file, found $N_SEG_BP"
  exit 1
fi
L_SYN=$(awk -F',' -v frac="$EXON_SYN_FRACTION" -v d="$DIST_MB" '
  $1 != "segment" && $2 == d { ex += $3; in_ += $4 }
  END { printf "%.0f\n", ex * frac + in_ }' "$bp_file")
echo "L_syn=$L_SYN"

python "$FIT_DEMOG_SCRIPT" "$SYN_DADI" --model_type "$MODEL" --L_syn "$L_SYN" \
  "$FOLD_FLAG" --pts $PTS "$DEMOG_OUTDIR/"
if ! grep -q "Best fit parameters:" "$DEMOG_FILE" 2>/dev/null; then
  echo "ERROR: demography fit failed -- see $DEMOG_OUTDIR/log.log"
  exit 1
fi

echo "=== Finished $RUN_NAME n=$SAMPLE_SIZE d=${DIST_MB}mb model=$MODEL (L_syn=$L_SYN) ==="
