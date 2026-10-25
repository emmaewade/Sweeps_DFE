#!/bin/bash
#$ -N dadi202610.$TASK_ID
#$ -o /u/scratch/e/eewade/logs/dadi202610.$TASK_ID.txt
#$ -j y
#$ -cwd
#$ -l highp,h_rt=23:59:00,h_vmem=16G
#$ -pe share 4
# 1440 tasks = 120 human-DFE run_names (6 conditions x 20 replicates)
#              x 2 sample sizes x 3 distances (1, 3, 5 Mb) x 2 demographic models
#$ -t 1-1440


. /u/local/Modules/default/init/modules.sh
module load anaconda3
conda activate dadi_env

FIT_DEMOG_SCRIPT=/u/project/awbigham/eewade/daily/202610/sweeps/fit_demo_fixedgrid_test.py
FIT_DFE_SCRIPT=/u/project/awbigham/eewade/daily/202610/sweeps/fit_dfe_poisson_test.py
ARGS_FILE=/u/project/awbigham/eewade/daily/202610/sweeps/single_task_args.txt
SLIM_OUTDIR=/u/home/e/eewade/project-awbigham/data/working-data/202610/tracksweep_output
ANALYSIS_OUTDIR=/u/home/e/eewade/project-awbigham/data/working-data/202610/dadi_analysis

# Folded vs unfolded. Set at submission: qsub -v FOLD_MODE=fold ...
FOLD_MODE=${FOLD_MODE:-nofold}
case "$FOLD_MODE" in
  nofold) FOLD_FLAG=--no-fold ;;
  fold)   FOLD_FLAG=--fold; ANALYSIS_OUTDIR=${ANALYSIS_OUTDIR}_folded ;;
  *) echo "ERROR: FOLD_MODE must be fold or nofold, got '$FOLD_MODE'"; exit 1 ;;
esac

EXON_SYN_FRACTION=0.302
PTS="1200 1600 2000"

SAMPLE_SIZES=(100 1000)
DISTANCES_MB=(1 3 5)
MODEL_TYPES=(two_epoch three_epoch)

# Human-DFE run_names (sweep and neutral, all replicates), in a fixed order.
# If RUN_NAME is passed (qsub -v RUN_NAME=...,FOLD_MODE=... -t 1-12), fit only that
# replicate: tasks 1-12 = its 2 sample sizes x 3 distances x 2 models. This is how
# dispatch_fits_202610.sh queues each replicate as soon as its 22 segments finish.
if [ -n "$RUN_NAME" ]; then
  RUN_NAMES=("$RUN_NAME")
else
  mapfile -t RUN_NAMES < <(awk '!seen[$4]++ {print $4}' "$ARGS_FILE" | grep '_dfehuman_' | sort -V)
fi

N_RUNS=${#RUN_NAMES[@]}
N_SIZES=${#SAMPLE_SIZES[@]}
N_DISTS=${#DISTANCES_MB[@]}
N_MODELS=${#MODEL_TYPES[@]}
N_PER_RUN=$((N_SIZES * N_DISTS * N_MODELS))

IDX=$((SGE_TASK_ID - 1))
RUN_IDX=$((IDX / N_PER_RUN))
REM=$((IDX % N_PER_RUN))
SIZE_IDX=$((REM / (N_DISTS * N_MODELS)))
REM2=$((REM % (N_DISTS * N_MODELS)))
DIST_IDX=$((REM2 / N_MODELS))
MODEL_IDX=$((REM2 % N_MODELS))

if [ "$RUN_IDX" -ge "$N_RUNS" ]; then
  echo "ERROR: task $SGE_TASK_ID out of range ($N_RUNS run_names x $N_PER_RUN combos = $((N_RUNS * N_PER_RUN)))"
  exit 1
fi

RUN_NAME="${RUN_NAMES[$RUN_IDX]}"
SAMPLE_SIZE="${SAMPLE_SIZES[$SIZE_IDX]}"
DIST_MB="${DISTANCES_MB[$DIST_IDX]}"
MODEL="${MODEL_TYPES[$MODEL_IDX]}"

OUTDIR="$ANALYSIS_OUTDIR/$RUN_NAME/n${SAMPLE_SIZE}_d${DIST_MB}mb/${MODEL}"
DEMOG_OUTDIR="$OUTDIR/demography"
DFE_OUTDIR="$OUTDIR/dfe"

echo "Task $SGE_TASK_ID -> $RUN_NAME n=$SAMPLE_SIZE d=${DIST_MB}mb model=$MODEL fold_mode=$FOLD_MODE"

if [ -s "$DFE_OUTDIR/dfe_results.csv" ]; then
  echo "Already complete ($DFE_OUTDIR/dfe_results.csv) -- skipping."
  exit 0
fi

syn_file="$SLIM_OUTDIR/syn_sfs_${RUN_NAME}_n${SAMPLE_SIZE}_d${DIST_MB}mb.tsv"
nonsyn_file="$SLIM_OUTDIR/nonsyn_sfs_${RUN_NAME}_n${SAMPLE_SIZE}_d${DIST_MB}mb.tsv"
bp_file="$SLIM_OUTDIR/genic_bp_${RUN_NAME}.csv"
for f in "$syn_file" "$nonsyn_file" "$bp_file"; do
  if [ ! -s "$f" ]; then echo "ERROR: missing $f"; exit 1; fi
done

# ------------------------------------------------------------------
# Sum segments and convert to dadi's unfolded Spectrum format (index 0
# padded with 0; dadi masks the 0 and n entries).
# ------------------------------------------------------------------
mkdir -p "$DEMOG_OUTDIR" "$DFE_OUTDIR"
to_dadi_sfs() {
  awk -v n="$SAMPLE_SIZE" '
    NF == 0 { next }
    { for (i = 1; i <= NF; i++) sum[i] += $i; ncols = NF }
    END {
      if (ncols != n) { print "ERROR: expected " n " columns, got " ncols > "/dev/stderr"; exit 1 }
      printf "%d unfolded\n0", n + 1
      for (i = 1; i <= ncols; i++) printf " %s", sum[i]
      printf "\n"
    }' "$1" > "$2"
}
SYN_DADI="$OUTDIR/${RUN_NAME}_syn_n${SAMPLE_SIZE}_d${DIST_MB}mb.sfs"
NONSYN_DADI="$OUTDIR/${RUN_NAME}_nonsyn_n${SAMPLE_SIZE}_d${DIST_MB}mb.sfs"
to_dadi_sfs "$syn_file" "$SYN_DADI" || exit 1
to_dadi_sfs "$nonsyn_file" "$NONSYN_DADI" || exit 1

# ------------------------------------------------------------------
# L_syn and the nonsyn/syn mutational-target ratio for the +/-DIST_MB window,
# summed over segments from genic_bp_<run>.csv (segment,dist_mb,exon_bp,
# intron_bp). Only genic sequence mutates in v4, so
#   L_syn    = exon_bp * 0.302 + intron_bp
#   L_nonsyn = exon_bp * 0.698
# ------------------------------------------------------------------
N_SEG_BP=$(awk -F',' -v d="$DIST_MB" '$1 != "segment" && $2 == d' "$bp_file" | wc -l)
if [ "$N_SEG_BP" -ne 22 ]; then
  echo "ERROR: expected 22 segments at d=${DIST_MB}mb in $bp_file, found $N_SEG_BP"
  exit 1
fi
read L_SYN NS_S_RATIO < <(awk -F',' -v frac="$EXON_SYN_FRACTION" -v d="$DIST_MB" '
  $1 != "segment" && $2 == d { ex += $3; in_ += $4 }
  END { printf "%.0f %.6g\n", ex * frac + in_, ex * (1 - frac) / (ex * frac + in_) }' "$bp_file")
echo "L_syn=$L_SYN  ns_s_ratio=$NS_S_RATIO"

# ------------------------------------------------------------------
# Demography, then DFE conditioned on it.
# ------------------------------------------------------------------
DEMOG_FILE="$DEMOG_OUTDIR/${MODEL}_demography.txt"
if ! grep -q "Best fit parameters:" "$DEMOG_FILE" 2>/dev/null; then
  python "$FIT_DEMOG_SCRIPT" "$SYN_DADI" --model_type "$MODEL" --L_syn "$L_SYN" \
    "$FOLD_FLAG" --pts $PTS "$DEMOG_OUTDIR/"
else
  echo "Reusing completed demography fit $DEMOG_FILE"
fi
if ! grep -q "Best fit parameters:" "$DEMOG_FILE" 2>/dev/null; then
  echo "ERROR: demography fit failed -- see $DEMOG_OUTDIR/log.log"
  exit 1
fi

python "$FIT_DFE_SCRIPT" "$SYN_DADI" "$NONSYN_DADI" "$DEMOG_FILE" "$MODEL" "$DFE_OUTDIR/" \
  --L_syn "$L_SYN" "$FOLD_FLAG" \
  --likelihood poisson --ns_s_ratio "$NS_S_RATIO" \
  --pts $PTS --cpus "${NSLOTS:-1}"
if [ ! -s "$DFE_OUTDIR/dfe_results.csv" ]; then
  echo "ERROR: DFE fit failed -- see $DFE_OUTDIR/log.log"
  exit 1
fi

echo "=== Finished $RUN_NAME n=$SAMPLE_SIZE d=${DIST_MB}mb model=$MODEL (L_syn=$L_SYN) ==="
