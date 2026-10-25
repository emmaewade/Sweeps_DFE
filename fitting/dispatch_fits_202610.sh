#!/bin/bash
#$ -N dispatch202610
#$ -o /u/scratch/e/eewade/logs/dispatch202610.$TASK_ID.txt
#$ -j y
#$ -cwd
#$ -l h_rt=1:00:00,h_data=1G
#$ -t 1-5280:22
# One task per replicate (run_name). Submit with an array dependency on the
# simulation array, so task k (SGE_TASK_ID = 22(k-1)+1) only starts after
# simulation tasks SGE_TASK_ID..SGE_TASK_ID+21 -- that replicate's 22 segments,
# which are consecutive rows of single_task_args.txt -- have finished:
#   qsub -hold_jid_ad <sim_job_id> dispatch_fits_202610.sh
#
# Each task:
#   1. checks that all 22 segments wrote complete output in scratch (a segment that
#      hit h_rt or crashed still releases the hold, so this check is required);
#   2. copies the replicate's files to data/working-data/202610/tracksweep_output;
#   3. for human-DFE replicates, queues the dadi -> fitdadi fits for that replicate
#      (fit_demo_dfe_human_202610.sh, 12 tasks unfolded + 12 folded).
# Incomplete replicates are listed in dispatch_incomplete.txt and get no fits;
# fix/rerun their segments, then rerun this task by hand:
#   qsub -t <first_task_id> dispatch_fits_202610.sh

# batch jobs don't get UGE's bin/ on PATH; this provides qsub on the compute node
. /u/systems/UGE8.6.4/hoffman2/common/settings.sh

SWEEPS=/u/project/awbigham/eewade/daily/202610/sweeps
ARGS_FILE=$SWEEPS/single_task_args.txt
SLIM_OUTDIR=/u/scratch/e/eewade/tracksweep_output_202610
LOGDIR=/u/home/e/eewade/project-awbigham/data/working-data/202610
N_SEGMENTS=22
mkdir -p "$LOGDIR"

RUN_NAME=$(awk -v id="$SGE_TASK_ID" '$1 == id {print $4}' "$ARGS_FILE")
N_ROWS=$(awk -v a="$SGE_TASK_ID" -v b=$((SGE_TASK_ID + N_SEGMENTS - 1)) -v r="$RUN_NAME" \
  '$1 >= a && $1 <= b && $4 == r' "$ARGS_FILE" | wc -l)
if [ -z "$RUN_NAME" ] || [ "$N_ROWS" -ne "$N_SEGMENTS" ]; then
  echo "ERROR: task $SGE_TASK_ID is not the first of $N_SEGMENTS segments of one run_name ($RUN_NAME, $N_ROWS rows)"
  exit 1
fi
echo "Task $SGE_TASK_ID -> $RUN_NAME"

# ------------------------------------------------------------------
# Completeness: 22 segments in every per-run file.
# ------------------------------------------------------------------
problems=()
check_rows() {  # file, expected line count
  local n=0
  [ -e "$1" ] && n=$(wc -l < "$1")
  [ "$n" -eq "$2" ] || problems+=("$(basename "$1"): $n lines (expected $2)")
}
check_rows "$SLIM_OUTDIR/summary_${RUN_NAME}.csv" $((N_SEGMENTS + 1))
check_rows "$SLIM_OUTDIR/genic_bp_${RUN_NAME}.csv" $((N_SEGMENTS * 5 + 1))
for pre in "" pre_; do
  for kind in syn nonsyn; do
    for n in 10 100 1000; do
      for d in 1 2 3 4 5; do
        check_rows "$SLIM_OUTDIR/${pre}${kind}_sfs_${RUN_NAME}_n${n}_d${d}mb.tsv" $N_SEGMENTS
      done
    done
  done
done
[ -e "$SLIM_OUTDIR/hitchhikers_${RUN_NAME}.tsv" ] || problems+=("hitchhikers_${RUN_NAME}.tsv missing")
for seg in $(seq 1 $N_SEGMENTS); do
  for tag in pre post; do
    [ -s "$SLIM_OUTDIR/vcf_${RUN_NAME}_seg${seg}_${tag}.vcf.gz" ] || problems+=("vcf seg${seg} ${tag} missing")
  done
done

if [ "${#problems[@]}" -gt 0 ]; then
  echo "INCOMPLETE: ${#problems[@]} problem(s):"
  printf '  %s\n' "${problems[@]}"
  echo -e "$(date +%F_%T)\t$SGE_TASK_ID\t$RUN_NAME\t${problems[0]}" >> "$LOGDIR/dispatch_incomplete.txt"
  exit 0
fi

# ------------------------------------------------------------------
# Copy, then queue fits (human DFE only).
# ------------------------------------------------------------------
bash "$SWEEPS/copy_completed_sim_files.sh" "$RUN_NAME" || { echo "ERROR: copy failed"; exit 1; }

if [[ "$RUN_NAME" != *_dfehuman_* ]]; then
  echo "Not a human-DFE replicate -- copied only, no fits queued."
  echo -e "$(date +%F_%T)\t$SGE_TASK_ID\t$RUN_NAME\tcopied_only" >> "$LOGDIR/dispatch_submitted.txt"
  exit 0
fi

REP=${RUN_NAME##*_rep}
COND=$(echo "$RUN_NAME" | sed -E 's/^humansim_([a-z]+)_dfehuman_mu[^_]+_bens([^_]+)_rec([^_]+)_.*/\1_s\2_r\3/')
for FOLD_MODE in nofold fold; do
  JID=$(qsub -terse -t 1-12 -v RUN_NAME="$RUN_NAME",FOLD_MODE="$FOLD_MODE" \
    -N "dfe202610_${COND}_rep${REP}_${FOLD_MODE}" \
    -o "/u/scratch/e/eewade/logs/dadi202610_${RUN_NAME}_${FOLD_MODE}.\$TASK_ID.txt" \
    "$SWEEPS/fit_demo_dfe_human_202610.sh")
  if [ -z "$JID" ]; then echo "ERROR: qsub of $FOLD_MODE fits failed"; exit 1; fi
  echo "Queued $FOLD_MODE fits: job $JID"
  echo -e "$(date +%F_%T)\t$SGE_TASK_ID\t$RUN_NAME\t$FOLD_MODE\t$JID" >> "$LOGDIR/dispatch_submitted.txt"
done
