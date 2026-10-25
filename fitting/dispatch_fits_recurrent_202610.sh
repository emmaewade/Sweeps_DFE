#!/bin/bash
#$ -N dispatchrec202610
#$ -o /u/scratch/e/eewade/logs/dispatchrec202610.$TASK_ID.txt
#$ -j y
#$ -cwd
#$ -l h_rt=1:00:00,h_data=1G
#$ -t 1-3520:22
# Recurrent-sweep counterpart of dispatch_fits_202610.sh. One task per replicate
# (run_name); task k has SGE_TASK_ID = 22(k-1)+1, the first of that replicate's 22
# consecutive rows in recurrent_task_args.txt. Each task:
#   1. checks that all 22 segments wrote complete output in scratch;
#   2. copies the replicate's files to data/working-data/202610/recurrentsweep_output;
#   3. queues the dadi -> fitdadi fits (fit_demo_dfe_recurrent_202610.sh,
#      4 tasks unfolded + 4 folded: n = 100, 1000 x two/three-epoch, +/-5 Mb).
# Incomplete replicates go to dispatch_recurrent_incomplete.txt and get no fits;
# rerun all 22 of their segments (SFS rows carry no segment id, so a partial
# replicate can't be patched), then rerun this task:
#   qsub -hold_jid <rerun_job> -t <first_task_id> dispatch_fits_recurrent_202610.sh

# batch jobs don't get UGE's bin/ on PATH; this provides qsub on the compute node
. /u/systems/UGE8.6.4/hoffman2/common/settings.sh

SWEEPS=/u/project/awbigham/eewade/daily/202610/sweeps
ARGS_FILE=$SWEEPS/recurrent_task_args.txt
SLIM_OUTDIR=/u/scratch/e/eewade/recurrentsweep_output_202610
LOGDIR=/u/home/e/eewade/project-awbigham/data/working-data/202610
DEST_DIR=$LOGDIR/recurrentsweep_output
N_SEGMENTS=22
mkdir -p "$DEST_DIR"

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
check_rows "$SLIM_OUTDIR/genic_bp_${RUN_NAME}.csv" $((N_SEGMENTS + 1))
check_rows "$SLIM_OUTDIR/fixed_positions_${RUN_NAME}.txt" $N_SEGMENTS
for kind in syn nonsyn nonsyn_del nonsyn_ben; do
  for n in 10 100 1000; do
    check_rows "$SLIM_OUTDIR/${kind}_sfs_${RUN_NAME}_n${n}_d5mb.tsv" $N_SEGMENTS
  done
done
for seg in $(seq 1 $N_SEGMENTS); do
  [ -s "$SLIM_OUTDIR/vcf_${RUN_NAME}_seg${seg}.vcf.gz" ] || problems+=("vcf seg${seg} missing")
done

if [ "${#problems[@]}" -gt 0 ]; then
  echo "INCOMPLETE: ${#problems[@]} problem(s):"
  printf '  %s\n' "${problems[@]}"
  echo -e "$(date +%F_%T)\t$SGE_TASK_ID\t$RUN_NAME\t${problems[0]}" >> "$LOGDIR/dispatch_recurrent_incomplete.txt"
  exit 0
fi

# ------------------------------------------------------------------
# Copy, then queue fits.
# ------------------------------------------------------------------
for f in "$SLIM_OUTDIR"/*_"${RUN_NAME}".* "$SLIM_OUTDIR"/*_"${RUN_NAME}"_n*_d5mb.tsv \
         "$SLIM_OUTDIR"/vcf_"${RUN_NAME}"_seg*.vcf.gz; do
  [ -e "$f" ] && { cp "$f" "$DEST_DIR"/ || { echo "ERROR: copy failed: $f"; exit 1; }; }
done

REP=${RUN_NAME##*_rep}
COND=$(echo "$RUN_NAME" | sed -E 's/^recurrentsweep_dfehuman_mu[^_]+_bens([^_]+)_benfrac([^_]+)_rec([^_]+)_.*/s\1_p\2_r\3/')
for FOLD_MODE in nofold fold; do
  JID=$(qsub -terse -t 1-4 -v RUN_NAME="$RUN_NAME",FOLD_MODE="$FOLD_MODE" \
    -N "dferec202610_${COND}_rep${REP}_${FOLD_MODE}" \
    -o "/u/scratch/e/eewade/logs/dadirec202610_${RUN_NAME}_${FOLD_MODE}.\$TASK_ID.txt" \
    "$SWEEPS/fit_demo_dfe_recurrent_202610.sh")
  if [ -z "$JID" ]; then echo "ERROR: qsub of $FOLD_MODE fits failed"; exit 1; fi
  echo "Queued $FOLD_MODE fits: job $JID"
  echo -e "$(date +%F_%T)\t$SGE_TASK_ID\t$RUN_NAME\t$FOLD_MODE\t$JID" >> "$LOGDIR/dispatch_recurrent_submitted.txt"
done
