#!/bin/bash
## the sweep of validate_maxmem_sweep.sh on a computer cluster with SLURM,
## where each run, for each number of cells and each memory budget <maxmem>,
## is a job of its own whose memory (--mem) is that budget, which SLURM
## enforces as a hard memory limit, as validate_maxmem_sweep.sh does with its
## option -l, killing the job when it exceeds it, only when its control groups
## constrain the memory of the jobs, i.e., with ConstrainRAMSpace=yes in its
## cgroup.conf; otherwise jobs exceeding their memory complete, and the sweep
## does not validate the memory model. a job that allocates more than its
## memory tells whether it is enforced, e.g., with
##   sbatch --mem=1G --wrap "Rscript -e 'x <- numeric(2^28); x[] <- 1; Sys.sleep(30)'"
## which ends OUT_OF_MEMORY when it is, and otherwise COMPLETED. each job runs
## validate_maxmem_sweep.sh without a memory profiler, see GSVA_MEMTIME in
## validate_maxmem_sweep.sh, because cgmemtime cannot create its control group
## in the jobs of SLURM, and SLURM samples their peak memory every second, see
## --acctg-freq below. a last job, after all of them, joins their results in
## results/maxmem_sweep_<label>.csv
## and writes their accounting by SLURM, with their state, e.g.,
## OUT_OF_MEMORY, and their peak memory (MaxRSS), in
## results/maxmem_sweep_<label>_sacct.txt. with the memory limit enforced, the
## memory model is validated when all the jobs complete, see README.md, while
## MaxRSS, which adds up the memory that forked workers share with the main R
## process once for each of them, overestimates their memory. the GSVA library used is the one of
## R_LIBS, which the jobs inherit. run it from the calibration directory, on
## the login node, e.g.:
##   ./validate_maxmem_sweep_slurm.sh -p normal -t 02:00:00 mylib \
##       data/TENxBrainDataMouseE18_processed_dense.sce \
##       data/m8.all.v2025.1.Mm.symbols.gmt 100000 100000 300000 no \
##       8G,12G,16G 6

usage() {
  echo "Usage: $0 [-p partition] [-t time] [-a sbatch_args] <label> <se_dir> <gmt_file> <first> <incr> <last> <ondisk:auto|yes|no> <maxmem:XG[,YG,...]> [ncores]"
  echo
  echo "  -p  partition of the jobs (default: the one of SLURM)"
  echo "  -t  time limit of each job (default: 02:00:00)"
  echo "  -a  other arguments of sbatch for each job, e.g., \"--account=myaccount\""
  echo "  <label> names the GSVA library used, which is the one of R_LIBS, and the"
  echo "  results. the arguments after it are those of validate_maxmem_sweep.sh,"
  echo "  where each memory budget is also the memory of the job."
}

partition=""
timelimit=02:00:00
sbatchargs=""
while getopts "hp:t:a:" opt; do
  case $opt in
    p) partition="--partition=$OPTARG" ;;
    t) timelimit=$OPTARG ;;
    a) sbatchargs=$OPTARG ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done
shift $((OPTIND - 1))

if [[ $# -lt 8 ]]; then
  usage
  exit 1
fi
label=$1
se_dir=$2
gmt_file=$3
first=$4
incr=$5
last=$6
ondisk=$7
maxmem=$8
ncores=${9:-1}

## the scripts next to this one, by their absolute path, since the jobs run
## in the calibration directory
here=$(cd "$(dirname "$0")" && pwd)
if ! command -v sbatch > /dev/null; then
  echo "sbatch not found, run it on a computer cluster with SLURM" 1>&2
  exit 1
fi
for mm in ${maxmem//,/ }; do
  if [[ $mm == "auto" ]]; then
    echo "each job requires an explicit memory budget, not 'auto'" 1>&2
    exit 1
  fi
done

outdir=results/maxmem_sweep_${label}
mkdir -p $outdir
echo "label $label, R_LIBS=${R_LIBS}"

## the command of each job, without a memory profiler, while SLURM enforces
## the memory limit of the job, when it does
jobcmd() {
  echo "GSVA_MEMTIME=none" \
       "$here/validate_maxmem_sweep.sh $se_dir $gmt_file $1 1 $1 $ondisk $2 $ncores"
}

## the numbers of cells in the format given explicitly, see
## validate_maxmem_sweep.sh. the peak memory of each job (MaxRSS) is sampled
## every second, instead of the default of SLURM, often 30 seconds, which can
## miss the peaks of shorter phases
jobids=()
for nc in $(seq -f %.0f $first $incr $last); do
  for mm in ${maxmem//,/ }; do
    run=${nc}_${mm}
    id=$(sbatch --parsable $partition $sbatchargs --time=$timelimit \
                --job-name=maxmemsweep_${label}_${run} \
                --cpus-per-task=$ncores --mem=$mm --acctg-freq=task=1 \
                --export=ALL \
                --output=$outdir/run_${run}.csv --error=$outdir/run_${run}.err \
                --wrap "$(jobcmd $nc $mm)")
    if [[ -z $id ]]; then
      echo "the job of $run could not be submitted" 1>&2
      exit 1
    fi
    echo "submitted job $id: ncells $nc, maxmem $mm"
    jobids+=(${id%%;*})
  done
done

## the last job, after all the others, whatever their state, joins their
## results, with the header once, and writes their accounting by SLURM, whose
## peak memory MaxRSS is given for the batch step of each job
ids=$(IFS=,; echo "${jobids[*]}")
joincmd="sleep 30;
head -n 1 \$(ls $outdir/run_*.csv | head -n 1) > results/maxmem_sweep_${label}.csv;
for f in $outdir/run_*.csv; do tail -n +2 \$f; done >> results/maxmem_sweep_${label}.csv;
sacct -j $ids --format=JobID%20,JobName%50,State,ExitCode,Elapsed,ReqMem,MaxRSS \
  > results/maxmem_sweep_${label}_sacct.txt"
id=$(sbatch --parsable $partition $sbatchargs --time=00:10:00 \
            --job-name=maxmemsweep_${label}_join --mem=1G \
            --dependency=afterany:${ids//,/:} \
            --output=$outdir/join.out --error=$outdir/join.err \
            --wrap "$joincmd")
echo "submitted job ${id%%;*}, which joins the results in results/maxmem_sweep_${label}.csv"
