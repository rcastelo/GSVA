#!/bin/bash
## peak memory of each step of GSVA with forked workers, against GSVA's
## estimate, on the first <ncells> cells of the CSC data held in main memory;
## each step runs in a new R process under cgmemtime (macmemtime on macOS),
## and the peak is taken from the memory samples between the start and the
## end of the step, which excludes the loading of its input. it calibrates
## the constants of the memory model of forked workers with the input in main
## memory, 'procextramem' and the 'forkcopy' of each step, from the increase
## of the peak with each worker at the same memory budget, and validates the
## step factors with the input in main memory, see README.md. run it from the
## calibration directory on a machine with as many cores as the largest number
## of workers and more memory than the limit, one run at a time, e.g.:
##   ./measure_step_memory_forked.sh <label> 100000 1,4,12,24 120G \
##       norm_ecdf,norm_clr,ranks,scores 8G,16G,32G,64G
## where <label> names the GSVA library used, which is the one of R_LIBS. with
## memory budgets other than "auto", each run is killed above its budget. the
## results are written to results/step_memory_forked_<label>.tsv, with the
## memory in GiB

usage() {
  echo "Usage: $0 <label> [ncells] [workers] [limit] [steps]"
  echo "  ncells   number of cells (default: 100000)"
  echo "  workers  comma-separated numbers of workers (default: 1,4,16,32)"
  echo "  limit    memory limit of each run, above which it is killed, as"
  echo "           a workload manager would do (default: 120G)"
  echo "  steps    comma-separated steps (default: norm_ecdf,norm_clr,ranks,scores)"
  echo "  maxmem   comma-separated memory budgets given to GSVA, e.g. 8G,16G,32G;"
  echo "           a run with a budget other than auto is killed above it,"
  echo "           instead of above the limit (default: auto)"
}

if [[ $# -lt 1 || $1 == "-h" ]]; then
  usage
  exit 1
fi
label=$1
ncells=${2:-100000}
workers=${3:-1,4,16,32}
limit=${4:-120G}
steps=${5:-norm_ecdf,norm_clr,ranks,scores}
maxmems=${6:-auto}

## the R script next to this one, and the memory profiler, built in the
## directory 'tools' next to them, see README.md
here=$(dirname "$0")
memtime=$here/tools/cgmemtime/cgmemtime
if [[ $(uname -s) == "Darwin" ]]; then
  memtime=$here/tools/macmemtime/macmemtime
fi
if [[ ! -x $memtime ]]; then
  echo "$memtime not found, build it first, see README.md" 1>&2
  exit 1
fi

out=results/step_memory_forked_${label}.tsv
mkdir -p results
echo "host $(hostname), cores $(nproc 2>/dev/null || sysctl -n hw.ncpu), R_LIBS=${R_LIBS}"

## the input of the steps "ranks" and "scores"
if [[ ! -f data/forked_ranks_${ncells}.rds ]]; then
  echo "preparing the input of the steps ranks and scores"
  Rscript $here/measure_step_memory_forked.R prep $ncells $(echo $workers | tr ',' '\n' | sort -n | tail -1) || exit 1
fi

printf "label\tstep\tworkers\tmaxmem\testimate\tpermem\tpeak_mem\tpeak_anon\tstart_anon\trmem\telapsed\twarned\tkilled\n" > $out
for step in $(echo $steps | tr ',' ' '); do
  for w in $(echo $workers | tr ',' ' '); do
  for maxmem in $(echo $maxmems | tr ',' ' '); do
    runlimit=$limit
    [[ $maxmem != auto ]] && runlimit=$maxmem
    log=$(mktemp "${TMPDIR:-/tmp}/stepmem.log.XXXXXX")
    samples=$(mktemp "${TMPDIR:-/tmp}/stepmem.samples.XXXXXX")
    $memtime -t -d ',' -i 50 -s $samples -M $runlimit \
      Rscript $here/measure_step_memory_forked.R $step $ncells $w $maxmem > $log 2>&1
    head -1 $log | grep '^GSVA' > /dev/null && [[ $step == ${steps%%,*} && $w == ${workers%%,*} && $maxmem == ${maxmems%%,*} ]] && head -1 $log
    start=$(awk '$1 == "PHASE" && $2 == "start" {print $3}' $log)
    end=$(awk '$1 == "PHASE" && $2 == "end" {print $3}' $log)
    ## the largest estimate of the step and the memory per R process it uses
    est=$(awk '$1 == "ESTIMATE" && $2 + 0 > m {m = $2 + 0; p = $3} END {if (m) print m "\t" p; else print "NA\tNA"}' $log)
    rmem=$(awk '$1 == "RMEM" {print $2}' $log)
    elapsed=$(awk '$1 == "ELAPSED" {print $2}' $log)
    ## whether GSVA warned that the step may need more memory than its budget
    warned=no
    grep -q '^Warning' $log && warned=yes
    killed=no
    [[ -z $end ]] && killed=yes
    ## the peak of the samples from the start of the step to its end, or to
    ## the last sample when it was killed, in GiB
    peaks=$(awk -F, -v s="$start" -v e="${end:-1e12}" '
      NR == 1 || $1 <= s {a0 = $3}
      $1 >= s && $1 <= e {if ($2 > m) m = $2; if ($3 > a) a = $3}
      END {printf "%.3f\t%.3f\t%.3f", m / 1048576, a / 1048576, a0 / 1048576}' $samples)
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $label $step $w $maxmem \
      "$est" "$peaks" "${rmem:-NA}" "${elapsed:-NA}" $warned $killed | tee -a $out
    [[ $killed == yes ]] && tail -5 $log
    rm -f $log $samples
  done
  done
done
