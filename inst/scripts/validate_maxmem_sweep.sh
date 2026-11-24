#!/bin/bash
## peak memory and time of the three steps of GSVA for each number of cells,
## from <first> to <last> by <incr>, and each memory budget <maxmem>, each run
## in a new R process, see validate_maxmem_sweep.R, under cgmemtime (macmemtime
## on macOS), which records the peak memory of the process and its workers,
## and of each step, from its memory samples. it calibrates 'memfractionR' from
## the peak of the memory allocated by R against the peak of all the memory,
## and validates the memory model, whose runs under a hard memory limit equal
## to their budget (-l) must all complete, see README.md. run it from the
## calibration directory, with the CSV output redirected to a file, e.g.:
##   ./validate_maxmem_sweep.sh -l data/TENxBrainDataMouseE18_processed_csc.sce \
##       data/m8.all.v2025.1.Mm.symbols.gmt 1000000 1 1000000 yes \
##       16G,32G,64G,128G 10 > results/maxmem_sweep_1M.csv
## the environment variable GSVA_MEMTIME can give another memory profiler, or
## 'none' to run without one, e.g., in the jobs of a workload manager that
## does not let cgmemtime create its control group, whose memory limit is
## then the hard one, see validate_maxmem_sweep_slurm.sh; the columns of the
## memory of the process group are then NA

usage() {
  echo "Usage: $0 [-l] <se_dir> <gmt_file> <first> <incr> <last> <ondisk:auto|yes|no> <maxmem:auto|XG[,YG,...]> [ncores]"
  echo
  echo "  <maxmem> can be a comma-separated list of memory budgets, e.g., 16G,32G,64G,"
  echo "  in which case each number of cells is run with each of these budgets."
  echo "  -l  run each case under a hard memory limit equal to its budget <maxmem>,"
  echo "      as a job scheduler would do, recording in OOM_KILLS the number of"
  echo "      processes killed for exceeding it."
}

limit=no
while getopts "hl" opt; do
  case $opt in
    l) limit=yes ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done
shift $((OPTIND - 1))

if [[ $# -lt 7 ]]; then
  usage
  exit 1
fi

hdf5_se_dir=$1
gmt_file=$2
first=$3
incr=$4
last=$5
ondisk=$6
maxmem=$7

ncores=1
if [[ $# -ge 8 ]]; then
  ncores=$8
fi

## the R script next to this one, and the memory profiler, cgmemtime on Linux
## and macmemtime on macOS, built in the directory 'tools' next to them, see
## README.md
here=$(dirname "$0")
memtime=$here/tools/cgmemtime/cgmemtime
if [[ $(uname -s) == "Darwin" ]]; then
  memtime=$here/tools/macmemtime/macmemtime
fi
memtime=${GSVA_MEMTIME:-$memtime}
if [[ $memtime == "none" && $limit == "yes" ]]; then
  echo "a hard memory limit (-l) requires a memory profiler" 1>&2
  exit 1
fi
if [[ $memtime != "none" && ! -x $memtime ]]; then
  echo "$memtime not found, build it first, see README.md" 1>&2
  exit 1
fi

## memory budget in KiB, from a number of bytes or with suffix K, M, G or T
to_kib() {
  awk -v x="$1" 'BEGIN {
    n = x + 0; u = toupper(substr(x, length(x)))
    f = (u == "K") ? 1 : (u == "M") ? 1024 : (u == "G") ? 1024^2 : (u == "T") ? 1024^3 : 1/1024
    printf "%.0f", n * f
  }'
}

## number of cells and memory budgets; the format of the numbers of cells is
## given explicitly because the seq of macOS writes numbers of one million or
## more in exponential notation with six significant digits, e.g., 1234567 as
## 1.23457e+06, and can repeat the last number, e.g., 'seq 1000000 1 1000000'
## gives 1e+06 twice
ncells=($(seq -f %.0f $first $incr $last))
maxmems=(${maxmem//,/ })

if [[ $limit == "yes" ]]; then
  for mm in ${maxmems[@]}; do
    if [[ $mm == "auto" ]]; then
      echo "a hard memory limit (-l) requires explicit memory budgets, not 'auto'" 1>&2
      exit 1
    fi
  done
fi

## columns reported by the memory profiler
nmemcols=5
if [[ $limit == "yes" ]]; then
  nmemcols=6
fi

## columns from the progress of the phases written by the R script, see
## validate_maxmem_sweep.R, in the file 'phasefile', and from the samples of the memory
## of the process group written by the memory profiler in the file
## 'samplefile', as two groups separated by '|': the columns of the R script,
## which are 'NA' for the phases that did not complete, and the peak of the
## sampled memory of each phase, the peak of the sampled anonymous memory,
## which excludes the page cache, overall and of each phase, together with
## the last phase that began
phase_cols() {
  awk -F, -v pf="$1" -v mmkib="$2" '
    FILENAME == pf && $1 == "begin" { b[$2] = $3; order[++n] = $2; last = $2 }
    FILENAME == pf && $1 == "end"   { e[$2] = $3; used[$2] = $4; maxused[$2] = $5 }
    FILENAME == pf && $1 == "info"  { info[$2] = $3 }
    FILENAME != pf {
      ## the phase that began last before the sample, or the first one
      p = order[1]
      for (i = 2; i <= n; i++)
        if ($1 >= b[order[i]]) p = order[i]
      if (p != "" && (!(p in peak) || $2 + 0 > peak[p])) peak[p] = $2 + 0
      if (p != "" && (!(p in apeak) || $3 + 0 > apeak[p])) apeak[p] = $3 + 0
      if ($3 + 0 > ahigh) { ahigh = $3 + 0; anyanon = 1 }
    }
    function secs(p) { return (p in e) ? sprintf("%.3f", e[p] - b[p]) : "NA" }
    function val(x)  { return (x == "" || x == "NA") ? "NA" : x }
    function grp(p)  { return (p in peak) ? sprintf("%.0f", peak[p]) : "NA" }
    function anon(p) { return (p in apeak) ? sprintf("%.0f", apeak[p]) : "NA" }
    END {
      rhigh = "NA"
      if (("norm" in e) && ("ranks" in e) && ("scores" in e)) {
        rhigh = maxused["norm"]
        if (maxused["ranks"] + 0 > rhigh + 0) rhigh = maxused["ranks"]
        if (maxused["scores"] + 0 > rhigh + 0) rhigh = maxused["scores"]
      }
      if ("scores" in e) last = "done"
      if (last == "") last = "init"
      printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s|%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
        val(info["ngsets"]), ("maxmem" in info) ? info["maxmem"] : mmkib,
        secs("start"), secs("norm"), secs("ranks"), secs("scores"),
        val(used["start"]), val(maxused["norm"]), val(maxused["ranks"]),
        val(maxused["scores"]), rhigh,
        grp("start"), grp("norm"), grp("ranks"), grp("scores"),
        anyanon ? sprintf("%.0f", ahigh) : "NA",
        anon("start"), anon("norm"), anon("ranks"), anon("scores"), last
    }' "$1" "$3"
}

echo "NCELLS,NGSETS,NCORES,MAXMEM,ESTART,ENORM,ERANKS,ESCORES,R_MEM_START,R_MEM_NORM,R_MEM_RANKS,R_MEM_SCORES,R_MEM_HIGH,USER,SYS,WALL,CHILD_RSS_HIGH,GRP_MEM_HIGH,GRP_MEM_START,GRP_MEM_NORM,GRP_MEM_RANKS,GRP_MEM_SCORES,GRP_ANON_HIGH,GRP_ANON_START,GRP_ANON_NORM,GRP_ANON_RANKS,GRP_ANON_SCORES,OOM_KILLS,LAST_PHASE,STATUS,LOAD_BLOCKS,LOAD_AVAIL,LOAD_BLOCKMEM,LOAD_PROCMEM,LOAD_WORKERS"
for nc in ${ncells[@]}; do {
  for mm in ${maxmems[@]}; do {
    echo NCELLS: $nc MAXMEM: $mm 1>&2
    limitopt=""
    if [[ $limit == "yes" ]]; then
      limitopt="-M $mm"
    fi
    phasefile=$(mktemp "${TMPDIR:-/tmp}/maxmemsweep.phases.XXXXXX")
    samplefile=$(mktemp "${TMPDIR:-/tmp}/maxmemsweep.samples.XXXXXX")
    ## the memory profiler writes its columns at the end of the last line of
    ## its output, also when the R script does not complete
    if [[ $memtime == "none" ]]; then
      ## the columns of the memory profiler, as NA
      out=$(GSVA_PHASEFILE=$phasefile Rscript $here/validate_maxmem_sweep.R $hdf5_se_dir $gmt_file $nc $ondisk $mm $ncores)
      status=$?
      out="$(echo "$out" | tail -n 1)NA,NA,NA,NA,NA"
    else
      out=$(GSVA_PHASEFILE=$phasefile $memtime -t -d ',' -s $samplefile $limitopt Rscript $here/validate_maxmem_sweep.R $hdf5_se_dir $gmt_file $nc $ondisk $mm $ncores)
      status=$?
    fi
    out=$(echo "$out" | tail -n 1)
    memcols=$(echo "$out" | awk -F, -v n=5 -v m=$nmemcols '{ s = $(NF-m+1); for (i = NF-m+2; i <= NF-m+n; i++) s = s "," $i; print s }')
    oomkills=NA
    if [[ $limit == "yes" ]]; then
      oomkills=$(echo "$out" | awk -F, '{ print $NF }')
    fi
    mmkib=NA
    if [[ $mm != "auto" ]]; then
      mmkib=$(to_kib $mm)
    fi
    cols=$(phase_cols $phasefile $mmkib $samplefile)
    ## the reading of the input into main memory by parallel workers, NA
    ## when it is not read in blocks of several columns, see
    ## validate_maxmem_sweep.R
    loadcols=$(awk -F, '$1 == "info" { info[$2] = $3 }
      END {
        n = split("loadblocks loadavail loadblockmem loadprocmem loadworkers", k, " ")
        s = ""
        for (i = 1; i <= n; i++)
          s = s (i > 1 ? "," : "") ((k[i] in info) ? info[k[i]] : "NA")
        print s
      }' $phasefile)
    rm -f $phasefile $samplefile
    echo "$nc,${cols%%|*},$memcols,${cols##*|}" | awk -F, -v OFS=, -v nc=$ncores -v oom=$oomkills -v st=$status -v lc=$loadcols '{
      ## NCORES after NGSETS, OOM_KILLS before LAST_PHASE, STATUS after it,
      ## and the reading of the input last
      s = $1 "," $2 "," nc
      for (i = 3; i < NF; i++) s = s "," $i
      print s "," oom "," $NF "," st "," lc
    }'
  } done
} done
