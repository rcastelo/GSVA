#!/bin/bash
## peak memory of reading into main memory, as sparse, one block of columns of
## the processed data stored on disk, as GSVA's loader of sparse on-disk data
## does with each of its blocks, see .load_sparse_by_blocks(), per nonzero
## value of the block; each block is read in a new R process under cgmemtime
## (macmemtime on macOS), and the peak is taken from the memory samples
## between the start and the end of the reading, relative to the memory at its
## start, which excludes starting R and opening the data. it calibrates
## 'sparsereadbytes', the memory per nonzero value of a block of the CSC
## layout while it is read, which includes the block read, from the largest
## number of bytes per nonzero value of the peak of the process ('bytes_nz')
## over the widths of the blocks, and 'densereadfactor', the memory of a block
## of the dense layout while it is read relative to the size of its dense
## form, from the largest ratio of 'bytes_nz' to 'dense_nz', the size of the
## dense form per nonzero value, both with a margin, see README.md. the dense
## layout is read through the dense form of its chunks, so that its memory
## depends on the size of that form, whatever the fraction of nonzero values.
## 'result_nz' gives the memory of the block read per nonzero value, and
## 'rpeak_nz' the peak memory allocated by R. the narrowest blocks
## take a fixed memory, such as the buffers of the HDF5 library, which adds
## more per nonzero value and matters less in the blocks of GSVA's loader,
## which are at least as wide as the default block of the DelayedArray
## package. run it from the calibration directory, one run at a time,
## e.g.:
##   ./measure_sparse_read_memory.sh <label> csc,dense 1000,5000,20000,50000 3
## where <label> names the GSVA library used, which is the one of R_LIBS, and
## widths of a few thousand columns or more give blocks whose reading lasts
## long enough to be sampled. the results are written to
## results/sparse_read_memory_<label>.tsv

usage() {
  echo "Usage: $0 <label> [layouts] [ncols] [reps] [limit]"
  echo "  layouts  comma-separated layouts of the data on disk, csc and/or"
  echo "           dense (default: csc,dense)"
  echo "  ncols    comma-separated numbers of columns of the block read"
  echo "           (default: 1000,5000,20000,50000)"
  echo "  reps     number of runs of each layout and width (default: 3)"
  echo "  limit    memory limit of each run, above which it is killed, as"
  echo "           a workload manager would do (default: 64G)"
}

if [[ $# -lt 1 || $1 == "-h" ]]; then
  usage
  exit 1
fi
label=$1
layouts=${2:-csc,dense}
ncolss=${3:-1000,5000,20000,50000}
reps=${4:-3}
limit=${5:-64G}

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

out=results/sparse_read_memory_${label}.tsv
mkdir -p results
echo "host $(hostname), cores $(nproc 2>/dev/null || sysctl -n hw.ncpu), R_LIBS=${R_LIBS}"

printf "label\tlayout\tncols\trep\tnnz\tpeak_anon\tbytes_nz\tresult_nz\tdense_nz\trpeak_nz\telapsed\tstatus\n" > $out
first=yes
for layout in $(echo $layouts | tr ',' ' '); do
  for ncols in $(echo $ncolss | tr ',' ' '); do
    for rep in $(seq 1 $reps); do
      log=$(mktemp "${TMPDIR:-/tmp}/sparseread.log.XXXXXX")
      samples=$(mktemp "${TMPDIR:-/tmp}/sparseread.samples.XXXXXX")
      $memtime -t -d ',' -i 20 -s $samples -M $limit \
        Rscript $here/measure_sparse_read_memory.R $layout $ncols > $log 2>&1
      [[ $first == yes ]] && grep '^GSVA' $log
      first=no
      start=$(awk '$1 == "PHASE" && $2 == "start" {print $3}' $log)
      end=$(awk '$1 == "PHASE" && $2 == "end" {print $3}' $log)
      nnz=$(awk '$1 == "NNZ" {print $2}' $log)
      result=$(awk '$1 == "RESULT" {print $2}' $log)
      dense=$(awk '$1 == "DENSE" {print $2}' $log)
      rpeak=$(awk '$1 == "RPEAK" {print $2}' $log)
      elapsed=$(awk '$1 == "ELAPSED" {print $2}' $log)
      ## a run that did not end was either killed above the limit or ended
      ## with an error
      status=ok
      if [[ -z $end ]]; then
        status=killed
        grep -q '^Error' $log && status=failed
      fi
      ## the increase of the peak of the samples from the start of the
      ## reading to its end, or to the last sample when it was killed, over
      ## the last sample before it starts, in GiB, and per nonzero value, in
      ## bytes
      peak=
      [[ -n $start ]] && peak=$(awk -F, -v s="$start" -v e="${end:-1e12}" '
        NR == 1 || $1 <= s {a0 = $3}
        $1 >= s && $1 <= e {if ($3 > a) a = $3}
        END {d = a - a0; if (d < 0) d = 0; printf "%.0f", d * 1024}' $samples)
      perz() {
        awk -v x="$1" -v n="$nnz" 'BEGIN {if (x == "" || n == "" || n == 0)
                                              print "NA"
                                          else printf "%.1f", x / n}'
      }
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $label $layout \
        $ncols $rep "${nnz:-NA}" \
        $(awk -v p="$peak" 'BEGIN {if (p == "") print "NA"
                                    else printf "%.3f", p / 1024^3}') \
        $(perz "$peak") $(perz "$result") $(perz "$dense") $(perz "$rpeak") "${elapsed:-NA}" \
        $status | tee -a $out
      [[ $status != ok ]] && tail -5 $log
      rm -f $log $samples
    done
  done
done

## the largest number of bytes per nonzero value of each layout and width
echo
echo "largest bytes per nonzero value of the peak of the process, and its"
echo "ratio to the size of the dense form of the block"
printf "layout\tncols\tbytes_nz\tdense_nz\tdense_ratio\n"
awk -F'\t' 'NR > 1 && $7 != "NA" {k = $2 "\t" $3
    if (!(k in m) || $7 + 0 > m[k]) {m[k] = $7 + 0; d[k] = $9}}
  END {for (k in m) printf "%s\t%.1f\t%s\t%.2f\n", k, m[k], d[k],
                           m[k] / d[k]}' $out |
  sort -k1,1 -k2,2n
