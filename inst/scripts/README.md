# Calibration of the memory model of GSVA

GSVA sizes the blocks of its calculations, and estimates the memory that they
require, so that they stay within the maximum main memory given by the
argument `maxmem`. This memory model relies on constants estimated from the
peak memory of the calculations measured on single-cell data, on macOS, on
Linux and on the nodes of a computer cluster with SLURM. They are defined in
`.memory_model_defaults`, in the file `R/memory_model.R` of the package, each
one with how it was measured.

The scripts in this directory measure those peaks again, to estimate the
constants anew when the hardware or the software, such as R itself or the
packages HDF5Array, SparseArray or DelayedArray, may have changed them. Any
constant can be replaced, without installing GSVA again, with the option
`GSVA.memorymodel`, a list with the same structure as `.memory_model_defaults`
or part of it, e.g.:

```r
options(GSVA.memorymodel=list(memfractionR=0.65,
                              stepfactors=list(colranks=list(sparsedisk=list(workfactor=15)))))
```

New constants should only be adopted after the sweep of memory budgets of
`validate_maxmem_sweep.sh`, run with a hard memory limit (`-l`), completes all
its runs.

## Setup

The scripts run from a calibration directory with a copy of this one, where
they write the data in `data` and the results in `results`:

```sh
cp -R "$(Rscript -e 'cat(system.file("scripts", package="GSVA"))')" gsvacalib
cd gsvacalib
mkdir -p results
make -C tools/cgmemtime cgmemtime   # on Linux
make -C tools/macmemtime            # on macOS
Rscript prepare_calibration_data.R --dense
```

`prepare_calibration_data.R` fetches the 1.3 million brain cells of the
TENxBrainData package through ExperimentHub, filters and normalizes them with
the scrapper package, and writes the log-counts in the compressed sparse column
(CSC) layout of 10x Genomics and, with `--dense`, also in a dense HDF5 dataset,
which takes about two hours and 26 GB of disk space more. It also downloads
the M8 collection of mouse cell type signatures of MSigDB.

The peak memory of a process and its parallel workers is measured with
`cgmemtime` on Linux, which needs control groups (cgroups) version 2, and with
`macmemtime` on macOS, see below. Runs should not share the machine with other
memory-intensive processes, and the scripts run one calibration run at a time.

## Constants and scripts

| Constant | Script | Estimation |
|---|---|---|
| `memfractionR` | `validate_maxmem_sweep.sh` | Ratio of the peak memory allocated by R (`R_MEM_*`, one core) to the peak memory of the process (`GRP_ANON_*`), with a margin |
| `procextramem` | `validate_maxmem_sweep.sh` | Memory outside R of each forked worker processing blocks read from disk (`ondisk=yes`): increase of the peak memory with the number of cores at the same budget, beyond the memory of the blocks |
| `stepfactors$*$forkcopy` | `measure_step_memory_forked.sh` | Increase of the peak memory with each forked worker at the same budget, with the input in main memory, beyond `procextramem`, relative to the memory of the main R process including its data |
| `stepfactors` with input in main memory | `measure_step_memory_forked.sh` | Peak memory of one worker at budgets just above GSVA's estimate, which must stay within them |
| `stepfactors$*$sparsedisk` | `measure_map_reduce_jobs.sh` | Peak memory of the jobs of `gsvaMap()` (`MaxRSS` of the accounting of SLURM) relative to the sparse size of their blocks |
| `maxparloaddensity` | `measure_input_load_time.R` | Fraction of nonzero values at which reading in parallel as sparse, whose time grows with the size of the sparse blocks, takes as long as reading as dense in one process |
| all | `validate_maxmem_sweep.sh -l` | Every run within a hard memory limit equal to its budget must complete |

`newprocmem`, `diskreadmem`, `sparsereadbytes`, `dgcassembly` and
`plagedelayedextra` were measured in ad hoc runs and have no script yet.

`plot_maxmem_sweep.R` plots the peak memory and the time of a sweep against
the memory budgets, and `plot_memory_time.R` against the number of cells.
`compare_map_reduce_jobs.R` compares the runs of `measure_map_reduce_jobs.sh`,
which must be adapted to the partitions of the cluster, see the header of each
script.

## Examples

```sh
# peaks of each step with 1 to 24 forked workers on 100,000 cells in memory
./measure_step_memory_forked.sh mylib 100000 1,4,12,24 120G \
    norm_ecdf,norm_clr,ranks,scores 8G,16G,32G,64G > results/step_memory_forked_mylib.out

# sweep of memory budgets on one million cells read from disk, with 10 cores
./validate_maxmem_sweep.sh -l data/TENxBrainDataMouseE18_processed_csc.sce \
    data/m8.all.v2025.1.Mm.symbols.gmt 1000000 1 1000000 yes 16G,32G,64G,128G 10 \
    > results/maxmem_sweep_1M.csv
Rscript plot_maxmem_sweep.R results/maxmem_sweep_1M.csv results/maxmem_sweep_1M.png

# time of loading the input in one process and with 6 workers
Rscript measure_input_load_time.R 100000 6 3
```

## Memory profilers

The directory `tools` holds the source code of the two memory profilers, which
report the peak memory of a command together with all its descendant
processes, such as the parallel workers of R, sample it to find the peak of
each step, and can kill the command when it exceeds a hard memory limit, as a
workload manager would do:

- `tools/cgmemtime` is a modified version of
  [cgmemtime](https://github.com/gsauthof/cgmemtime), © 2022 Georg Sauthoff,
  for Linux, which adds the hard memory limit (`-M`) and the memory samples
  (`-s`, `-i`), among other changes listed at the start of `cgmemtime.c`;
  `README.upstream.md` is the documentation of the original version.
- `tools/macmemtime` is a tool for macOS written for GSVA, with the interface
  and output of the modified `cgmemtime`, which samples the physical footprint
  of the process tree, since macOS has no cgroups.

Both are distributed under the GNU General Public License version 3 or any
later version, see the file `COPYING` in each directory, unlike the rest of
GSVA, see the file `COPYRIGHTS` of the package. They are separate programs,
which are neither built nor used by the package itself.
