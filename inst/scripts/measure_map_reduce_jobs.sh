#!/bin/bash
## driver job for measure_map_reduce_jobs.R: runs each configuration below,
## one after another, each of which submits its own jobs to SLURM through
## gsvaMap(), and writes after each run the accounting of the jobs started
## during that run, to find the time and memory of each job. submit it from
## the calibration directory, after adapting the partition below to the
## cluster, with:
##   sbatch measure_map_reduce_jobs.sh
## the driver only waits for the jobs of gsvaMap(), so it needs little memory
## and one CPU core, but a time limit long enough for all the configurations.
## its partition may need to differ from the one of the jobs of gsvaMap(),
## given by the environment variable GSVA_PARTITION, see
## measure_map_reduce_jobs.R, when that one is too short for the whole driver
#SBATCH --job-name=measure_map_reduce_jobs
#SBATCH --partition=normal
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=measure_map_reduce_jobs_%j.out
#SBATCH --error=measure_map_reduce_jobs_%j.err

## <layout> <ncells> <nodes> <ncpus> <mem>: 200,000 cells in chunks of 6 GB
## jobs with one or two CPU cores, on the dense and on the CSC layout, and
## with jobs of 16 GB, whose blocks of columns take whole chunks of the
## dense layout
configs=(
    "dense 200000 10 2 6G"
    "csc   200000 10 2 6G"
    "dense 200000 10 1 6G"
    "csc   200000 10 1 6G"
    "dense 200000 10 2 16G"
    "csc   200000 10 2 16G"
)

mkdir -p results
for cfg in "${configs[@]}"; do
    set -- $cfg
    tag="${1}_${2}_${3}x${4}x${5}"
    start=$(date +%Y-%m-%dT%H:%M:%S)
    Rscript measure_map_reduce_jobs.R "$@" > "results/map_reduce_jobs_${tag}.out" \
        2> "results/map_reduce_jobs_${tag}.err"
    status=$?
    end=$(date +%Y-%m-%dT%H:%M:%S)
    echo "$tag: exit status $status, from $start to $end"
    sleep 30 ## for the accounting of the last jobs to be recorded
    sacct -u "$USER" -S "$start" -E "$end" -X \
          --format=JobID,JobName%20,State,Submit,Start,End,Elapsed,AllocCPUS,ReqMem,MaxRSS \
          > "results/map_reduce_jobs_${tag}_sacct.txt"
    ## MaxRSS is only given for the steps of each job
    sacct -u "$USER" -S "$start" -E "$end" \
          --format=JobID,JobName%20,Elapsed,MaxRSS \
          > "results/map_reduce_jobs_${tag}_sacct_steps.txt"
done

Rscript compare_map_reduce_jobs.R
