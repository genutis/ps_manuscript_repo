#!/usr/bin/env bash
# this has to be set for process limit on login nodes for our hpc system.
export OPENBLAS_NUM_THREADS=1

# path of directory to the snakemake file to run
snakemake_dir=$(dirname "$0")
# log file dir setup
mkdir -p "$snakemake_dir"/results/logs/snakemake


# mem limit for runs on our shared hpc in mb:
# for endeavour the limit is 500gb
# for discovery i have it as 4tb

if [[ "$HOSTNAME" == *endeavour* ]]; then
    #mem_limit="262114" # 256gb
    mem_limit="786432" # 768gb
    profile="$snakemake_dir"/slurm.hpc.endeavour
else
    mem_limit="4200000" # ~4tb
    profile="$snakemake_dir"/slurm.hpc.discovery
fi
# max number of local jobs to run for local rules on the head node
localjobs="6"

# hpc profile in snakemake dir for use with our systems.
# loosely templated off of https://github.com/jdblischak/smk-simple-slurm

# activating snakemake conda env

source activate "$snakemake_dir"/resources/snakemake-5.32.2.yaml

date
echo -e "Using snakemake profile $profile\nMaximum total memory requested to slurm: $mem_limit mb\nStarting snakemake work in $snakemake_dir"

initial_dir=$(pwd)


# beginning snakemake work
cd "$snakemake_dir"

# unlocks directory if it gets locked from crash
snakemake --unlock

# first build a dag pdf to look at before anything is ran
snakemake --conda-frontend mamba --forceall --dag | dot -Tpdf > "$snakemake_dir"/results/dag.pdf

# now run snakemake
# --use-envmodules to use software installed in our module load system on hpc
#snakemake -np --profile "$profile" --use-envmodules --resources mem_mb="$mem_limit" localjobs="$localjobs"
# -n flag can now just be supplied along with any other command line argumetns to this command with $@
snakemake --conda-frontend mamba -p --profile "$profile" --resources mem_mb="$mem_limit" localjobs="$localjobs" "$@"

date
echo "Finished launching the snakemake jobs."

cd "$initial_dir"
