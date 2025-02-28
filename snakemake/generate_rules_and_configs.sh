#!/bin/bash/

# 11/29/21

# script to pull my slurm scripts and convert to rules files with config information
# to be tweaked up by hand for each one used
# make snakemake directory structure
mkdir -p "$HOME"/scripts/hic/snakemake/{config,results/logs,resources,workflow/{rules,envs,scripts,report}}

find "$HOME"/scripts/hic/ -name "*.s*" -maxdepth 1 | while read -r file; do
    # pull script contents
    grep -v "SBATCH" "$file" > "$HOME"/scripts/hic/snakemake/workflow/rules/$(basename "$file" | cut -f 1 -d ".").smk
    # pull out sbatch parameters
    echo -e "\n$(basename "$file" | cut -f 1 -d ".")" >> "$HOME"/scripts/hic/snakemake/config/cluster.yaml
    grep "SBATCH" "$file" >> "$HOME"/scripts/hic/snakemake/config/cluster.yaml
done

# make nascent snakefile
find "$HOME"/scripts/hic/snakemake/workflow/rules/ -exec basename {} \; | awk '{FS="."} {print "rule " $1 ":"}' > "HOME"/scripts/hic/snakemake/snakefile
