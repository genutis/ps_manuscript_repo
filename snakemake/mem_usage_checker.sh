#!/bin/bash/

# slurm logs to parse for snakemake, jobs should each have their own folder in this structure
snakedir="$HOME/scripts/hic/snakemake/results/logs"

for x in $(find "$snakedir" -maxdepth 1 -exec basename {} \; | tail -n +2)
    do
        echo "$x top 10 memory consuming jobs"
        for i in $(find "$snakedir"/"$x" -exec basename {} .out \;| cut -d "-" -f 2 | tail -n +2)
            do
                seff "$i" | grep "Memory Utilized" | sed "s/$/ $i/"
        done | sort -nr
done
