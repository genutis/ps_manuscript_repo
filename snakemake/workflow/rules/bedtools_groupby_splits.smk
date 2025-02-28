#!/bin/bash
# 12/06/21
# groups split bed files (basically removes duplicate rows and puts results in columns)
rule bedtools_groupby_splits:
    input:
        "{out}/tmp/beds/split/{genome}/{sample}"
    output:
        "{out}/tmp/beds/split/{genome}/grouped/{sample}.coordgrouped"
    log:
        "{out}/logs/bedtools_groupby_splits/{sample}.{genome}.log"
    threads:
        config['bedtools_groupby_splits']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['bedtools_groupby_splits']['resources']['mem_mb'],
        time = config['bedtools_groupby_splits']['resources']['time']
#    envmodules:
#        "gcc/8.3.0",
#        "bedtools2/2.27.1"
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        bedtools groupby \
        -c "$(seq -s "," 1 86)" \
        -o distinct \
        -i {input} > {output}
        """

