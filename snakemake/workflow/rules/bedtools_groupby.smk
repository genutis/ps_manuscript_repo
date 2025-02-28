# 12/06/21
# updated 1/22/25
# when combining split bed files this leaves an extra row sometimes so this fixes that.
rule bedtools_groupby:
    input:
        bed="{out}/tmp/beds/intersects.{genome}.coordgrouped.bed",
        log="{out}/tmp/unsplitsdone.{genome}.log",
    output:
        "{out}/tmp/beds/intersects.{genome}.coordregrouped.bed"
    log:
        "{out}/logs/bedtools_groupby/{genome}.log"
    threads:
        config['bedtools_groupby']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['bedtools_groupby']['resources']['mem_mb'],
        time = config['bedtools_groupby']['resources']['time']
#    envmodules:
#        "gcc/8.3.0",
#        "bedtools2/2.27.1"
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        bedtools groupby \
        -c "$(seq -s "," 1 $(awk '{{print NF; exit}}' {input.bed}))" \
        -o distinct \
        -i {input.bed} > {output}
        """

