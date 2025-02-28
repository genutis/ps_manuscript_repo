# 05/22/22
# updated 1/22/24
# makes nuleotide statistic scaffold bed with ps file for later intersecting
#        ps = "{out}/tmp/PS_CisS_PnM_Abed_Erceg_Golobordko2019.dm3.sorted.bed",
#        sizes = "resources/{genome}.chrom.sizes",

rule bedtools_nuc:
    input:
        bed = "{out}/tmp/beds/intersects.{genome}.coordregrouped.bed",
        fasta = "resources/genomes/{genome}.fasta",
        fai = "resources/genomes/{genome}.fasta.fai",
#    params:
#        windows = "4000"
    output:
        temp("{out}/tmp/beds/int.nuc.{genome, [A-Za-z0-9]+}.bed")
    log:
        "{out}/logs/bedtools_nuc/{genome}.log"
    benchmark:
        "{out}/logs/bedtools_nuc/{genome}.benchmark.tsv"
    threads:
        config['bedtools_nuc']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['bedtools_nuc']['resources']['mem_mb'],
        time = config['bedtools_nuc']['resources']['time']
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        bedtools nuc -fi {input.fasta} -bed {input.bed} > {output}
        """
