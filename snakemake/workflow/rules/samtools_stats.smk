# 8/17/22
# run samtools_stats on bam files
rule samtools_stats:
    input:
        sample="{out}/data/bam/{sample}.{genome}.sorted.rg.bam",
        reference="{out}/tmp/{genome}.fasta"
    output:
        "{out}/qc/samtools_stats/{sample}.{genome}.txt",
    log:
        "{out}/logs/samtools_stats/samtools_stats_{sample}.{genome}.log"
    benchmark:
        "{out}/logs/samtools_stats/{sample}.{genome}.benchmark.tsv"
    threads:
        config['samtools_stats']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['samtools_stats']['resources']['mem_mb'],
        time = config['samtools_stats']['resources']['time']
    conda:
       "envs/samtools.yaml"
    shell:
        """
        samtools stats --threads {threads} --ref-seq {input.reference} {input.sample} > {output}
        """


