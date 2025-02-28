# 8/24/22
# samtools index rule updated for conda
rule samtools_index:
    input:
        "{out}/data/bam/{sample}.{genome}.sorted.rg.bam"
    output:
        "{out}/data/bam/{sample, [A-Za-z0-9_]+}.{genome, [A-Za-z0-9]+}.sorted.rg.bam.bai"
    log:
        "{out}/logs/samtools_index/{sample}.{genome}.log"
    benchmark:
        "{out}/logs/samtools_index/{sample}.{genome}.benchmark.tsv"
    threads:
        config['samtools_index']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['samtools_index']['resources']['mem_mb'],
        time = config['samtools_index']['resources']['time']
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools index -@ {threads} {input} {output}
        """
