# 1/14/21
# updated to use sra conda env on 10/11/23
# fastq dump on sra files
rule fastq_dump:
    input:
        "{out}/data/ncbi/public/sra/{sample}/{sample}.sra"
    output:
        "{out}/data/fastq/{sample, [A-Z0-9]+}_1.fastq.gz",
        "{out}/data/fastq/{sample, [A-Z0-9]+}_2.fastq.gz"
    log:
        "{out}/logs/fastq_dump/fastq_dump_{sample, [A-Z0-9]+}.log"
    benchmark:
        "{out}/logs/fastq_dump/{sample, [A-Z0-9]+}.benchmark.tsv"
    threads:
        config['fastq_dump']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['fastq_dump']['resources']['mem_mb'],
        time = config['fastq_dump']['resources']['time']
    conda:
        "envs/sra.yaml"
#    envmodules:
#        "gcc/8.3.0",
#        "intel/19.0.4",
#        "sra-toolkit/2.9.6"
    shell:
        """
        fastq-dump \
            --outdir {wildcards.out}/data/fastq/ \
            --split-3 \
            --gzip \
            --origfmt \
            {input}
        """

