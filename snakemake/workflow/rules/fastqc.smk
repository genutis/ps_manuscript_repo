# run fastqc on fastq files
rule fastqc:
    input:
       "{out}/data/fastq/{sample}_1.fastq.gz",
       "{out}/data/fastq/{sample}_2.fastq.gz"
    output:
        "{out}/qc/fastqc/{sample}_1_fastqc.html",
        "{out}/qc/fastqc/{sample}_2_fastqc.html",
        "{out}/qc/fastqc/{sample}_1_fastqc.zip",
        "{out}/qc/fastqc/{sample}_2_fastqc.zip",
    log:
        "{out}/logs/fastqc/fastqc_{sample, [A-Z0-9]+}.log"
    benchmark:
        "{out}/logs/fastqc/{sample, [A-Z0-9]+}.benchmark.tsv"
    threads:
        config['fastqc']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['fastqc']['resources']['mem_mb'],
        time = config['fastqc']['resources']['time']
    conda:
       "envs/fastqc.yaml"
    shell:
        """
        mkdir -p "{wildcards.out}/qc/fastqc" &&
        fastqc -t {threads} -o "{wildcards.out}/qc/fastqc" {input}
        """


