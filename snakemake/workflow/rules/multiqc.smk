# run multiqc on fastq files
rule multiqc:
    input:
        zip=expand("{{out}}/qc/fastqc/{sample}_{reads}_fastqc.zip", sample = samples['sample'], reads = [ 1, 2 ]),
        sstats=expand("{{out}}/qc/samtools_stats/{sample}.{{genome}}.txt", sample = samples['sample']),
    output:
        "{out}/qc/multiqc.{genome}.html"
    log:
        "{out}/logs/multiqc/multiqc.{genome}.log"
    benchmark:
        "{out}/logs/multiqc/benchmark.{genome}.tsv"
    threads:
        config['multiqc']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['multiqc']['resources']['mem_mb'],
        time = config['multiqc']['resources']['time']
    conda:
       "envs/multiqc.yaml"
    shell:
        """
        # force flag to force overwrites
        output_dir=$(dirname {output})
        output_name=$(basename {output})
        multiqc --force -o "$output_dir" -n "$output_name" {input}
        """



