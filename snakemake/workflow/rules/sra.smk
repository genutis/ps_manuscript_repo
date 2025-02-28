# runs sra prefetch on sample wildcard (expects an SRR number here), if there is no sample_1.fastq sample_2.fastq upstream (as fastq_dump.smk will need this job as a requirement if there is no fastqs for the sample id).
# this can pull single ended data but the pipeline is not set to anticipate that sort of output (as it will lack _1.fastq and _2.fastq files)
rule sra:
    output:
        "{out}/data/ncbi/public/sra/{sample, [A-Z0-9]+}/{sample}.sra"
    log:
        "{out}/logs/sra/sra_{sample}.log"
    conda:
        "envs/sra_3.1.1.yaml"
    benchmark:
        "{out}/logs/sra/sra_{sample}.benchmark.tsv"
    threads:
        config['sra']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['sra']['resources']['mem_mb'],
        time = config['sra']['resources']['time']
    params:
        outpath = "{out}/data/ncbi/public/sra",
    shell:
        """
        prefetch -X 9999999999999 {wildcards.sample} -O {params.outpath}
        """

