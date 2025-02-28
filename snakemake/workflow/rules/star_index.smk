# 12/02/21
# rule to create star index for alignment of rna
rule star_index:
    input:
        fa="{out}/tmp/{genome}.fasta",
        gtf="{out}/tmp/{genome}.refGene.gtf"
    output:
        directory("{out}/STAR/{genome}")
    log:
        "{out}/logs/star_index/{genome}.log"
    benchmark:
        "{out}/logs/star_index/{genome}.benchmark.tsv"
    threads:
        config['star_index']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['star_index']['resources']['mem_mb'],
        time = config['star_index']['resources']['time']
    conda:
        "envs/star.yaml"
    shell:
        """
            mkdir -p {wildcards.out}/STAR/{wildcards.genome} &&
            STAR \
            --runMode genomeGenerate \
            --genomeDir {wildcards.out}/STAR/{wildcards.genome} \
            --genomeFastaFiles {input.fa} \
            --sjdbGTFfile {input.gtf} \
            --genomeChrBinNbits 16 \
            --runThreadN {threads}
        """

