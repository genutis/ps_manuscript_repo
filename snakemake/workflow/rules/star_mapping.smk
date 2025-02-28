# 8/25/22
# rule to map rna reads with star

rule star_mapping:
    input:
        r1 = "{out}/data/fastq/{sample}_1.fastq.gz",
        r2 = "{out}/data/fastq/{sample}_2.fastq.gz",
        fa = "{out}/tmp/{genome}.fasta",
        gtf = "{out}/tmp/{genome}.refGene.gtf",
        indx = rules.star_index.output
    output:
        "{out}/data/bam/{sample}.{genome, [A-Za-z0-9]+}.sorted.rg.bam"
    wildcard_constraints:
        sample = "|".join(samples.loc[samples['type'] == 'rna']['sample'])
    log:
        "{out}/logs/star_mapping/star_{sample}.{genome}.log"
    benchmark:
        "{out}/logs/star_mapping/star_{sample}.{genome}.benchmark.tsv"
    threads:
        config['star_mapping']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['star_mapping']['resources']['mem_mb'],
        time = config['star_mapping']['resources']['time']
    conda:
        "envs/star.yaml"
    shell:
        """
        # remove temp files generated from any previous runs
        rm -rf {wildcards.out}/tmp/{wildcards.sample}.{wildcards.genome}

        # star mapping
        STAR \
        --runThreadN {threads} \
        --genomeDir {wildcards.out}/STAR/{wildcards.genome} \
        --readFilesCommand zcat \
        --sjdbGTFfile {input.gtf} \
        --readFilesIn {input.r1} {input.r2} \
        --outFileNamePrefix {wildcards.out}/data/bam/{wildcards.sample}.{wildcards.genome}.\
        --outTmpDir {wildcards.out}/tmp/{wildcards.sample}.{wildcards.genome} \
        --outFilterMismatchNoverLmax 0.05 \
        --outFilterMismatchNoverReadLmax 0.05 \
        --alignSJoverhangMin  8  \
        --outSAMtype BAM SortedByCoordinate \
        --alignSJDBoverhangMin 1  \
        --alignIntronMin 20 \
        --alignIntronMax 30000 \
        --alignMatesGapMax 30000 \
        --twopassMode None

        # rename output file for downstream work
        mv {wildcards.out}/data/bam/{wildcards.sample}.{wildcards.genome}.Aligned.sortedByCoord.out.bam \
        {wildcards.out}/data/bam/{wildcards.sample}.{wildcards.genome}.sorted.rg.bam
        """

