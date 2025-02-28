# 1/14/21
# rule for mapping reads with bwa mem

rule bwa_mem:
    input:
        indx=rules.bwa_index.output,
        r1="{out}/data/fastq/{sample}_1.fastq.gz",
        r2="{out}/data/fastq/{sample}_2.fastq.gz"
    output:
         "{out}/data/bam/{sample}.{genome, [A-Za-z0-9]+}.sorted.rg.bam"
#         "{out}/data/bam/{sample}.{genome}.sorted.rg.bam"
    wildcard_constraints:
        sample = "|".join(samples.loc[samples['type'] != 'rna']['sample'])
    log:
        "{out}/logs/bwa_mem/{sample}.{genome}.log"
    benchmark:
        "{out}/logs/bwa_mem/{sample}.{genome}.benchmark.tsv"
    threads:
        config['bwa_mem']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['bwa_mem']['resources']['mem_mb'],
        time = config['bwa_mem']['resources']['time']
    conda:
        "envs/bwa.yaml"
    shell:
        """
        mkdir -p {wildcards.out}/tmp &&
        bwa mem -t {threads} \
        -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA\\tPU:HFKGFBCXY.2" \
        {wildcards.out}/bwa/indices/{wildcards.genome} \
        <(seqtk trimfq {input.r1}) \
        <(seqtk trimfq {input.r2}) |
        samtools sort -n \
        -@ {threads} \
        -T {wildcards.out}/tmp |
        samtools fixmate -m \
        -@ {threads} - - |
        samtools sort \
        -@ {threads} \
        -T {wildcards.out}/tmp |
        samtools markdup \
        -r \
        -@ {threads} \
        -T {wildcards.out}/tmp \
        --output-fmt BAM \
        - {output}
        """
