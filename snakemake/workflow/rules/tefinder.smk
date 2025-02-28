# tefinder rule
rule tefinder:
    input:
        bam="{out}/data/bam/{sample}.{genome}.sorted.rg.bam",
        bai="{out}/data/bam/{sample}.{genome}.sorted.rg.bam.bai",
        ref="{out}/tmp/{genome}.fasta",
        gtf="resources/genomes/{genome}_rmsk_TE.gtf",
        tes="{out}/tmp/list_of_tes_{genome}.txt"
    output:
        bam="{out}/data/tes/{sample}/{sample}.{genome, [A-Za-z0-9]+}.DiscordantReads.bam",
        bed="{out}/data/tes/{sample}/{sample}.{genome, [A-Za-z0-9]+}.TEinsertions.bed"
    wildcard_constraints:
        sample = "|".join(samples.loc[samples['type'] != 'chip']['sample'])
    log:
        "{out}/logs/tefinder/{sample}.{genome}.log"
    benchmark:
        "{out}/logs/tefinder/{sample}.{genome}.benchmark.tsv"
    threads:
        config['tefinder']['threads']
    params:
        path = "workflow/scripts/TEfinder",
        mem_withheld_mb = 114688
    resources:
        mem_mb = config['tefinder']['resources']['mem_mb'],
        time = config['tefinder']['resources']['time'],
        partition = config['tefinder']['resources']['partition'],
    conda:
        "envs/tefinder.yaml"
    shell:
        """
        rm -rf {wildcards.out}/data/tes/{wildcards.sample}
        PICARD=$(find ".snakemake/conda/" -name "picard.jar") # this will find the picard installed from the tefinder yaml
        echo "picard found at: $PICARD"
        echo "java found at: $(which java)"
        {params.path} \
        -threads {threads} \
        -picard "$PICARD" \
        -maxHeapMem $(expr {resources.mem_mb} - {params.mem_withheld_mb}) \
        -alignment {input.bam} \
        -fa {input.ref} \
        -gtf {input.gtf} \
        -te {input.tes} \
        -workingdir {wildcards.out}/data/tes/{wildcards.sample} \
        -outname {wildcards.sample}.{wildcards.genome}.
        """
