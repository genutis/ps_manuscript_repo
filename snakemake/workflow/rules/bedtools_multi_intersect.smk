# 02/03/22
# this takes the intersects of the tes and gtf file and uses bedtools to interlay them
# also computes gc content with bedtools nuc
rule bedtools_multi_intersect:
    input:
        ps = "{out}/tmp/PS_CisS_PnM_Abed_Erceg_Golobordko2019.{genome}.sorted.bed",
        pi = "{out}/tmp/JJg14_S1_L002.{genome}.windowed.sorted.pi",
        sizes = "resources/{genome}.chrom.sizes",
        gtf = "{out}/tmp/{genome}.refGene.sorted.gtf",
        dna_bed = expand("{{out}}/data/tes/{sample}/{sample}.{{genome}}.TEinsertions.sorted.bed" , sample = samples.loc[samples['type'] == 'dna']['sample']),
        rna_bed = expand("{{out}}/data/tes/{sample}/{sample}.{{genome}}.TEinsertions.sorted.bed" , sample = samples.loc[samples['type'] == 'rna']['sample']),
        fasta = "resources/genomes/{genome}.fasta",
        fai = "resources/genomes/{genome}.fasta.fai",
        chip_bed = expand("{{out}}/data/peaks/{protein}.{{genome}}_peaks.sorted.xls.noheader", protein = proteins['chip1'])
    params:
        windows = "4000"
    output:
        "{out}/tmp/beds/intersect.{genome, [A-Za-z0-9]+}.bed"
    log:
        "{out}/logs/bedtools_multi_intersect/{genome}.log"
    benchmark:
        "{out}/logs/bedtools_multi_intersect/{genome}.benchmark.tsv"
    threads:
        config['bedtools_multi_intersect']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['bedtools_multi_intersect']['resources']['mem_mb'],
        time = config['bedtools_multi_intersect']['resources']['time']
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        bedtools makewindows -g {input.sizes} -w {params.windows} |
        bedtools nuc -fi {input.fasta} -bed "stdin" |
        bedtools intersect -sorted -wa -wb -filenames  -a "stdin" -b \
        {input.ps} \
        {input.pi} \
        {input.dna_bed} \
        {input.rna_bed} \
        {input.chip_bed} > {output} # header might have been throwing off bedtools for these files
        """
