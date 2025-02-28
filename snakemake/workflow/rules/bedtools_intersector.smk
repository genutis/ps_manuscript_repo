# 12/06/21
# updated 1/22/25
# this will just run the multicov for bam coverage counts now
# this takes the intersects of the tes and gtf file and uses bedtools to interlay them
# also computes gc content
# not used in pipeline
rule bedtools_intersect:
    input:
        ps = "{out}/tmp/PS_CisS_PnM_Abed_Erceg_Golobordko2019.{genome}.sorted.bed",
        dnabam = "{out}/data/bam/JJg14_PnM_S1_L002.{genome}.sorted.rg.bam",
        rna1bam = "{out}/data/bam/SRR8058303.{genome}.sorted.rg.bam",
        rna2bam = "{out}/data/bam/SRR8058304.{genome}.sorted.rg.bam",
        rna3bam = "{out}/data/bam/SRR8058305.{genome}.sorted.rg.bam",
        rna4bam = "{out}/data/bam/SRR8058306.{genome}.sorted.rg.bam",
        pi = "{out}/tmp/JJg14_S1_L002.{genome}.windowed.sorted.pi",
        sizes = "resources/{genome}.chrom.sizes",
        gtf = "{out}/tmp/{genome}.refGene.sorted.gtf",
        pnm_bed = "{out}/data/tes/JJg14_PnM_S1_L002/JJg14_PnM_S1_L002.{genome}.TEinsertions.sorted.bed",
        mat_bed = "{out}/data/tes/JJg21_DGRP057mat_S1_L001/JJg21_DGRP057mat_S1_L001.{genome}.TEinsertions.sorted.bed",
        pat_bed = "{out}/data/tes/JJg22_DGRP439pat_S2_L001/JJg22_DGRP439pat_S2_L001.{genome}.TEinsertions.sorted.bed",
        rna1_bed = "{out}/data/tes/SRR8058303/SRR8058303.{genome}.TEinsertions.sorted.bed",
        rna2_bed = "{out}/data/tes/SRR8058304/SRR8058304.{genome}.TEinsertions.sorted.bed",
        rna3_bed = "{out}/data/tes/SRR8058305/SRR8058305.{genome}.TEinsertions.sorted.bed",
        rna4_bed = "{out}/data/tes/SRR8058306/SRR8058306.{genome}.TEinsertions.sorted.bed",
        fasta = "resources/genomes/{genome}.fasta",
        fai = "resources/genomes/{genome}.fasta.fai"
    params:
        windows = "4000"
    output:
        int1 = temp("{out}/tmp/beds/intersect1.{genome, [A-Za-z0-9]+}.bed"),
        int2 = temp("{out}/tmp/beds/intersect2.{genome, [A-Za-z0-9]+}.bed"),
        int3 = temp("{out}/tmp/beds/intersect3.{genome, [A-Za-z0-9]+}.bed"),
        int4 = temp("{out}/tmp/beds/intersect4.{genome, [A-Za-z0-9]+}.bed"),
        int5 = temp("{out}/tmp/beds/intersect5.{genome, [A-Za-z0-9]+}.bed"),
        int6 = temp("{out}/tmp/beds/intersect6.{genome, [A-Za-z0-9]+}.bed"),
        int7 = temp("{out}/tmp/beds/intersect7.{genome, [A-Za-z0-9]+}.bed"),
        int8 = temp("{out}/tmp/beds/intersect8.{genome, [A-Za-z0-9]+}.bed"),
        int9 = temp("{out}/tmp/beds/intersect9.{genome, [A-Za-z0-9]+}.bed"),
        int10 = "{out}/tmp/beds/intersect10.{genome, [A-Za-z0-9]+}.bed",
    log:
        "{out}/logs/bedtools_intersect/{genome}.log"
    threads:
        config['bedtools_intersect']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['bedtools_intersect']['resources']['mem_mb'],
        time = config['bedtools_intersect']['resources']['time']
#    envmodules:
#        "gcc/8.3.0",
#        "bedtools2/2.27.1"
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        bedtools makewindows -g {input.sizes} -w {params.windows} |
#        bedtools nuc -fi {input.fasta} -bed "stdin" |
        bedtools multicov -bams {input.dnabam} {input.rna1bam} {input.rna2bam} {input.rna3bam} {input.rna4bam} -bed "stdin" |
        bedtools intersect -sorted -a "stdin" -b {input.ps} -wao > {output.int1}
        bedtools intersect -sorted -a {output.int1} -b {input.gtf} -wao > {output.int2}
        bedtools intersect -sorted -a {output.int2} -b {input.pi} -wao > {output.int3}
        bedtools intersect -sorted -a {output.int3} -b {input.pnm_bed} -wao > {output.int4}
        bedtools intersect -sorted -a {output.int4} -b {input.mat_bed} -wao > {output.int5}
        bedtools intersect -sorted -a {output.int5} -b {input.pat_bed} -wao > {output.int6}
        bedtools intersect -sorted -a {output.int6} -b {input.rna1_bed} -wao > {output.int7}
        bedtools intersect -sorted -a {output.int7} -b {input.rna2_bed} -wao > {output.int8}
        bedtools intersect -sorted -a {output.int8} -b {input.rna3_bed} -wao > {output.int9}
        bedtools intersect -sorted -a {output.int9} -b {input.rna4_bed} -wao > {output.int10}
        """
