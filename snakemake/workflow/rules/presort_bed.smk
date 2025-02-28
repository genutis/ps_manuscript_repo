# 1/14/21
# this takes the intersects of the tes and gtf file and uses bedtools to interlay them
# note, rna stuff has to happen first
# one day i will make these not so hardcoded.
# windowed pi needs a rule to be added to this.
rule presort_bed:
    input:
        ps = "resources/PS_CisS_PnM_Abed_Erceg_Golobordko2019.{genome}.bed",
        pi = "resources/JJg14_S1_L002.{genome}.windowed.pi",
        gtf = "resources/genomes/{genome}.refGene.gtf.gz",
        pnm_bed = "{out}/data/tes/JJg14_PnM_S1_L002/JJg14_PnM_S1_L002.{genome}.TEinsertions.bed",
        mat_bed = "{out}/data/tes/JJg21_DGRP057mat_S1_L001/JJg21_DGRP057mat_S1_L001.{genome}.TEinsertions.bed",
        pat_bed = "{out}/data/tes/JJg22_DGRP439pat_S2_L001/JJg22_DGRP439pat_S2_L001.{genome}.TEinsertions.bed",
        rna1_bed = "{out}/data/tes/SRR8058303/SRR8058303.{genome}.TEinsertions.bed",
        rna2_bed = "{out}/data/tes/SRR8058304/SRR8058304.{genome}.TEinsertions.bed",
        rna3_bed = "{out}/data/tes/SRR8058305/SRR8058305.{genome}.TEinsertions.bed",
        rna4_bed = "{out}/data/tes/SRR8058306/SRR8058306.{genome}.TEinsertions.bed",

    output:
        pi = "{out}/tmp/JJg14_S1_L002.{genome, [A-Za-z0-9]+}.windowed.sorted.pi",
        ps = "{out}/tmp/PS_CisS_PnM_Abed_Erceg_Golobordko2019.{genome, [A-Za-z0-9]+}.sorted.bed",
        gtf = "{out}/tmp/{genome, [A-Za-z0-9]+}.refGene.sorted.gtf",
        pnm_bed = "{out}/data/tes/JJg14_PnM_S1_L002/JJg14_PnM_S1_L002.{genome, [A-Za-z0-9]+}.TEinsertions.sorted.bed",
        mat_bed = "{out}/data/tes/JJg21_DGRP057mat_S1_L001/JJg21_DGRP057mat_S1_L001.{genome, [A-Za-z0-9]+}.TEinsertions.sorted.bed",
        pat_bed = "{out}/data/tes/JJg22_DGRP439pat_S2_L001/JJg22_DGRP439pat_S2_L001.{genome, [A-Za-z0-9]+}.TEinsertions.sorted.bed",
        rna1_bed = "{out}/data/tes/SRR8058303/SRR8058303.{genome, [A-Za-z0-9]+}.TEinsertions.sorted.bed",
        rna2_bed = "{out}/data/tes/SRR8058304/SRR8058304.{genome, [A-Za-z0-9]+}.TEinsertions.sorted.bed",
        rna3_bed = "{out}/data/tes/SRR8058305/SRR8058305.{genome, [A-Za-z0-9]+}.TEinsertions.sorted.bed",
        rna4_bed = "{out}/data/tes/SRR8058306/SRR8058306.{genome, [A-Za-z0-9]+}.TEinsertions.sorted.bed"
    log:
        "{out}/logs/presort_bed/{genome}.log"
    benchmark:
        "{out}/logs/presort_bed/{genome}.benchmark.tsv"
    threads:
        config['presort_bed']['threads']

    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['presort_bed']['resources']['mem_mb'],
        time = config['presort_bed']['resources']['time']
    shell:
        """
#        sort -k1,1 -k2,2n {input.ps} > {output.ps}
#        sort -k1,1 -k2,2n {input.pi} > {output.pi}
#        sort -k1,1 -k2,2n {input.pnm_bed} > {output.pnm_bed}
#        sort -k1,1 -k2,2n {input.mat_bed} > {output.mat_bed}
#        sort -k1,1 -k2,2n {input.pat_bed} > {output.pat_bed}
#        sort -k1,1 -k2,2n {input.rna1_bed} > {output.rna1_bed}
#        sort -k1,1 -k2,2n {input.rna2_bed} > {output.rna2_bed}
#        sort -k1,1 -k2,2n {input.rna3_bed} > {output.rna3_bed}
#        sort -k1,1 -k2,2n {input.rna4_bed} > {output.rna4_bed}
#        sort -k1,1 -k4,4n <(zcat {input.gtf}) > {output.gtf} # lamda function could make these work in parallel perhaps except if file is gtf then it does this sort -k4,4n thing instead of the 2,2n
        # this is another approach below:
        sort -k1,1 -k2,2n {input.ps} > {output.ps} &
        sort -k1,1 -k2,2n {input.pi} > {output.pi} &
        sort -k1,1 -k2,2n {input.pnm_bed} > {output.pnm_bed} &
        sort -k1,1 -k2,2n {input.mat_bed} > {output.mat_bed} &
        sort -k1,1 -k2,2n {input.pat_bed} > {output.pat_bed} &
        sort -k1,1 -k2,2n {input.rna1_bed} > {output.rna1_bed} &
        sort -k1,1 -k2,2n {input.rna2_bed} > {output.rna2_bed} &
        sort -k1,1 -k2,2n {input.rna3_bed} > {output.rna3_bed} &
        sort -k1,1 -k2,2n {input.rna4_bed} > {output.rna4_bed} &
        sort -k1,1 -k4,4n <(zcat {input.gtf}) > {output.gtf}
        """
