# 05/24/22
# this takes the intersects of the tes and gtf file and uses bedtools to interlay them
# also computes gc content with bedops nuc
#        fasta = "resources/genomes/{genome}.fasta",
#        fai = "resources/genomes/{genome}.fasta.fai",
rule bedops_multi_intersect:
    input:
        nuc = "{out}/tmp/beds/nuc.{genome}.bed",
        ps = "{out}/tmp/PS_CisS_PnM_Abed_Erceg_Golobordko2019.{genome}.sorted.bed",
        pi = "{out}/tmp/JJg14_S1_L002.{genome}.windowed.sorted.pi",
        dna_bed = expand("{{out}}/data/tes/{sample}/{sample}.{{genome}}.TEinsertions.sorted.bed" , sample = samples.loc[samples['type'] == 'dna']['sample']),
        rna_bed = expand("{{out}}/data/tes/{sample}/{sample}.{{genome}}.TEinsertions.sorted.bed" , sample = samples.loc[samples['type'] == 'rna']['sample']),
        gtf = "{out}/tmp/{genome}.refGene.sorted.gtf",
        chip_bed = expand("{{out}}/data/peaks/{protein}.{{genome}}_peaks.sorted.xls", protein = proteins['chip1'])
    output:
        union="{out}/tmp/beds/union.{genome, [A-Za-z0-9]+}.bed",
        merged="{out}/tmp/beds/merged.{genome, [A-Za-z0-9]+}.bed",
        fillna="{out}/tmp/beds/intersect.nas.{genome, [A-Za-z0-9]+}.bed",
        nofillna="{out}/tmp/beds/intersect.nonas.{genome, [A-Za-z0-9]+}.bed",
    log:
        "{out}/logs/bedops_multi_intersect/{genome}.log"
    threads:
        config['bedops_multi_intersect']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['bedops_multi_intersect']['resources']['mem_mb'],
        time = config['bedops_multi_intersect']['resources']['time']
    conda:
        "envs/bedops.yaml"
    shell:
        """
        # take the union and their merged
        bedops -u {input} > {output.union}
        bedops -m {input} > {output.merged}

        # map to ids but don't fill nas here
        bedmap --echo --echo-map-id --delim '\\t' {output.merged} {output.union} > {output.nofillna}

        # fill with nas with this loop, this may be better downstream than the above formatting...
        for file in {input}; do
            echo "$file"
            bedops -n 1 {output.merged} "$file" |
            awk '{{ print $0"\\tNA" }}' |
            bedops -u - "$file" |
            cut -f4 > "$file".map
        done

        # combine merged files and maps
        paste {output.merged} {wildcards.out}/tmp/beds/*.bed.map > {output.fillna}
        """
