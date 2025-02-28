# 08/19/22
# adds bam coverages to the intersect bed file
# 1/31/25 updated chip equality to actually load chips...
rule bedtools_multicov:
    input:
#        bed = "{out}/tmp/beds/nuc.{genome}.bed",
        bed =  "{out}/tmp/beds/int.nuc.{genome}.bed",
        bams = expand("{{out}}/data/bam/{sample}.{{genome}}.sorted.rg.bam", sample = samples.loc[samples['type'] == 'chip']['sample']),
        bais = expand("{{out}}/data/bam/{sample}.{{genome}}.sorted.rg.bam.bai", sample = samples.loc[samples['type'] == 'chip']['sample'])
    output:
        "{out}/tmp/beds/multicov.{genome}.bed"
    log:
        "{out}/logs/bedtools_multicov/{genome}.log"
    benchmark:
        "{out}/logs/bedtools_multicov/{genome}.benchmark.tsv"
    threads:
        config['bedtools_multicov']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['bedtools_multicov']['resources']['mem_mb'],
        time = config['bedtools_multicov']['resources']['time']
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        bedtools multicov -bams {input.bams} -bed {input.bed} > {output}
  #| #apply headers next with sed
#         sed '1ichr\tstart\tstop\tpairing_score\tcis_score\tpercent_AT\tpercent_GC\tnumber_A\tnumber_C\tnumber_G\tnumber_N\tnumber_other\t'> {output}
        """
