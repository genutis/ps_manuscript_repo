# 12/6/21
# not used
# puts the split intersects back together again.
rule unsplit_beds:
    input:
        splits=expand("{{out}}/tmp/beds/split/{{genome}}/grouped/intersects.bed.{id}.coordgrouped", id = ["%.3d" % i for i in list(range(128))]),
        log="{out}/tmp/splitsdone.{genome}.log", # to make sure splits are there potentially
    output:
        unsplit="{out}/tmp/beds/intersects.{genome, [A-Za-z0-9]+}.coordgrouped.bed",
        touched=touch("{out}/tmp/unsplitsdone.{genome}.log"),
    log:
        "{out}/logs/unsplit_beds/{genome}.log"
    benchmark:
        "{out}/logs/prepare_header_files/{genome}.benchmark.tsv"
    threads:
        config['unsplit_beds']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['unsplit_beds']['resources']['mem_mb'],
        time = config['unsplit_beds']['resources']['time']
    shell:
        "cat {input.splits} > {output.unsplit}"

