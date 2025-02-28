# 12/06/21
# updated 1/22/25
# splits intersect file up into more manageable chunks
rule split_intersect_bed:
    input:
#        "{out}/tmp/beds/intersect.{genome}.bed"
        "{out}/tmp/beds/intersect10.{genome}.bed"
    output:
        touch("{out}/tmp/splitsdone.{genome}.log"), # this is to try and force all jobs to finish first before next, not sure if split stores in memory or writes in place...
        expand("{{out}}/tmp/beds/split/{{genome, [A-Za-z0-9]+}}/intersects.bed.{id}", id = ["%.3d" % i for i in list(range(128))]),
    log:
        "{out}/logs/split_intersect_bed/{genome}.log"
    benchmark:
        "{out}/logs/split_intersect_bed/{genome}.benchmark.tsv"
    threads:
        config['split_intersect_bed']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['split_intersect_bed']['resources']['mem_mb'],
        time = config['split_intersect_bed']['resources']['time']
    shell:
        """
        mkdir -p {wildcards.out}/tmp/beds/split/{wildcards.genome}
#        split --verbose -n l/128 -a 3 -d {input} {wildcards.out}/tmp/beds/split/{wildcards.genome}/intersects.bed.
        # note snakemake formatting for awk to escape curly bracket. "" are needed as well as {{}}
#        split --verbose -l $(expr $(awk "END {{printNR}}" {input}) /128) -a 3 -d {input} {wildcards.out}/tmp/beds/split/{wildcards.genome}/intersects.bed.
#        split --verbose -n 128 -a 3 -d {input} {wildcards.out}/tmp/beds/split/{wildcards.genome}/intersects.bed. casuing problems with uneaven fields
        split --verbose -l $(expr $(wc -l < {input}) / 128) -a 3 -d {input} {wildcards.out}/tmp/beds/split/{wildcards.genome}/intersects.bed.
        """


