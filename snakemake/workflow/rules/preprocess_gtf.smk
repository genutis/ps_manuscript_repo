#12/2/21
# create list of TEs from TE GTF for TEfinder if its not already generated

rule preprocess_gtf:
    input:
        "resources/genomes/{genome}_rmsk_TE.gtf"
    output:
        "{out}/tmp/list_of_tes_{genome}.txt"
    threads:
        config['preprocess_gtf']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['preprocess_gtf']['resources']['mem_mb'],
        time = config['preprocess_gtf']['resources']['time']
    shell:
        """
        awk -F '\t' '{{print $9}}' {input} | awk -F '"' '{{print $2}}' | \
        sort | uniq > {output}
        """
