# 11/30/21
# rule for building bwa index

rule bwa_index:
    input:
        "resources/genomes/{genome}.fasta.gz"
    output:
        "{out}/bwa/indices/{genome, [A-Za-z0-9]+}.amb",
        "{out}/bwa/indices/{genome, [A-Za-z0-9]+}.ann",
        "{out}/bwa/indices/{genome, [A-Za-z0-9]+}.bwt",
        "{out}/bwa/indices/{genome, [A-Za-z0-9]+}.pac",
        "{out}/bwa/indices/{genome, [A-Za-z0-9]+}.sa"
    log:
        "{out}/logs/bwa_index/{genome}.log"
    benchmark:
        "{out}/logs/bwa_index/{genome}.benchmark.tsv"
    threads:
        config['bwa_index']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['bwa_index']['resources']['mem_mb'],
        time = config['bwa_index']['resources']['time']
    conda:
        "envs/bwa.yaml"
    shell:
        """
        bwa index -p {wildcards.out}/bwa/indices/{wildcards.genome} {input}
        """
