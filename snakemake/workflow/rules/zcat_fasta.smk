#12/02/21

# star needs unzipped fasta so this makes one temporarily from reference fasta
rule zcat_fasta:
    input:
        fa="resources/genomes/{genome}.fasta.gz",
        gtf="resources/genomes/{genome}.refGene.gtf.gz"
    output:
        fa=temp("{out}/tmp/{genome}.fasta"),
        gtf=temp("{out}/tmp/{genome}.refGene.gtf")
    log:
        "{out}/logs/zcat_fasta/{genome}.log"
    benchmark:
        "{out}/logs/zcat_fasta/{genome}.benchmark.tsv"
    threads:
        config['zcat_fasta']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['zcat_fasta']['resources']['mem_mb'],
        time = config['zcat_fasta']['resources']['time']
    shell:
        """
        zcat {input.fa} > {output.fa}
        zcat {input.gtf} > {output.gtf}
        """
