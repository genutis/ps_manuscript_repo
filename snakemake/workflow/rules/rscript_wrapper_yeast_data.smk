
# script to call up data and run R

rule rscript_wrapper_yeast_data:
    input:
        cfs = "resources/GSM2355767_Y12_DBVPG6044_exponential_Sau3AI.Y12xDBVPG6044.32000.matrix.txt",
        bins = "resources/GSE88952_Y12xDBVPG6044.32000.bed",
        tes = "resources/other_features_genomic_R64-2-1_20150113.fasta",
    output:
        DBVPG6044 = "{out}/data/DBVPG6044_ps_means_medians.csv",
        Y12 = "{out}/data/Y12_ps_means_medians.csv",
    benchmark:
        "{out}/logs/rscript_wrapper_yeast_data/benchmark.tsv"
    threads:
        config['rscript_wrapper_yeast_data']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['rscript_wrapper_yeast_data']['resources']['mem_mb'],
        time = config['rscript_wrapper_yeast_data']['resources']['time']
    conda:
        "envs/r.yaml"
    log:
        "{out}/logs/rscript_wrapper_yeast_data/log.log"
    params:
        script = "workflow/scripts/pairing_score_from_matrix_yeast_output.R",
        cf_col = 4, # column specifying contact frequencies
    shell:
        """
        Rscript {params.script} {input.cfs} {params.cf_col} {input.tes} {input.bins} {output.DBVPG6044} {output.Y12}

        """
