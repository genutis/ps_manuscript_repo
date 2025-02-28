
# script to call up data and run R

rule rscript_wrapper_mouse_data:
    input:
        cfs_mat = "resources/GSE82185_ICM_rep123_maternal_40000_iced.matrix",
        cfs_pat = "resources/GSE82185_ICM_rep123_paternal_40000_iced.matrix",
        bins_mat = "resources/GSE82185_ICM_rep123_maternal_40000_abs.bed",
        bins_pat = "resources/GSE82185_ICM_rep123_paternal_40000_abs.bed",
        tes = "resources/mm9_rmsk_TE.gtf.locInd.locations",
    output:
        du_mat_mean_median = "{out}/data/du_mat_ps_means_medians.tsv",
        du_pat_mean_median= "{out}/data/du_pat_ps_means_medians.tsv",
        te_counts_mat = "{out}/data/te_counts_du_mat.tsv",
        te_counts_pat = "{out}/data/te_counts_du_pat.tsv",
    threads:
        config['rscript_wrapper_mouse_data']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['rscript_wrapper_mouse_data']['resources']['mem_mb'],
        time = config['rscript_wrapper_mouse_data']['resources']['time']
    conda:
        "envs/r.yaml"
    log:
        "{out}/logs/rscript_wrapper_mouse_data/log.log"
    benchmark:
        "{out}/logs/rscript_wrapper_mouse_data/benchmark.tsv"
    params:
        script = "workflow/scripts/pairing_score_from_matrix_mouse_output.R",
        cf_col = 3, #column specifying contact frequencies
    shell:
        """
        Rscript {params.script} {input.cfs_mat} {input.cfs_pat} {params.cf_col} {input.bins_mat} {input.bins_pat} {input.tes} {output.du_mat_mean_median} {output.du_pat_mean_median} {output.te_counts_mat} {output.te_counts_pat}

        """
