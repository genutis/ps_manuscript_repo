
# script to call up data and run R

rule rscript_wrapper_dmel_data:
    input:
        intersects = "{out}/tmp/beds/multicov.{genome}.bed",
    output:
        pdf = "{out}/figs/pairing_score_plots_horizontal.{genome, [A-Za-z0-9]+}.pdf",
    benchmark:
        "{out}/logs/rscript_wrapper/{genome}.benchmark.tsv"
    threads:
        config['rscript_wrapper']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config['rscript_wrapper']['resources']['mem_mb'],
        time = config['rscript_wrapper']['resources']['time']
    conda:
        "envs/r.yaml"
    log:
        "{out}/logs/rscript_wrapper/{genome}.log"
    benchmark:
        "{out}/logs/rscript_wrapper/{genome}.benchmark.tsv"
    params:
        script="workflow/scripts/pairing_score_from_intersects_dmel_output.R",
        figpath="{out}/figs",
        genome="{genome}",
    shell:
        """
        Rscript {params.script} {input.intersects} {params.figpath} {params.genome}
        """
