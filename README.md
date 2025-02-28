# ps_manuscript_repo
Git repo for final version of pairing score analysis manuscript.



## System requirements

- 128gb memory per node.
- Nodes must be internet connected for sra rule (although this can be configured to run on login nodes through appending the local rule in the Snakefile with sra).


## Software Requirements

- Slurm or another job scheduler
- Snakemake
- Conda


## Usage

Snakemake pipeline must be configured for your system's resource policies. This is set as an example for a slurm configuration. The pipeline can be configured and ran through the snakemake_launcher.sh wrapper executable script. Drosophila melanogaster pnm line component of pipeline requires PnM DNA sequencing fastq files, maternal DNA sequencing fastq files, and paternal DNA sequencing fastq files, providided by Dr. Jumana AlHaj Abed. PnM RNA sequencing files and other data are pulled from the National Center for Biotechnology Information Sequencing Read Archive (NCBI SRA). All other required software dependencies are provided or built via conda channels in the snakemake pipeline. A directed acyclic graph describing the pipeline is provided in the results directory. For the mouse pairing analysis, contact frequency matrices for maternal and paternal samples must be downloaded from NCBI Gene expression omnibus (GSE82185) as these files are too large to provide in this repo. Please place the files in resources/GSE82185_ICM_rep123_maternal_40000_iced.matrix and resources/GSE82185_ICM_rep123_paternal_40000_iced.matrix. These files can also be provided by request.


For a dry-run (omit flag for full run):

snakemake_launcher.sh --dry-run
