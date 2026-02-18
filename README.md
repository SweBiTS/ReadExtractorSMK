<p align="center">
  <img src="assets/logo.png" alt="ReadExtractorSMK Logo" width="400">
</p>

# ReadExtractorSMK

Extract or filter reads from fastq files based on Kraken 2 classifications.

This pipeline was created to run on a cluster with the SLURM Workload Manager and Conda package manager, but can be run locally as well.

Currently this is tailored to work on the PDC cluster Dardel or HPC2N. You just need to update the SLURM details in `slurm/config.yaml` if you want to use this pipeline on another cluster.

You can set the pipeline to extract reads classified to a specific taxon/clade, or to filter those reads and output all other reads.

## Usage

### Before you run the pipeline

#### Snakemake

Importantly, this pipeline was made to work with **snakemake v8+**.

Here is a command for installing snakemake and the required slurm executor in its own conda environment:

`conda create -c conda-forge -c bioconda -n snakemake snakemake snakemake-executor-plugin-slurm`

Before running the pipeline, activate the environment:

`conda activate snakemake`

#### Taxonomy

You need to supply the NCBI-style taxonomy files `names.dmp` and `nodes.dmp` that were used to create the Kraken 2 database used for classification.
Either supply paths to your taxonomy files in the `config.yaml`, or copy or create links to them in the folder `supporting_files`.

#### Taxonomic IDs

Supply the taxonomic ID(s) for which you want to extract/filter reads. This is a simple text file with one ID per line.

The pipeline defaults to `supporting_files/tax_ids.txt`. Alternatively, update config.yaml to point to another file. The file needs to have 1 tax ID per line.

#### Input files

1) Populate `input/classifications` with the kraken 2 output classification files (symlinks are fine).
2) Populate `input/fastq` with FASTQ files from which you want to extract/filter reads from (symlinks are fine).

#### Configs

1) **Pipeline** `config.yaml`: ensure paths and parameters (like mode and include) are correct.
2) **SLURM** `slurm/config.yaml`: ensure project number and partition suits your slurm context.

### Run the pipeline

#### Cluster execution

To handle varying memory requirements across different samples (especially for abundant clades), use the --restart-times flag. This allows the pipeline to automatically resubmit jobs with doubled memory if they hit a SLURM memory limit.

Run the following command (hint: use screen or tmux):

`snakemake --executor slurm --profile slurm --cores 1 --use-conda --restart-times 4`

#### Local execution

If you are running this on a standalone server instead of a cluster, you can bypass the SLURM executor. Snakemake will manage the cores and memory directly, just tell it how many cores and how much memory it has to work with.

To run locally, use the following command:

`snakemake --cores 64 --use-conda --restart-times 4 --resources mem_mb=1500000`
