# ReadExtractorSMK

Extract or filter reads from fastq files based on Kraken 2 classifications.

This pipeline was created to run on a cluster with the SLURM Workload Manager and Conda package manager. Therefore this README is tailored towards that kind of environment.

Currently this is tailored to work on the PDC cluster Dardel. You just need to update the SLURM details in `slurm/config.yaml` if you want to use this pipeline on another cluster.

You can set the pipeline to extract reads classified to a specific taxon/clade, or to filter those reads and output all other reads.

## Usage

### Before you run the pipeline

#### Snakemake

Importantly, this pipeline was made to work with snakemake v8+.

Here is a command for installing snakemake and the required slurm executor in its own conda environment:

`conda create -c conda-forge -c bioconda -n snakemake snakemake snakemake-executor-plugin-slurm`

Before running the pipeline, activate the environment:

`conda activate snakemake`

#### Taxonomy

You need to supply the NCBI-style taxonomy files `names.dmp` and `nodes.dmp` that were used to create the Kraken 2 database used for classification.
Either supply paths to your taxonomy files in the `config.yaml`, or copy or create links to them in the folder `supporting_files`.

#### Taxonomic IDs

You need to supply the taxonomic ID(s) of the taxa or taxons for which you want to extract reads for. This is done in a simple text file.
The pipeline is set up to look for taxonomic IDs in `supporting_files/tax_ids.txt`, so easiest for you is to just create that file and fill
it with tax IDs. The file needs to have 1 tax ID per line. Alternatively, update the `config.yaml` to point to another file.

#### Input files

1) You need to populate the `input/classifications` folder with the kraken 2 output classification files (symlinks are fine).
2) You need to populate the `input/fastq` folder with the FASTQ files from which you want to extract reads from (symlinks are fine).

#### Configs

1) Make sure to go over the pipeline `config.yaml` in the root directory to make sure it fits your setup.
2) Make sure to go over the SLURM `config.yaml` in the slurm directory to make sure it fits your setup.

### Run the pipeline

Paste this into the terminal and run (hint: use screen):

`snakemake --executor slurm --profile slurm --cores 1 --use-conda`
