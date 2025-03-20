# ReadExtractorSMK
Extract reads from fastq files based on Kraken 2 classifications.

Currently this is tailored to work on **Dardel**. You just need to update the SLURM details in slurm/config.yaml if you want to use this pipeline somewhere else.

# Usage
## Before you run the pipeline
### Taxonomy
You need to supply the NCBI-style taxonomy files `names.dmp` and `nodes.dmp` that were used to create the Kraken 2 database used for classification.
Either supply paths to your taxonomy files in the `config.yaml`, or copy or create links to them in the folder `supporting_files`.

### Taxonomic IDs
You need to supply the taxonomic ID(s) of the taxa or taxons for which you want to extract reads for.
Either update the file `supporting_files/tax_ids.txt` or supply your own file. The file needs to have 1 tax ID per line.

### Input files
1) You need to populate the `input/classification` folder with the kraken 2 output classification files (symlinks are fine).
2) You need to populate the `input/fastq` folder with the FASTQ files from which you want to extract reads from (symlinks are fine).

### Configs
1) Make sure to go over the pipeline `config.yaml` in the root directory to make sure it fits your setup.
2) Make sure to go over the SLURM `config.yaml` in the slurm directory to make sure it fits your setup.

## Run the pipeline
Paste this into the terminal and run (hint: use screen):

`snakemake --executor slurm --profile slurm --cores 1 --use-conda`
