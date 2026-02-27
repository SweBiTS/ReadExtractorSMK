# ReadExtractorSMK - Project Context

`ReadExtractorSMK` is a Snakemake workflow designed to extract or filter reads from FASTQ files based on Kraken 2 taxonomic classifications. It is optimized for use on SLURM-managed clusters but can also run locally.

## Project Overview

- **Purpose:** Efficiently isolate or remove reads corresponding to specific taxonomic IDs or entire clades from large sequencing datasets.
- **Workflow Engine:** Snakemake (v8+).
- **Key Components:**
    - `scripts/selmeout.py`: A Python script that parses Kraken 2 classification files and identifies read IDs belonging to target taxa or clades.
    - `rules/extract_reads.smk`: Snakemake rules defining the logic for read ID extraction and FASTQ filtering.
    - `SeqKit`: Used for high-performance filtering of FASTQ files based on read IDs.
    - `stringmeup`: A Python library used for handling NCBI taxonomy trees.

## Directory Structure

- `Snakefile`: The main entry point for the Snakemake workflow.
- `config.yaml`: Configuration for input/output paths, sample patterns, taxonomy file paths, and filtering modes.
- `rules/`: Contains the Snakemake rule definitions (`extract_reads.smk`).
- `scripts/`: Python scripts used by the workflow, notably `selmeout.py`.
- `envs/`: Conda environment definitions for reproducibility.
- `slurm/`: SLURM configuration for cluster execution, including a job script and resource profiles.
- `input/`: Directory for input data, organized into `classifications/` (Kraken 2 output) and `fastq/` (sequencing data).
- `supporting_files/`: Store for taxonomy files (`names.dmp`, `nodes.dmp`) and target taxonomic IDs (`tax_ids.txt`).

## Key Configuration (config.yaml)

- `mode`: `clade` (extract/filter all reads in the subtree) or `single` (exact matches only).
- `include`: `True` to **extract** matching reads; `False` to **filter out** matching reads.
- `sample_pattern` & `fastq_file_pattern`: Glob patterns to identify samples and their corresponding FASTQ files.

## Building and Running

### Prerequisites
- Conda or Mamba.
- Snakemake v8+ with the SLURM executor plugin (if running on a cluster).

### Environment Setup
Create the base Snakemake environment:
```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake snakemake-executor-plugin-slurm
conda activate snakemake
```

### Execution

#### Local Run
```bash
snakemake --cores <N> --use-conda --restart-times 4 --resources mem_mb=<MB>
```

#### Cluster Run (SLURM)
Ensure `slurm/config.yaml` is updated with your account and partition details.
```bash
snakemake --executor slurm --profile slurm --cores 1 --use-conda --restart-times 4
```

## Development Conventions

- **Snakemake Rules:** Defined in `rules/*.smk` and included in the main `Snakefile`.
- **Environment Management:** Each rule specifies its environment in `envs/`.
- **Resource Management:** Rules use dynamic resource functions (e.g., `get_mem_mb`, `get_runtime`) to automatically scale resources on retry attempts.
- **Performance:** `SeqKit` is the primary engine for FASTQ manipulation due to its speed and memory efficiency.
- **Logging:** All rules pipe stdout/stderr to files in the `logs/` directory.
- **Reproducibility:** Use the `--use-conda` flag to ensure consistent tool versions across different environments.
