# -*- coding: utf-8 -*-

"""
Extract reads classified to certain tax_id(s) or clade(s) by Kraken 2.
"""

def get_mem_mb(wildcards, attempt):
    # Start at 4GB, double each time: 4, 8, 16, 32
    return 4000 * (2 ** (attempt - 1))

# Function to triple runtime: 30, 90, 270, 810 minutes
def get_runtime(wildcards, attempt):
    return 30 * (3 ** (attempt - 1))

# Parameters from the configuration that we'll use
tax_ids_file = config['tax_ids_file']               # The file with all tax_ids that we should extract reads from
mode = config["mode"]                               # Should we get reads hitting only the specified tax_id, or its clade?
names_file = config["taxonomy_names"]               # Path to the names file
nodes_file = config["taxonomy_nodes"]               # Path to the nodes file
k2_classification_file = config["sample_pattern"]   # The pattern of the Kraken 2 classification files. Just want to rename it for clarity here in the workflow
fastq_file_pattern = config["fastq_file_pattern"]   # The pattern of the fastq files

# Reads tax IDs from a file and returns them as a list
def read_tax_ids(file):
    with open(file) as f:
        return [line.strip() for line in f]

# Read tax IDs from the file
tax_ids = read_tax_ids(tax_ids_file)

# Fail-safe: Exit if no tax IDs are provided
if not tax_ids:
    raise ValueError(f"No tax IDs found in {tax_ids_file}. Check the file contents.")

# Determine the status string based on the include flag
# include: True -> "extracted"
# include: False -> "filtered"
status = "extracted" if config.get("include", True) else "filtered"

# The fastq output files containing taxonomic ID-specific reads
extracted_reads = expand(
    str(OUTDIR/"mode_{mode}/{tax_id}/{sample}_taxID-{tax_id}_{status}_{direction}.fastq.gz"),
    mode=mode,
    sample=SAMPLES,
    tax_id=tax_ids,
    status=status,
    direction=['R1', 'R2'])
all_outputs.extend(extracted_reads)

rule getAll_taxID_readIDs:
    conda:
        "../envs/extract_reads.yaml"
    input:
        kraken2 = lambda wildcards: INPUTDIR/k2_classification_file.format(sample=wildcards.sample)
    output:
        temp(OUTDIR/"mode_{mode}/{sample}_allTaxIDs_readIDs.txt")
    params:
        names = names_file,
        nodes = nodes_file,
        extract_mode = mode,
        tax_id_file = tax_ids_file
    log:
        LOGDIR / "getAll_taxID_readIDs_{sample}_{mode}.log"
    shell:
        """
        /usr/bin/time -v scripts/selmeout.py \
            --mode {params.extract_mode} \
            --nodes {params.nodes} \
            --names {params.names} \
            --output {output} \
            --input {input.kraken2} \
            --tax_id_file {params.tax_id_file} \
            2>&1 | tee {log}
        """

rule extract_taxon_readIDs:
    conda:
        "../envs/extract_reads.yaml"
    input:
        allReadIDs = OUTDIR/"mode_{mode}/{sample}_allTaxIDs_readIDs.txt"
    output:
        taxon_readIDs = temp(OUTDIR/"mode_{mode}/{tax_id}/{sample}_taxID-{tax_id}_readIDs.txt")
    log:
        LOGDIR / "extract_taxon_readIDs_{sample}_taxID-{tax_id}_{mode}.log"
    threads: 1
    resources:
        mem_mb = 1000,
        runtime = "15m"
    shell:
        """
        awk -v taxid={wildcards.tax_id} -F '\t' '$2 == taxid {{print $3}}' {input.allReadIDs} > {output.taxon_readIDs} 2> {log}
        """

rule extract_taxID_reads:
    conda:
        "../envs/extract_reads.yaml"
    input:
        fastq = lambda wildcards: INPUTDIR/fastq_file_pattern.format(sample=wildcards.sample, direction=wildcards.direction),
        taxon_readIDs = OUTDIR/"mode_{mode}/{tax_id}/{sample}_taxID-{tax_id}_readIDs.txt"
    output:
        fastq = str(OUTDIR/"mode_{mode}/{tax_id}/{sample}_taxID-{tax_id}_") + status + "_{direction}.fastq.gz"
    params:
        # include: True -> "" (keep matches)
        # include: False -> "-v" (invert matches)
        grep_flags = "" if config.get("include", True) else "-v",
    log:
        LOGDIR / "extract_taxID_reads_{sample}_taxID-{tax_id}_{mode}_{direction}.log"
    threads: 4
    resources:
        mem_mb = get_mem_mb,
        runtime = get_runtime
    shell:
        """
        /usr/bin/time -v seqkit grep \
            {params.grep_flags} \
            -f {input.taxon_readIDs} \
            {input.fastq} \
            -o {output.fastq} \
            -j {threads} \
            2>&1 | tee {log}
        """
