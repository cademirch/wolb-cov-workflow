"""
Utility functions for workflow.
"""
import pandas as pd
from pathlib import Path
from snakemake.exceptions import WorkflowError


samples = pd.read_json(config["sample_sheet"], orient="records")  # noqa: F821
sampleids = samples["SampleID"].unique().tolist()
mixed = samples[samples['infection'].apply(lambda x: isinstance(x, list) and len(x) > 1)]
mixed_ids = mixed["SampleID"].unique().tolist()

def get_combined_genome_name(wc):
    combined = []
    sample = samples[samples["SampleID"] == wc.sample].iloc[0]
    infection = sample["infection"]
    host = sample["host"]
    combined.append(host)
    if infection:
        combined.extend(infection)
    combined = sorted(combined)
    combined = "-".join(combined)
    return combined


def get_combined_genome(wc):
    combined = get_combined_genome_name(wc)
    d = {
        "ref": f"data/combined_genomes/{combined}.fa",
        "idx": f"data/combined_genomes/{combined}.fa.fai",
    }
    return d

def get_pileup_bed(wc):
    combined = get_combined_genome_name(wc)
    return f"data/bed/{combined}.bed"

def get_single_references(wc):
    genomes = wc.genome.split("-")
    return [config["genome_paths"][g] for g in genomes]  # noqa: F821


def get_mappability_bed(wc):
    genome = get_combined_genome_name(wc)
    return expand(  # noqa: F821
        "results/mappability/{genome}/{e}/filtered.bed", genome=genome, e=wc.e
    )


def get_reads(wc):
    """
    Get read files for a sample. If multiple read pairs exist,
    return the path to the merged files, otherwise return the direct paths.
    """
    read1_files = samples.loc[samples["SampleID"] == wc.sample]["read1_files"].iloc[0]
    read2_files = samples.loc[samples["SampleID"] == wc.sample]["read2_files"].iloc[0]
    
    if len(read1_files) == 1 and len(read2_files) == 1:
        # If only one pair of reads, return them directly
        return {"read_1": read1_files[0], "read_2": read2_files[0]}
    else:
        # If multiple pairs, return the path to the merged files
        return {
            "read_1": f"results/merged_fastqs/{wc.sample}_R1.fq.gz", 
            "read_2": f"results/merged_fastqs/{wc.sample}_R2.fq.gz"
        }

def get_reads_to_merge(wc):
    """
    Get all read files for a sample to be merged.
    Each read1 file corresponds to the read2 file at the same index.
    """
    read1_files = samples.loc[samples["SampleID"] == wc.sample]["read1_files"].iloc[0]
    read2_files = samples.loc[samples["SampleID"] == wc.sample]["read2_files"].iloc[0]
    
    return {"read_1": read1_files, "read_2": read2_files}

def get_octomom(wc):
    s = samples.loc[samples["infection"] == "wMel"]

    return expand(  # noqa: F821
        "results/octomom/{sample}.regions.bed.gz", sample=s["seqID"].tolist()
    )
