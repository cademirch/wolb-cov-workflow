"""
Utility functions for workflow.
"""
import pandas as pd
from pathlib import Path
from snakemake.exceptions import WorkflowError


samples = pd.read_json(config["sample_sheet"], orient="records")  # noqa: F821
sampleids = samples["SampleID"].unique().tolist()


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


def get_single_references(wc):
    genomes = wc.genome.split("-")
    return [config["genome_paths"][g] for g in genomes]  # noqa: F821


def get_mappability_bed(wc):
    genome = get_combined_genome_name(wc)
    return expand(  # noqa: F821
        "results/mappability/{genome}/{e}/filtered.bed", genome=genome, e=wc.e
    )


def get_reads(wc):
    read_1 = samples.loc[samples["SampleID"] == wc.sample]["read1"].tolist()
    read_2 = samples.loc[samples["SampleID"] == wc.sample]["read2"].tolist()
    if len(read_1) == 1 and len(read_2) == 1:
        return {"read_1": read_1, "read_2": read_2}
    else:
        return {"read_1": "results/merged_fastqs/{sample}_R1.fq.gz", "read_2": "results/merged_fastqs/{sample}_R1.fq.gz"}

def get_reads_to_merge(wc):
    read_1 = samples.loc[samples["SampleID"] == wc.sample]["read1"].tolist()
    read_2 = samples.loc[samples["SampleID"] == wc.sample]["read2"].tolist()
    return {"read_1": read_1, "read_2": read_2}

def get_octomom(wc):
    s = samples.loc[samples["infection"] == "wMel"]

    return expand(  # noqa: F821
        "results/octomom/{sample}.regions.bed.gz", sample=s["seqID"].tolist()
    )
