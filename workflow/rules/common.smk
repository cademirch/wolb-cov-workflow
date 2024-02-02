"""
Utility functions for workflow.
"""
import pandas as pd
from pathlib import Path
from yaml import safe_load
from snakemake.exceptions import WorkflowError


def validate_config():
    for k, v in config["coverage_groups"].items():
        if len(v) < 1:
            raise WorkflowError(
                f"Coverage group: `{k}` has no contigs listed, please correct this. "
            )
        role_values = []
        for j, r in config["roles"].items():
            if isinstance(r, str):
                role_values.append(r)
            elif isinstance(r, list):
                role_values.extend(r)
            else:
                raise TypeError(
                    f"Wrong type for role '{j}', must be list or str. Found '{type(r)}'"
                )
        if k not in role_values:
            raise WorkflowError(f"Coverage group: `{k}` has no role.")


validate_config()

# print(config["coverage_groups"].items())
## Read in sample sheets and resources yaml ##
samples = pd.read_table(config["sample_sheet"], sep=",", dtype=str)


with open(config["resource_config"], "r") as f:
    resource_config = safe_load(f)


def get_reads(wc):
    read_1 = samples.loc[samples["seqID"] == wc.sample]["read1"].tolist()
    read_2 = samples.loc[samples["seqID"] == wc.sample]["read2"].tolist()
    if len(read_1) == 1 and len(read_2) == 1:
        return {"read_1": read_1, "read_2": read_2}
    else:
        return {"read_1": "results/merged_fastqs/{sample}_R1.fq.gz", "read_2": "results/merged_fastqs/{sample}_R1.fq.gz"}

def get_reads_to_merge(wc):
    read_1 = samples.loc[samples["seqID"] == wc.sample]["read1"].tolist()
    read_2 = samples.loc[samples["seqID"] == wc.sample]["read2"].tolist()
    return {"read_1": read_1, "read_2": read_2}

def get_octomom(wc):
    s = samples.loc[samples["infection"] == "wMel"]
    
    return expand("results/octomom/{sample}.regions.bed.gz", sample=s["seqID"].tolist())