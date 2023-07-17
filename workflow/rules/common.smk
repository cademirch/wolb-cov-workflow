"""
Utility functions for workflow.
"""
import pandas as pd
from pathlib import Path
from yaml import safe_load
from snakemake.exceptions import WorkflowError

for k, v in config["coverage_groups"].items():
    if len(v) < 1:
        raise (
            WorkflowError(
                f"Coverage group: `{k}` has no contigs listed, please correct this. "
            )
        )

## Read in sample sheets and resources yaml ##
samples = pd.read_table(config["sample_sheet"], sep=",", dtype=str).replace(
    " ", "_", regex=True
)

with open(config["resource_config"], "r") as f:
    resource_config = safe_load(f)


def get_reads(wc):
    read_1 = samples.loc[samples["sample_id"] == wc.sample]["read_1"].item()
    read_2 = samples.loc[samples["sample_id"] == wc.sample]["read_2"].item()
    return {"read_1": read_1, "read_2": read_2}