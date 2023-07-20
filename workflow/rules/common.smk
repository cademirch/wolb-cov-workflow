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
samples = pd.read_table(config["sample_sheet"], sep=",", dtype=str).replace(
    " ", "_", regex=True
)

with open(config["resource_config"], "r") as f:
    resource_config = safe_load(f)


def get_reads(wc):
    read_1 = samples.loc[samples["sample_id"] == wc.sample]["read_1"].item()
    read_2 = samples.loc[samples["sample_id"] == wc.sample]["read_2"].item()
    return {"read_1": read_1, "read_2": read_2}
