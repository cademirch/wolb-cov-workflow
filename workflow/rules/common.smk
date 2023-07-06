"""
Utility functions for workflow.
"""
import pandas as pd
from yaml import safe_load

## Read in sample sheets and resources yaml ##
samples = pd.read_table(config["sample_sheet"], sep=",", dtype=str).replace(' ', '_', regex=True)
with open(config["resource_config"], "r") as f:
    resources = safe_load(f)