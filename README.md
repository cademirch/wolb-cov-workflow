# Wolbachia Coverage Workflow
This Snakemake workflow generates average coverage information from Illumina sequencing data.
## Usage
### Environment Setup
First, clone this repo to wherever you are working:
```bash
$ git clone https://github.com/cademirch/wolb-cov-workflow.git
```
Then, setup your mamba env. I recommend following [Snakemake's instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba). 

Briefly:
```bash
$ conda activate base
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
$ mamba activate snakemake
```

### Workflow Configuration
The workflow requires a few files in order to be run.
    