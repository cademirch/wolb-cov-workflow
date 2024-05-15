Table of Contents
=================
* [Wolbachia Coverage Workflow](#wolbachia-coverage-workflow)
   * [Usage](#usage)
      * [Environment Setup](#environment-setup)
      * [Workflow Configuration](#workflow-configuration)
         * [Reference genome](#reference-genome)
         * [Config.yaml](#configyaml)
            * [Coverage groups](#coverage-groups)
         * [Sample sheet](#sample-sheet)
      * [Running the workflow](#running-the-workflow)
         * [Running locally](#running-locally)
         * [Running on SLURM](#running-on-slurm)
      * [Workflow output](#workflow-output)
# Wolbachia Coverage Workflow
This Snakemake workflow generates average coverage information from Illumina sequencing data.
## Usage
### Environment Setup
I suggest forking this repo and cloning your version. That way any changes you might make can be pulled into this repo.
```bash
$ git clone https://github.com/<youraccount>/wolb-cov-workflow.git
```
Then, setup your mamba env. I recommend following [Snakemake's instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba). 

Briefly:
```bash
$ conda activate base
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
$ mamba activate snakemake
```

### Workflow Configuration
In order to run the workflow, you must first setup the following files:
#### Reference genome
You must specify a path to a uncompressed reference genome you would like to map against. Typically this will be a 'merged' reference genome of host and symbiont(s) genomes. After downloading your desired genomes, you can merge them like so:
```{bash}
$ cat genome1.fa genome2.fa genome3.fa > merged.fa
```
#### Config.yaml
The file `config/config.yaml` controls the core workflow options. It contains the following fields:
| Option | Description | Type |
| ---- | -------------| ------ |
| sample_sheet | Path to CSV sample sheet.| `str` |
| resource_config | Path to resource.yaml.| `str` |
| genome_path | Path to merged genome. MUST be uncompressed.| `str` |
| coverage_groups | Each subkey represents a coverage group, value must be list of contigs belonging to that group. See below for more info.| `dict[str:list[str]]` |
| roles | Defines host and symbiont relationships. Values must match keys from `coverage_groups`.| `dict[str:list[str]]` |
| kmer_size | Kmer size for genmap (-k).| `int` |
| e_vals | Allowed mismatches for genmap (-e).| `list[int]` |
| score_cutoff | Minimum mappability score for coverage calculation.| `int` |
| min_mapq | Minimum mapq for filtering alignments.| `int` |

##### Coverage groups
We define coverage groups as contig(s) in the reference genome that we want to calculate the average coverage of. For example, we typically only care about the dmel6 autosomes in our analysis, so we would specify that like so:
```yaml
# config.yaml

coverage_groups:
    # Dmel autosomes
    dmel: ["NT_033779.5","NT_033778.4","NT_037436.4","NT_033777.3","NC_004353.4",]
    wmel: ["NC_002978.6"]
    wri: ["NC_012416.1"]
```
In this example, the average coverage coverage across the 5 dmel autosomes will be calculated and reported as `dmel` in the coverage CSV output file.

**It is important to check that the contigs you specify exist in the reference genome.**

##### Roles
We use roles to define host/symbiont(s) relationships for titer calculation. Your role definitions must have only one host and as many symbionts as you'd like. The names of each must match the names of your coverage groups. Continuing our example:
```yaml
# config.yaml

## Roles to define host and symbiont relationships. Only one host is allowed.
roles:
  hosts: ["dmel", "dsim"]
  symbionts: ["wmel","wri"]
```

#### Sample sheet
In order to define your samples and the path to their reads, you must construct JSON file with the following fields:
| Field                | Description                                              | Type                |
|----------------------|----------------------------------------------------------|---------------------|
| SampleID             | Unique identifier for the sample                         | `str`                 |
| Lane                 | Sequencing lane number                                   | `int`                 |
| sequencing_order_id  | Identifier for the sequencing order                      | `str`                 |
| read1                | Path to the first sequencing read file                   | `str`                 |
| read2                | Path to the second sequencing read file                  | `str`                 |
| host                 | Host genome name                                         | `str`                 |
| infection            | List of infections found in the sample                   | `List[str]`           |

There is a python script, `prepare_sample_sheet.py` that makes this JSON given a Duke Sequencing demulitplexing report. 

### Running the workflow
After completing the configuration steps, the workflow is ready to be run. 

I recommend creating a directory to store your data and results outside of the workflow repo. This creates better separation and prevents git clutter. You can create your directory where ever you'd like and then copy the `config/` from the workflow repo there. 

Example directory structure:
```
.
├── wolb-cov-workflow/
│   ├── config/
│   │   ├── config.yaml
│   │   └── samples.csv
│   └── workflow/
│       ├── envs
│       ├── rules
│       ├── scripts
│       └── Snakefile
└── your-directory/
    ├── your-data
    └── config/
        ├── config.yaml
        └── samples.csv
```


#### Running locally
```
snakemake -s <path/to/wolb-cov-workflow/workflow/Snakefile -d <path/to/your/data/directory> --use-conda --cores <num cores to use>
```

#### Running on SLURM
```
srun snakemake -s <path/to/wolb-cov-workflow/workflow/Snakefile -d <path/to/your/data/directory> --use-conda --slurm --jobs <num jobs to run in parallel> 
```

### Workflow output
The workflow produces a CSV file reporting the mean depth in each `coverage_group` specified in the config file. As well as the titer of each symbiont denoted in `roles`. The output file has the following fields:
| Field | Description | Type |
| ---- | -------------| ------ |
| `sample_id` | Name of the sample| `str` |
| `Name of coverage_group 1` | Mean depth of coverage_group 1| `float` |
| `Name of coverage_group N` | Mean depth of coverage_group N| `float` |
| `Symbiont_1_titer` | Titer of Symbiont 1| `float` |
| `Symbiont_N_titer` | Titer of Symbiont N| `float` |




