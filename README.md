Table of Contents
=================
* [Wolbachia Coverage Workflow](#wolbachia-coverage-workflow)
   * [Usage](#usage)
      * [Environment Setup](#environment-setup)
      * [Workflow Configuration](#workflow-configuration)
         * [Reference genome](#reference-genome)
         * [Config.yaml](#configyaml)
            * [Coverage groups](#coverage-groups)
         * [Resources.yaml](#resourcesyaml)
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
| `sample_sheet` | Path to CSV sample sheet.| `str` |
| `resource_config` | Path to resource.yaml.| `str` |
| `genome_path` | Path to merged genome. MUST be uncompressed.| `str` |
| `coverage_groups` | Each subkey represents a coverage group, value must be list of contigs belonging to that group. See below for more info.| `dict[str:list[str]]` |
| `kmer_size` | Kmer size for genmap (-k).| `int` |
| `mismatches` | Allowed mismatches for genmap (-e).| `int` |
| `score_cutoff` | Minimum mappability score for coverage calculation.| `int` |
| `min_mapq` | Minimum mapq for filtering alignments.| `int` |

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

#### Resources.yaml
The `resource_config.yaml` specifes the computational resources to be used by each rule. It is set with reasonable defaults for the phoenix.prism cluster, but should be checked and adjusted if running elsewhere. For more info about specifying resources see [Snakemake's documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources).

#### Sample sheet
In order to define your samples and the path to their reads, you must construct CSV file with the following fields:
| Field | Description | Type |
| ---- | -------------| ------ |
| `sample_id` | Name of the sample (must be unique).| `str` |
| `read_1` | Path to read_1.| `str` |
| `read_2` | Path to read_2.| `str` |

### Running the workflow
After completing the configuration steps, the workflow is ready to be run. 

From the root directory of the repo and in your snakemake env:

#### Running locally
```
snakemake --cores <num cores to use>
```

#### Running on SLURM
```
srun snakemake --slurm --jobs <num jobs to run in parallel> 
```

### Workflow output
The workflow produces a CSV file reporting the average coverage in each `coverage_group` specified in the config file. The output file has the following fields:
| Field | Description | Type |
| ---- | -------------| ------ |
| `sample_id` | Name of the sample| `str` |
| `Name of coverage_group 1` | Average coverage of coverage_group 1| `float` |
| `Name of coverage_group N` | Average coverage of coverage_group N| `float` |




