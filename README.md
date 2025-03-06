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

**For UCSC Phoenix users, these values have mostly been filled in for you**

#### Reference genome
You must specify a path to a uncompressed reference genome(s) you would like to map against.

#### Config.yaml
The file `config/config.yaml` controls the core workflow options. It contains the following fields:

| Option | Description | Type |
| ---- | -------------| ------ |
| sample_sheet | Path to JSON sample sheet created by prepare_samplesheet.py | `str` |
| recombination | Whether to run recombination analysis steps | `bool` |
| genome_paths | Dictionary mapping genome IDs to their respective FASTA file paths | `dict[str:str]` |
| genome_contigs_for_depth_calculation | Dictionary mapping genome IDs to lists of contigs to include in depth calculations | `dict[str:list[str]]` |
| roles | Defines host and symbiont relationships for titer calculations | `dict[str:list[str]]` |
| mappability | Options for genmap read mappability analysis | `dict` |
| mappability.e_vals | Allowed mismatches for genmap | `list[int]` |
| mappability.kmer_size | Kmer size for genmap (-k) | `int` |
| mappability.mismatches | Number of mismatches allowed | `int` |
| mappability.score_cutoff | Minimum mappability score for coverage calculation | `int` |
| filter_bam | Options for BAM file filtering | `dict` |
| filter_bam.min_mapq | Minimum mapping quality for filtering alignments | `int` |

##### genome_contigs_for_depth_calculation
This dictionary defines which contigs from each reference genome should be included in depth calculations. This is useful to focus analysis on main chromosomes and exclude smaller contigs like mitochondria or unplaced scaffolds:

```yaml
genome_contigs_for_depth_calculation:
  dmel: ["NT_033779.5","NT_033778.4","NT_037436.4","NT_033777.3","NC_004353.4"]
  wmel: ["NC_002978.6"]
  # Add more genomes as needed
```

##### roles
This dictionary defines the relationships between hosts and symbionts for calculating titers:

```yaml
roles:
  hosts: ["dmel", "dsim"]  # List of host genome IDs
  symbionts: ["wmel", "wri", "wwil"]  # List of symbiont genome IDs
```

The genome IDs must match those used in `genome_paths` and `genome_contigs_for_depth_calculation`.

#### Sample Sheet
In order to define your samples and the path to their reads, you must construct a JSON file with the following fields:

| Field        | Description                                      | Type         |
|--------------|--------------------------------------------------|--------------|
| SampleID     | Unique identifier for the sample                 | `str`        |
| host         | Host genome name (e.g., "dmel", "dsim")          | `str`        |
| infection    | List of infections found in the sample           | `List[str]`  |
| read1_files  | List of paths to first sequencing read files     | `List[str]`  |
| read2_files  | List of paths to second sequencing read files    | `List[str]`  |

### Creating the Sample Sheet

There is a Python script, `updated_sample_sheet.py`, that generates this JSON given a sample information CSV file. The script can search multiple read directories to find matching FASTQ files for each sample.

#### Usage:

```bash
python updated_sample_sheet.py -i samples.csv -d /path/to/reads_dir1 -d /path/to/reads_dir2 -o sample_sheet.json
```

#### Arguments:

- `-i, --input`: Path to the input CSV file with sample information
- `-d, --directory`: Path to directory containing read files (can be specified multiple times)
- `-o, --output`: Output JSON file path (defaults to stdout)
- `-v, --verbose`: Enable verbose logging

#### Input CSV Format:

The input CSV should contain at least the following columns:
- `Sample ID`: Unique sample identifier
- `Cell Line`: Host cell line (will be mapped to appropriate host genome)
- `Infection`: Infection status or type
- `Date Collected`: Date when sample was collected (rows with empty dates will be skipped)

The script will:
1. Automatically find all read files matching each sample ID in the specified directories
2. Map the Cell Line values to appropriate host genomes
3. Process infection information
4. Generate a properly formatted JSON file for use with the workflow

#### Example:

Sample CSV content:
```
Experiment ID,Sample ID,Sample #,Initials,Date Collected,Cell Line,Infection
Exp001,JJ_MW_1,1,JJ,2023-01-15,S2,Wmel
Exp001,JJ_MW_2,2,JJ,2023-01-15,JW18,Wmel:Wwil
```

Generated JSON (abbreviated):
```json
[
  {
    "SampleID": "JJ_MW_1",
    "host": "dmel",
    "infection": ["wmel"],
    "read1_files": ["/path/to/JJ_MW_1_S1_L001_R1_001.fastq.gz"],
    "read2_files": ["/path/to/JJ_MW_1_S1_L001_R2_001.fastq.gz"]
  },
  {
    "SampleID": "JJ_MW_2",
    "host": "dmel",
    "infection": ["wmel", "wwil"],
    "read1_files": [
      "/path/to/JJ_MW_2_S2_L001_R1_001.fastq.gz",
      "/path/to/JJ_MW_2_S2_L002_R1_001.fastq.gz"
    ],
    "read2_files": [
      "/path/to/JJ_MW_2_S2_L001_R2_001.fastq.gz",
      "/path/to/JJ_MW_2_S2_L002_R2_001.fastq.gz"
    ]
  }
]
```

### Running the workflow
After completing the configuration steps, the workflow is ready to be run. 

I recommend creating a directory to store your data and results outside of the workflow repo. This creates better separation and prevents git clutter. You can create your directory where ever you'd like and then copy the `config/` from the workflow repo there. 

Example directory structure:
```
.
├── wolb-cov-workflow/
│   ├── config/
│   └── workflow/
│       ├── envs
│       ├── rules
│       ├── scripts
│       └── Snakefile
└── your-directory/
    └── config/
        ├── config.yaml
        └── samples.json
```


#### Running locally
```
conda activate snakemake
snakemake -s <path/to/wolb-cov-workflow/workflow/Snakefile -d <path/to/your/data/directory> --use-conda --cores <num cores to use>
```

#### Running on SLURM
Make sure to install the slurm executor (only need to be done once):
```
conda activate snakemake
pip install snakemake-executor-plugin-slurm
```
Then run,
```
conda activate snakemake
snakemake -s <path/to/wolb-cov-workflow>/workflow/Snakefile -d <path/to/your/data/directory> --profile <path/to/wolb-cov-workflow>/workflow-profiles/slurm
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




