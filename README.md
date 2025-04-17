[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.1-23aa62.svg)](https://www.nextflow.io/) [![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/) ![Nf-test](https://img.shields.io/badge/NFtest-%E2%89%A50.8.4-23aa62.svg?labelColor=0000)

# Viral Lens

The **Viral Lens** is a bioinformatic pipeline deal with short-read sequencing data generated from the bait-capture protocols for enrichment designed under the context of the [RVI project](https://www.sanger.ac.uk/group/respiratory-virus-and-microbiome-initiative/) by [GSU](https://www.sanger.ac.uk/collaboration/genomic-surveillance-unit/).
This pipeline will generate, if possible, high quality consensus sequences for every virus identified according to a kraken database.

---
## Contents
- [Contents](#contents)
- [Pipeline Summary](#pipeline-summary)
- [How to Cite](#how-to-cite)
- [Quick Start](#quick-start)
- [Installation](#installation)
  - [dependencies](#dependencies)
  - [build containers](#build-containers)
- [Usage](#usage)
- [Inputs](#inputs)
  - [Manifest](#manifest)
  - [Kraken Database](#kraken-database)
- [Outputs](#outputs)
  - [Main output files](#main-output-files)
  - [Intemediate output files](#intemediate-output-files)
- [Configuration Sections](#configuration-sections)
  - [Parameters](#parameters)
    - [input](#input)
    - [Pipeline Flow Controls](#pipeline-flow-controls)
    - [Kraken Database Parameters](#kraken-database-parameters)
    - [Kraken2Ref Handling](#kraken2ref-handling)
    - [Kraken2ref Report Filter](#kraken2ref-report-filter)
    - [iVar Parameters](#ivar-parameters)
    - [Virus Subtyping](#virus-subtyping)
    - [Consensus Generation](#consensus-generation)
    - [Output Directory](#output-directory)
    - [General](#general)
  - [Containerization](#containerization)
  - [Process Settings](#process-settings)
  - [Profiles](#profiles)
    - [`sanger_standard` Profile](#sanger_standard-profile)
- [Unit Tests](#unit-tests)
- [Pipeline components documentation](#pipeline-components-documentation)
  - [Processes](#processes)
      - [run\_kraken](#run_kraken)
      - [run\_k2r\_sort\_reads](#run_k2r_sort_reads)
      - [run\_k2r\_dump\_fastqs\_and\_pre\_report](#run_k2r_dump_fastqs_and_pre_report)
      - [concatenate\_fqs\_parts](#concatenate_fqs_parts)
      - [get\_taxid\_references](#get_taxid_references)
      - [bwa\_alignment\_and\_post\_processing](#bwa_alignment_and_post_processing)
      - [run\_ivar](#run_ivar)
      - [run\_pagolin](#run_pagolin)
      - [run\_qc\_script](#run_qc_script)
      - [write\_classification\_report](#write_classification_report)
  - [Workflow](#workflow)
    - [SORT\_READS\_BY\_REF](#sort_reads_by_ref)
    - [GENERATE\_CONSENSUS](#generate_consensus)
    - [COMPUTE\_QC\_METRICS](#compute_qc_metrics)
    - [SCOV2\_SUBTYPING](#scov2_subtyping)
    - [GENERATE\_CLASSIFICATION\_REPORT](#generate_classification_report)
  - [Custom pipeline scripts](#custom-pipeline-scripts)
    - [Kraken2ref JSON to TSV Report Script](#kraken2ref-json-to-tsv-report-script)
    - [QC Script for BAM and FASTA Files](#qc-script-for-bam-and-fasta-files)
- [Licence](#licence)

---

## Pipeline Summary

The pipeline takes a manifest containing  **fastq pairs file** paths and a **kraken detabase** as inputs (check [Inputs section](#inputs) for more details) and outputs a **classification report**, **consensus sequences** and a **collection of intermediate files**. Here is an broad overview of the pipeline logic

0. **Preprocessing** (Optional): An optional preprocessing workflow is provided on this pipeline. The Preprocessing pipeline remove adapters (via `trimmomatic`), tandem repeats (via `TRF`) and remove human reads (via `sra-human-scrubber`) from fastq files.

1. **Sort Reads**: The initial step is sort reads using `kraken2` for each fastq pairs according to the database provided. The classified reads is used as input to [kraken2ref](https://github.com/genomic-surveillance/kraken2ref) which will generate one pair of fastq files per taxid found.
   - An option to split big files is provided (check [Parameter section](#parameters)).

2. **Generate Consensus**: After all samples been classified, all references observed for that samples batch are fetch from the `kraken database` (or an arbitrary fasta file provided by the user). The classfied reads are aligned to their respective references (via `bwa`). The alignment is used as input for `ivar` to obtain a consensus sequence.

3. **Compute QC**: QC metrics are computed via `samtools` and a custom script (`qc.py`)


4. **SARS-CoV-2 Subtyping**: SARS-CoV-2 subtyping can be done if present on the sample

>**NOTE**: An optional preprocessing workflow is provided 

---

[**(&uarr;)**](#contents)


## How to Cite

This  software will be published soon. Until it is, please provide the URL to this GitHub repository when you use the software in your own work.

[**(&uarr;)**](#contents)

---
## Quick Start

Assuming [dependencies](#dependencies) are installed on the system:

1. Setup the pipeline:

```bash
# clone the repo
git clone --recursive https://github.com/genomic-surveillance/rvi-viral-lens
cd viral_pipeline/
# (optional) build containers
cd containers/
sudo singularity build base_container.sif baseContainer.sing
sudo singularity build ivar.sif ivarContainer.sing
sudo singularity build pangolin.sif pangolinContainer.sing
sudo singularity build kraken.sif krakenContainer.sing
sudo singularity build kraken2ref.sif kraken2ref.sing

# (optional) Download database
cd ../
wget https://rvi_kraken2_dbs.cog.sanger.ac.uk/refseq_ncbiFlu_kfv2_20241027.tar.gz
tar -xf  refseq_ncbiFlu_kfv2_20241027.tar.gz
```

2. Run pipeline

You will need a manifest and a kraken database (check [Inputs section](#inputs) and [Usage section](#usage) for more details)

```bash
PIPELINE_CODES=<path to viral lens repo>
MANIFEST=<path to my manifest>
kraken_db_path=<path to my kraken DB>
PIPELINE_CONTAINERS=<path to my containers dir>

nextflow run ${PIPELINE_CODES}/main.nf --manifest ${MANIFEST} \
    --db_path ${kraken_db_path} \
    --outdir ./output/ \
    --containers_dir ${PIPELINE_CONTAINERS} \
    -with-trace -with-report -with-timeline \
    -profile sanger_local \
    -resume
```

[**(&uarr;)**](#contents)

---
## Installation

### dependencies

- [Nextflow](https://www.nextflow.io) (tested on `23.10.1`)
- [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) (required to use Singularity Containers, tested on ``ce version 3.11.4``)

> We strongly recommend to run the pipeline using the containers provided at [quay.io gsu-pipeline](https://quay.io/organization/gsu-pipelines).

If not using containers, all the software needs to be available at run time. Here is a list of all the softwares versions currently in use in the pipeline as set on each container.

- **Base container**

  - Samtools = `1.21`
  - BWA = `0.7.17`
  - biopython = `1.79`
  - pysam = `0.22.0`
  - pandas = `1.1.5`
  - matplotlib = `3.3.4`

- **Ivar container**
  - samtools = `1.21`
  - BWA = `0.7.17`
  - iVar = `1.4.3`

- **Kraken2Ref cotainer**:
  - pytest = `6.2.2`
  - importlib-resources = `5.1.0`
  - flake8 = `7.0.0`
  - pandas = `2.1.4`
  - cached-property = `1.5.2`
  - scipy = `1.12.0`
  - kraken2ref = `v2.1.0`

- **Kraken containers**
  - kraken2 = `v2.1.3`
  - kraken_tools = `v1.2`
  - biopython = `v1.81`

> NOTE: If not running under the provided containers, is strongly recommended to use the same versions described here.

### build containers

By default, the pipeline will use the containers provided at [quay.io gsu-pipeline](https://quay.io/organization/gsu-pipelines).

Singularity and Docker recipes for the containers used on this pipeline are available on this repository at `containers/` dir. To build the containers, run the commands bellow.

```{bash}
cd containers/
sudo singularity build base_container.sif baseContainer.sing
sudo singularity build ivar.sif ivarContainer.sing
sudo singularity build pangolin.sif pangolinContainer.sing
sudo singularity build kraken.sif krakenContainer.sing
sudo singularity build kraken2ref.sif kraken2ref.sing
```

> NOTE: To use local containers on the pipeline set the parameters `use_local_containers` to `true` and `use_registry_containers` to `false` either at the config file or at the CLI.

[**(&uarr;)**](#contents)

## Usage

1. **Generate manifest**
For convenience, a script to generate the manifest is provided on this repo:

```{bash}
python write_manifest.py ./path/to/my/fastqs_dir/ -fq1_ext my_r1_ext -fq2_ext my_r2_ext
```

2. **Run pipeline**

```{bash}
nextflow run /path/to/rvi_consensus_gen/main.nf --manifest /path/to/my/manifest.csv \
        --db_path /path/to/my/kraken_db \
        --outdir outputs/ \
        --containers_dir /path/to/my/containers_dir/ \
        -profile sanger_stantard -resume -with-trace -with-report -with-timeline
```

[**(&uarr;)**](#contents)

## Inputs

This pipeline relies on two **main inputs**:

- **`manifest`** : CSV Manifest of input fastq file pairs.
  - Must have `sample_id`,`reads_1` and `reads_2` columns
  - If you have your set of fastq pairs in a single dir, a script (`write_manifest.py`) is provided to facilitate this process.

- **`db_path`** : Path of a valid [kraken2 database](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#kraken-2-databases)

### Manifest


The pipeline require as input a manifest containing a unique sample id (`sample_id`) and paths to each of the fastq pair file (`reads_1` and `reads_3`)

```csv
sample_id,reads_1,reads_2
sample1,/path/to/output/sample1_R1.fq,/path/to/output/sample1_R2.fq
sample2,/path/to/output/sample2_R1.fq,/path/to/output/sample2_R2.fq
sample3,/path/to/output/sample3_R1.fq,/path/to/output/sample3_R2.fq
```

> NOTE: The collumn **`sample_id`** must be **unique, alphanumerics (non consecutive "_" are accepted) and cannot be empty**. Pipeline will fail if any of these conditions are not met.

**Write Manifest Script**

For user convenience, a script to write a manifest (`write_manifest.py`). This script generates a CSV manifest file from a directory of FASTQ files.

1. It takes a glob pattern to locate files, checks for paired-end FASTQ files (based on provided extensions)
2. extracts metadata such as sample IDs (optionaly taxonomic IDs) from the filenames.
    1. The user can specify if taxonomic IDs are included in the filenames, along with the filename separator. 
    2. If reference files are provided for taxonomic IDs, they are linked to the samples in the output. 

The script handles errors such as missing paired-end files or improperly formatted filenames, and generates a CSV manifest that includes sample information, file paths, and reference data (if applicable).

Example Command:

```bash
python write_manifest.py "output/*/reads_by_taxon/*.extracted_{1,2}.fq" \
    -fq1_ext "_R1.fq" \
    -fq2_ext "_R2.fq" \
    -filename_sep "." \
    -manifest_out "/path/to/output/manifest.csv"
```

**Breakdown**:

- `"output/*/reads_by_taxon/*.extracted_{1,2}.fq"` - A glob pattern to locate the FASTQ files in the specified directory structure.
- `-fq1_ext "_R1.fq"` - Specifies the naming pattern for the forward FASTQ files (e.g., `_R1.fq`).
- `-fq2_ext "_R2.fq"` - Specifies the naming pattern for the reverse FASTQ files (e.g., `_R2.fq`).
- `-filename_sep "."` - The separator used in the filename to extract sample id (e.g., `.`), anything before the separator will be considered as the value which should be the `sample_id`
- `-manifest_out "/path/to/output/manifest.csv"` - The output path for the generated CSV manifest file.

### Kraken Database

This pipeline was developed under the RVI project and a modified Kraken2 database was developed to better account the phyogenetic structure of Flu and RSV.

For the RVI project, this pipeline was run with a custom database built using [kraken_flu](https://github.com/genomic-surveillance/krakenflu). 
This database is based on the NCBI viral taxonomy with the following modifications to the taxonomy structure:
- The Influenza taxonomy is modified below the level of the species such that each segment is represented as its own distinct branch. Further, for specifically segments 4 (HA) and 6 (NA), those branches of the custom taxonomy have an additional level so that nodes representing subtypes H1, H2, H3... are directly below the segment 4 node, and similarly, nodes representing subtypes N1, N2, N3... fall directly under the segment 6 node. Essentially, the segments of the flu genome are considered as distinct sequences for the purposes of the kraken2 classification step
- The number of Influenza sequences in the database is also expanded as compared to those found in viral RefSeq
- The number of RSV A and B sequences has also been increased as compared to those found in viral RefSeq

To download the database used on the RVI project, run the `wget` command bellow or check [here](https://rvi_kraken2_dbs.cog.sanger.ac.uk/refseq_ncbiFlu_kfv2_20241027.tar.gz
):

```
wget https://rvi_kraken2_dbs.cog.sanger.ac.uk/refseq_ncbiFlu_kfv2_20241027.tar.gz
```

[**(&uarr;)**](#contents)

## Outputs

The output file tree should look like the tree bellow:

```bash
<output_dir>/
├── <sample_id>
│   ├── <taxid>
│   │   ├── <sample_id>.<taxid>.consensus.fa
│   │   ├── <sample_id>.<taxid>.qc.csv
│   │   ├── <sample_id>.<taxid>.sorted.bam
│   │   └── <sample_id>.<taxid>.sorted.bam.bai
│   ├── [...]
│   ├── <sample_id>.class_seqs_1.fq
│   ├── <sample_id>.class_seqs_2.fq
│   ├── <sample_id>.kraken.output
│   ├── <sample_id>.report.txt
│   ├── <sample_id>.unclass_seqs_1.fq
│   ├── <sample_id>.unclass_seqs_2.fq
├── [...]
├── classification_report.csv
```

### Main output files

- **Report CSV**: A csv file summarizing consensus sequences information and qc metrics obtained from each input samples.
  - location: `<output_dir>/classification_report.csv`
  - output of write_classification_report process

- **Consensus sequences per taxid**: the consensus sequence obtained for a given `taxid` and `sample_id` combination.
  - location: `<output_dir/<sample_id>/<sample_id>.<taxid>.consensus.fa`
  - output of [run_ivar process](#runivar)

### Intemediate output files

**Kraken Output files** generated by run_kraken process

- `<sample_id>.class_seqs_{1,2}.fq`: The pair of fastq files containing all reads which were associated to a `taxid` in the database.
- `<sample_id>.unclass_seqs_{1,2}.fq`: The pair of fastq files containing all reads which were not classified to any `taxid` in the database
- `<sample_id>.report.txt`: A `tsv` file sumarizing the number of reads associated to a given item in the taxonomic tree of the kraken database. For more details, check [this file format kraken2 documentation](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#sample-report-output-format)

  - Here is an example of what the content of this file should look like:

```tsv
81.03	3954901	3954901	U	0	unclassified
18.97	926041	0	R	1	root
18.97	926041	29	D	10239	  Viruses
18.76	915890	0	D1	2559587	    Riboviria
18.76	915890	0	K	2732396	      Orthornavirae
17.42	850471	0	P	2732408	        Pisuviricota
17.42	850445	0	C	2732506	          Pisoniviricetes
17.42	850445	0	O	76804	            Nidovirales
17.42	850445	0	O1	2499399	              Cornidovirineae
17.42	850445	0	F	11118	                Coronaviridae
[...]
```

- `<sample_id>.kraken.output`: a `tsv` file containing classification of each reads id. For more details, check [this kraken2 output file specific documentation](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#standard-kraken-output-format)
  - here is an example of its content:

```
U	A00948:701:HJHWGDSX7:2:1101:4218:1031	0	151|151	A:1 0:116 |:| 0:117
U	A00948:701:HJHWGDSX7:2:1101:9498:1031	0	151|151	A:1 0:116 |:| 0:117
U	A00948:701:HJHWGDSX7:2:1101:7699:1047	0	151|151	A:1 0:116 |:| 0:117
C	A00948:701:HJHWGDSX7:2:1101:15058:1047	2697049	151|151	A:1 2697049:116 |:| 2697049:117
U	A00948:701:HJHWGDSX7:2:1101:30481:1047	0	151|151	A:1 0:116 |:| 0:117
```

The **BWA output BAM files** generated by the `bwa_alignment_and_post_processing` process

- `<sample_id>.<taxid>.sorted.bam`: A [Binary Alignment Map (BAM)](https://en.wikipedia.org/wiki/Binary_Alignment_Map) file resulted from the alignment between the `<sample_id>.<taxid>.fastq` pair outputed by kraken2ref output (`run_kraken2ref_dump_fastq`).
- `<sample_id>.<taxid>.sorted.bai`: the index file of the `sorted.bam` file. the BAI file provide fast retrieval of alignments overlapping a specified region without going through the whole alignments, check [section 5 of the samtools documenation](https://samtools.github.io/hts-specs/SAMv1.pdf) for a more detailed description.

The **QC report** which is the output of `run_qc_script`.

- `<sample_id>.<taxid>.qc.csv`: a csv containing QC metrics of the consensus sequence obtained (`<sample_id.<taxid>.consensus.fa`).

```csv
sample_name,pct_N_bases,pct_covered_bases,longest_no_N_run,num_aligned_reads,fasta,bam,qc_pass,total_mapped_reads,ivar_md,total_unmapped_reads
<sample_id.<taxid>,0.76,99.19,3543,129724,<sample_id.<taxid>.consensus.fa,<sample_id.<taxid>.sorted.bam,TRUE,129837,10,0
```

[**(&uarr;)**](#contents)

## Configuration Sections

### Parameters

The `params` block defines key user-modifiable settings for the workflow.

#### input
- `manifest`: Path to the manifest file (default: `null`).

#### Pipeline Flow Controls

- `entry_point` [Default = `"sort_reads"`]: Entry point for the workflow. It can be set to `"sort_reads"` or `"consensus_gen"`.
- `do_scov2_subtyping` [Default = `true`] : Switch SARS-CoV-2 on and off.
- `scv2_keyword` [Default = `"Severe acute respiratory syndrome coronavirus 2"`] : keyword to identified samples with SARS-CoV-2. This string is obtained from the kraken2 database, therefore, it should be in line with the database in use.

#### Kraken Database Parameters

- `db_path`: Path to the Kraken database.
- `db_library_fa_path` (OPTIONAL): Path to the Kraken database library FASTA file.
  - By default, it assumes there is a `${params.db_path}/library/library.fna`.

#### Kraken2Ref Handling

- `k2r_fq_load_mode`: Loading mode for Kraken2 fastq files (either `full` or `chunks`).
  - Default: `"full"`.
- `k2r_max_total_reads_per_fq`: Maximum number of reads to process per fastq file.
  - Default: `10,000,000`.
- `k2r_dump_fq_mem`: Memory allocated for dumping fastq files.
  - Default: `"6 GB"`.

#### Kraken2ref Report Filter

- `min_reads_for_taxid`: Minimum number of reads required to assign a taxonomic ID.
  - Default: `100`.

#### iVar Parameters

- `ivar_min_depth`: Minimum depth for consensus calling.
  - Default: `10`.
- `ivar_freq_threshold`: Sets base frequency threshold for consensus calling.
  - Default: `0.75`.

#### Virus Subtyping

- `scv2_keyword`: Keyword to identify SARS-CoV-2 sequences. Any taxid name equal to the string set by this parameter will be considered as SCOV2 and subjected to specific SARS-CoV-2 subtyping.
  - Default: `"Severe acute respiratory syndrome coronavirus 2"`.

- `do_scov2_subtyping`: Boolean flag to enable or disable SARS-CoV-2 subtyping via Pangolin.
  - Default: `true`.

#### Consensus Generation

- `consensus_mnf` [OPTIONAL]: Path to the consensus manifest file (default: `null`).

#### Classification report

- `min_coverage_percent`: minimum coverage percent threshold, fractional values between 0.0 and 1.0 are accepted (default: `0.1`) 

#### Output Directory

- `outdir`: Directory where output files will be published.
  - Default: `"$launchDir/output/"`.

#### General

- `containers_dir` [DEFAULT =  `containers/` dir of this repository] : By default, the pipeline relies on Singularity containers and __assumes__ all containers are present on this directory and were named on a specific manner
- `outdir` [DEFAULT = `$launchDir/output/`] : set where output files should be published. By default, it will write files to an `output/` dir (if not existent, it will be created) at pipeline launch directory.

### Containerization

Currently, the pipeline only provide Singularity containers.

**Docker**

- `enabled`: Flag to enable or disable Docker.
  - Default: `false`.

**Singularity**

- `enabled`: Flag to enable or disable Singularity.
  - Default: `true`.

By default, the current version of the pipeline assumes all singularity containers are available on specific paths with specific names defined at `./conf/containers.config`.

- Currently, `"$projectDir/containers/"` is the default location for the containers, it can be changed by the user using `containers_dir` parameter.

All containers used by this pipeline recipes can be found at `./containers/` dir of this repository. The current containers in use are:

- `base_container.sif`: the default container for all processes unless overridden.
- `ivar.sif`: container to be used for processes using the `ivar` label
- `kraken.sif`: container to be used for processes using the `kraken` label
- `pangoling.sif`: container to be used for processes using the `pangolin` label
- `kraken2ref.sif`: container to be used for processes using the `kraken2ref` label.

### Process Settings

- `cache='lenient'`: Defines the cache behavior for processes, allowing cached results to be reused even when minor changes occur.
- `executor='local'`: The default executor is set to local, meaning processes will run on the local machine unless otherwise specified.

### Profiles

Profiles define environment-specific configurations. Currently, there is a single predefined profile: `sanger_standard`. User should consider to write their own, check [profiles Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles) for more details

#### `sanger_standard` Profile

This profile is optimized for running jobs on the Sanger Institute's infrastructure and includes settings for using Singularity, job execution on the LSF scheduler, and handling of specific processes.

**Singularity Settings**

- `autoMounts`: Automatically mounts paths required for job execution.
- `cacheDir`: Cache directory set to "`$PWD`" (the current working directory).
- `runOptions`: Singularity run options to bind necessary paths (`/lustre`, `/nfs`, `/software`, `/data/`).

**Process-Specific Settings**
Inherited from the global process configuration
- `cache='lenient'`: this makes `resume` more tolerant to changes in files attributes, such as timestamps. This is handy when using distributed file system which uses [Network File System protocols](https://en.wikipedia.org/wiki/Network_File_System), check [Nextflow documentation](https://www.nextflow.io/docs/latest/cache-and-resume.html#inconsistent-file-attributes) for more details.
- `executor='local'`: Default executor is local unless overridden.

**LSF Job Execution**
Certain processes are configured to run on the LSF scheduler with specific resource allocations:

- `run_k2r_sort_reads`:
  - Executor: `lsf`
  - CPU cores: `1`

- `run_k2r_dump_fastqs_and_pre_report`:
  - Executor: `lsf`
  - CPU cores: `1`

- `run_kraken`:
  - Executor: `lsf`
  - CPU cores: `16`
  - Memory: `1 GB`

**Job Naming and Memory Management**

- `jobName`: Custom job name format for LSF jobs, based on the task name and tag (`"RVI-preprocess - $task.name - $task.tag"`).

- `perJobMemLimit=true`: Ensures that memory limits are set per job.

[**(&uarr;)**](#contents)

---

## Unit Tests

The workflow & process unit tests for this pipeline are written in the [nf-test](https://www.nf-test.com/) (`v0.8.4`) Nextflow testing framework.

**Running Tests**

The `nf-test` looks for an environment variable (`CONTAINER_DIR`) to set the containers directory. Therefore, set this variable before running `nf-test`.

```{bash}
export CONTAINER_DIR=<my/container/dir/path>
```

The following command if entered from the repository top-level directory can be used to execute all of the per-process & per-workflow unit tests:

**Run all tests**

```{bash}
nf-test test ./
```

**Run individual module/workflow test**

```{bash}
nf-test test tests/<modules or workflows>/<module_to_test>.nf.test
```

[**(&uarr;)**](#contents)

---
## Pipeline components documentation

### Processes


##### run_kraken

Executes Kraken2 on paired-end FASTQ files, producing outputs that include the classification results, classified and unclassified reads, and a summary report.

##### run_k2r_sort_reads

This process runs Kraken2Ref to parse Kraken reports and sort reads by taxonomic classification. It generates JSON files that map taxonomic IDs to read IDs and performs sorting based on the decomposed taxonomy tree if available.

##### run_k2r_dump_fastqs_and_pre_report

Extract classified reads into FASTQ files and generate a preliminary report based on the taxonomic classification data. It processes classified reads and produces a detailed report for further analysis.

##### concatenate_fqs_parts
This process concatenates FASTQ files from multiple parts into final combined FASTQ files for each taxonomic classification. This process ensures that all parts corresponding to the same taxonomic ID are merged into single files.

##### get_taxid_references

Retrieves sequences for a given taxid from a source FASTA file and indexes them for further analysis.

##### bwa_alignment_and_post_processing

Mapping sequencing reads to a reference genome using BWA (Burrows-Wheeler Aligner), followed by post-processing steps including conversion to BAM format, sorting, and indexing. At current version of this pipeline, this process generates high-quality, sorted BAM files that are essential for consensus sequence generation.

##### run_ivar

Generates a consensus sequence from a sorted BAM file using samtools mpileup and ivar consensus.

##### run_pagolin

The process runs the Pangolin tool on a consensus FASTA file to determine the SARS-CoV-2 lineage and extracts relevant metadata from the output.

##### run_qc_script

This process runs a QC analysis on the input BAM, FASTA, and reference files, outputting QC metrics.

##### write_classification_report

This process creates a classification report from provided report lines, ensuring the output is correctly formatted and free from common formatting issues such as unexpected character encodings.

[**(&uarr;)**](#contents)

### Workflow

#### SORT_READS_BY_REF

The SORT_READS_BY_REF workflow processes paired-end sequencing reads by sorting them according to taxonomic classifications obtained from Kraken2. This workflow uses a manifest file to process multiple samples and produces sorted by taxid FASTQ files for each sample and classification reports.

#### GENERATE_CONSENSUS

The `GENERATE_CONSENSUS` workflow performs read alignment and consensus sequence generation for sequencing data. It processes paired-end reads by aligning them to reference genomes using BWA, followed by consensus calling with iVar. This workflow is designed to take in sequencing data for different samples and taxonomic IDs, process them, and produce consensus sequences.

#### COMPUTE_QC_METRICS

The `COMPUTE_QC_METRICS` workflow is designed to compute quality control (QC) metrics for consensus sequences generated from sequencing data. The workflow processes each sample's data to evaluate the quality and coverage of the generated consensus sequences. The QC metrics include the percentage of bases covered, the percentage of N bases, the longest segment without N bases, and read alignment statistics (total reads aligned, unmapped and mapped).

#### SCOV2_SUBTYPING

The `SCOV2_SUBTYPING` workflow is designed to determine the SARS-CoV-2 lineage (subtype) of consensus sequences using the PANGOLIN tool. This workflow takes in a channel of consensus sequences along with their metadata, runs the PANGOLIN lineage classification, and outputs updated metadata with the assigned lineage.

#### GENERATE_CLASSIFICATION_REPORT

The GENERATE_CLASSIFICATION_REPORT workflow generates a classification report based on metadata associated with sequencing samples. This workflow collects metadata from each sample, formats the data into a report line, and aggregates these lines into a final classification report file.

### Custom pipeline scripts

#### Kraken2ref JSON to TSV Report Script

The `kraken_report.py` script is designed to convert a Kraken2ref JSON output into a tab-separated values (TSV) report file. The script reads two main input files: the Kraken2ref JSON file and a corresponding Kraken2 taxonomic report. It then extracts relevant information about selected reference taxa, including virus subtypes and the number of reads per taxon, and generates a TSV report summarizing this data. The report provides detailed information on the sample, viruses, selected taxonomic IDs, flu segments, and subtyping data, specifically for influenza viruses if present.

**Key Features**

- **Taxonomic and Virus Information Extraction**: The script identifies the virus taxonomic IDs and names from both the Kraken2ref JSON and Kraken2 report, allowing detailed annotation of selected viruses and reference sequences.
- **Influenza Subtyping**: If influenza is detected, the script extracts subtype information (e.g., H and N subtypes) and flu segment numbers. These details are specifically captured for influenza viruses, helping to identify the sample subtype.
- **Reads Per Taxon**: The script reports the number of reads assigned to each selected taxonomic ID, providing insights into the abundance of each virus in the sample.

Customizable Output: Users can specify the output file suffix, giving flexibility in naming the report files.

**Example Command**

```bash
<path/to/viral_pipeline>/bin/kraken2ref_to_tsv.py -i kraken2ref_output.json -r kraken2_report.txt --out_suffix ".custom_report.tsv"
```

This command reads the JSON and report files, processes the data, and outputs a report named `<sample_id>.custom_report.tsv`.

#### QC Script for BAM and FASTA Files

This Python script generates a quality control (QC) summary report for a sample by analyzing a BAM file, a consensus FASTA file, a reference FASTA, and a per-position depth file. The QC metrics include the percentage of N bases in the consensus, the largest contiguous gap of N bases, and the percentage of reference bases covered at a minimum depth. The script also integrates alignment statistics from a SAMtools `flagstat` output to provide additional insights on the quality of the aligned reads.

**Inputs:**

1. **Consensus FASTA File** (`--fasta`): The consensus sequence file.

2. **Reference FASTA File** (`--ref`): The reference genome against which the reads were aligned and used for depth calculations.

3. **BAM File** (`--bam`): The aligned and filtered BAM file, which contains the reads aligned to the reference.

4. **Per-position Depths File** (`--depths_file`): A tab-delimited file listing read depth per position in the alignment.

5. **SAMtools Flagstat File** (`--flagstat_file`): The output from the `samtools flagstat` command, providing alignment statistics.

6. **Sample Name** (`--sample`): The name of the sample being processed.

7. **Output File** (`--outfile`): The path where the QC summary report (in CSV format) will be written.

8. **Minimum Depth** (`--minimum_depth`, optional): The minimum depth threshold used when calculating covered bases. Default is `10`.

9. **Ivar Minimum Depth** (`--ivar_md`, optional): The minimum depth used by ivar when generating the consensus FASTA.

**Outputs:**

The script generates a CSV file containing various QC metrics for the sample. The output columns include:

- `sample_name`: Name of the sample.
- `pct_N_bases`: Percentage of N bases in the consensus FASTA.
- `pct_covered_bases`: Percentage of the reference genome covered at or above the minimum depth threshold.
- `longest_no_N_run`: Length of the largest contiguous region without N bases in the consensus.
- `num_aligned_reads`: Number of aligned reads (from `SAMtools flagstat`).
- `bam`: Path to the input BAM file.
- `total_mapped_reads`: Total number of mapped reads (from `SAMtools flagstat`).
- `total_unmapped_reads`: Total number of unmapped reads (from `SAMtools flagstat`).
- `qc_pass`: Indicates whether the sample passed QC based on N content and gap size criteria.
- `ivar_md`: If applicable, the minimum depth used by ivar.

**Key Features:**

- **N Content Analysis**: The script calculates the percentage of bases that are 'N' in the consensus FASTA and finds the longest contiguous stretch of N bases.

- **Depth-Based Coverage**: It calculates the percentage of reference genome positions covered by reads at or above a specified depth threshold, using the per-position depths file.

- **Alignment Statistics**: The script reads the `SAMtools flagstat` file to report the number of aligned, mapped, and unmapped reads.

- **QC Pass/Fail**: The QC status is determined by evaluating N content (e.g., `TRUE` if the largest N gap exceeds `10,000` or if `N` bases constitute less than 50% of the sequence).

**Workflow:**

1. **QC Metric Calculation**:

- It calculates the percentage of `N` bases and identifies the longest `N` gap in the consensus sequence.
- It counts the positions in the depth file where the depth exceeds or equals the minimum threshold and calculates the coverage percentage.

2. **Flagstat Integration**: The script reads alignment statistics from the SAMtools flagstat file.

3. **QC Report Generation**: It assembles the calculated QC metrics into a dictionary, formats them, and writes them as a CSV report.

**Example Command**

```bash
<path/to/viral_pipeline>/bin/generate_qc.py --outfile sample123.qc.csv --sample sample123 --ref ref.fasta --bam sample123.bam --fasta sample123.consensus.fasta --depths_file sample123.depths.txt --flagstat_file sample123.flagstat.txt --minimum_depth 10 --ivar_md 5
```

This command processes the provided files and generates a QC summary report in `sample123.qc.csv`. The minimum depth for coverage calculations is set to 10, and ivar was run with a minimum depth of 5.

[**(&uarr;)**](#contents)

---

## Licence

[GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html)

[**(&uarr;)**](#contents)