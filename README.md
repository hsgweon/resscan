# ResScan v.1.0.0

A comprehensive pipeline for identifying antimicrobial resistance (AMR) genes and variants from metagenomic sequencing data.

**ResScan** is a robust and user-friendly pipeline designed to process raw sequencing reads and provide a detailed, normalised report of AMR content. It leverages the CARD database and uses a dual-analysis approach:

1.  **HomScan:** Detects the presence and abundance of AMR genes based on homology.
2.  **VarScan:** Detects known resistance-conferring mutations (e.g., SNPs) in target genes.

The pipeline normalises results against universal single-copy genes (USCGs) and produces summary tables and rich, interactive HTML visualisations for data exploration.

## Features

-   **Dual Detection:** Simultaneously screens for both AMR gene presence (homology) and known resistance variants.
-   **Flexible Database Management:** Decouples the CARD database from the software, allowing users to easily download and build the latest version of the database. The SCG database is pre-packaged for convenience.
-   **Robust Normalisation:** Normalises AMR gene abundance using a suite of metrics (RPK, RPKG, RPKPC, FPK, etc.) for well-informed, multi-level sample comparisons.
-   **Advanced Ambiguity Resolution:** Implements a **Maximum A Posteriori (MAP)** iterative algorithm to statistically resolve ambiguous reads, identifying the most probable gene within a homologous family and reporting the proportion of evidence supporting each candidate.
-   **Prior-Guided Analysis:** Allows the incorporation of external knowledge (e.g., clinical prevalence data) via a user-provided priors file to improve the accuracy of the MAP resolution.
-   **Interactive Visualisations:** Automatically generates detailed HTML reports for each AMR gene family, showing read coverage, identity, and uniqueness.
-   **Safe & Smart Reruns:** Protects against accidental overwrites and allows for efficient re-analysis by skipping the time-consuming mapping step.
-   **Batch Processing:** Includes Nextflow-based parallel execution capabilities to process multiple samples simultaneously via a simple samplesheet.
-   **Result Aggregation:** Provides automated tools to pivot results from multiple samples into wide-format, zero-filled tables for downstream statistical analysis.
-   **Easy Installation:** Packaged for simple installation into a Conda environment.

## Installation

Installation is handled via Conda to ensure all dependencies are managed correctly.

**Prerequisites:**
*   You must have `conda` (or `miniconda`/`mamba`) installed.
*   You must have `git` installed to clone the repository.

Follow these three steps:

**1. Clone the Repository**: First, clone this repository to your local machine.

```bash
git clone https://github.com/hsgweon/resscan.git
cd resscan
```

**2. Create the Conda Environment**: Use the provided `environment.yml` file to create a self-contained environment with all necessary software (BWA, Diamond, Samtools, etc.). This command also uses `pip` to install the `resscan` scripts.

```bash
conda env create -f environment.yml
```

**3. Activate the Environment**: Activate the newly created environment. You must do this every time you want to use the pipeline.

```bash
conda activate resscan-env
```

## Database Preparation

Before running the pipeline, you must download the CARD database and prepare it using the included `resscan_build_db` tool.

**1. Download CARD Data**: Download the "Broadstreet" dataset from the [CARD Website](https://card.mcmaster.ca/).

```bash
# Example for CARD v4.0.1
wget https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2
mkdir -p CARD_raw_v4.0.1
tar -xvjf broadstreet-v4.0.1.tar.bz2 -C CARD_raw_v4.0.1
```

**2. Build the ResScan Database**: Run the builder script, pointing it to the folder where you extracted the files.

```bash
resscan_build_db -i CARD_raw_v4.0.1 -d resscan_CARD_v4.0.1
```

This will create a new directory `resscan_CARD_v4.0.1` containing the formatted FASTA and metadata files required by the pipeline.


**3. Audit Ignored Sequences**: The build process automatically filters CARD references that cannot be parsed or contain unsupported mutation types. By default, the tool produces three audit files in the output directory:

- `resscan_DB_CARD_metadata_ignored_non_variant.txt`: Details on failed homolog/knockout models.
- `resscan_DB_CARD_metadata_ignored_variant.txt`: Log of failed snps.txt entries.
- `ignored_sequences_report.html`: Interactive HTML report summarising exclusion.

## Optional: Metadata Curation

Gene family names in the raw database can be curated using `resscan_curate_metadata`.

Because ResScan utilises a hierarchical classification system, certain gene families with complex or non-standardised labeling—most notably the OXA beta-lactamases—can present challenges during assignment. The inherent complexity of these labels can result in a significant number of OXA-associated reads remaining unassigned to a specific family.

To improve taxonomic resolution and reporting accuracy, a level of manual metadata curation is highly recommended. This is particularly critical for the OXA gene family until such time as it undergoes formal reclassification within the underlying reference databases.

### 1. Create a Rules File (CSV)
Create a comma-separated file (e.g., an example file `curation_rules.csv`, is in scripts directory) defining renaming rules.

```csv
# Name_Substring,New_Family_Name
OXA-48,OXA-48-like Carbapenemase
OXA,Other OXA Beta-lactamase
NDM,NDM Metallo-beta-lactamase
```

### 2. Run the Curation Tool
```bash
# Backup metadata file just in case!
cp resscan_CARD_v4.0.1/resscan_DB_CARD_NR_metadata.txt resscan_CARD_v4.0.1/resscan_DB_CARD_NR_metadata.txt.backup

resscan_curate_metadata \
    -m resscan_CARD_v4.0.1/resscan_DB_CARD_NR_metadata.txt \
    -r scripts/curation_rules.csv
```

## Quick Start

To test your installation and run a basic analysis on paired-end reads using your prepared database and 16 threads. a test data (`test.fastq.gz`) is provided in `test_data` directory. Remember to use the correct PATH (can be either relative or absolute) of your CARD DB you have just built:

```bash
resscan -i test_data/test.fastq.gz \
        -o test_run \
        -t 16 \
        --card-db-dir resscan_CARD_v4.0.1
```
This will create a directory named `test_run` containing all the results. 

## Pipeline Workflow

The `resscan` command executes a series of scripts in a coordinated workflow:

1.  ***Mapping & QC (BWA & Diamond)***: Reads are mapped against CARD (BWA) and Universal Single-Copy Genes (DIAMOND).
2.  ***Normalisation Scaffolding (scgscan_*)***: USCG abundance is used to calculate average coverage for per-cell normalisation.
3.  ***VarScan (varscan_*)***: Processes BWA alignments to identify point mutations (SNPs) for rRNA and protein models. Requires all specified mutations for a model to be present on a single read. The pipeline applies a strict confirmation step:
    -   **Point Mutations Only:** Currently, VarScan detects single nucleotide polymorphisms (SNPs) for rRNA models and amino acid substitutions for protein models. **Frameshifts, deletions, and insertions are not yet implemented.**
    -   **Multi-Mutation Requirement:** If a resistance model defined by CARD requires multiple mutations (e.g., "Mutation A AND Mutation B"), VarScan requires *all* specified mutations to be present on a **single read**. If a read covers Mutation A but ends before reaching Mutation B, it is not counted. *Note: Using short-read data, this requirement is not always achievable if mutations are distant.*
    -   Confirmed variant hits are tabulated, normalised, and visualisations are generated showing the alignment and specific mutation sites.
4.  ***HomScan (homscan_*)***: The BWA alignments are processed again, this time to identify any read that maps to a homology-type AMR gene above a specified identity cutoff.
    -   **Tabulation & Normalisation:** All passing hits are aggregated. This step produces two initial reports: a _homscan.tsv file where ambiguous reads are grouped into 'multiple' categories.
    -   **MAP Resolution:** A final, more sophisticated Maximum A Posteriori (MAP) solver (homscan_resolve_MAP.py) is run. It uses an iterative algorithm to distribute the abundance from ambiguous reads among candidate genes based on the evidence from uniquely mapped reads and optional user-provided priors. This produces the final, most accurate quantitative report.
5.  ***Consolidation***: The final summary tables from all steps are copied from the temporary directory into the main output directory for easy access.

### Understanding the Normalisation Metrics

Raw read counts are not directly comparable between genes or samples due to variations in gene length and sequencing depth. To address this, **ResScan** calculates several normalised abundance metrics.

The pipeline uses two fundamental units for these calculations: **Reads** and **Fragments**.
* A **Read** is a single sequence from a FASTQ file.
* A **Fragment** represents the original piece of DNA sequenced. For paired-end data, the two reads (R1 and R2) from the same DNA fragment are counted as a single fragment. 

The primary metric for interpreting AMR abundance in this pipeline is **FPKPMC**, which estimates the copy number of a gene per million bacterial cells. This is the **signature metric** of the pipeline and the default used by the MAP algorithm for the most accurate abundance estimation (configurable).

| Metric | Calculation | Interpretation |
| :--- | :--- | :--- |
| **RPK / FPK** | Reads / (Gene Length / 1000) | **Reads/Fragments Per Kilobase.** Normalises for gene length. |
| **RPKG / FPKG** | RPK / (Total Sample Bases / 1e9) | **Reads/Fragments Per Kilobase per Gigabase.** Normalises for sequencing depth. |
| **RPKPC / FPKPC** | RPK_amr / RPK_uscg_avg | **Reads/Fragments Per Kilobase Per single-copy gene Copy.** Estimates average copy number per bacterium. |
| **RPKPMC / FPKPMC** | RPKPC * 1,000,000 | **RPKPC/FPKPC Multiplied by one Million.** Similar to "parts per million". |

## Understanding the MAP Resolver
The most challenging aspect of quantifying AMR genes from metagenomes is handling ambiguous reads—reads that map with high identity to multiple different but highly similar reference genes.

ResScan's ***Maximum A Posteriori (MAP)*** resolver offers a statistical solution to this problem by determining which gene is the most likely source of the ambiguous reads using an iterative expectation-maximisation-like algorithm.

### Incorporating External Knowledge with Priors
Sometimes, the evidence from unique reads is sparse. The MAP resolver can be guided by external knowledge using a priors file.

The `--map-priors-file` option takes a simple tab-separated file where you assign a numeric "weight" to specific AROs. An example priors file is provided in `scripts` directory.

```tsv
# Format: ARO_ID    Weight
ARO_3002312	15.0	# KPC-2
ARO_3000589	12.0	# NDM-1
```

## Usage

```bash
resscan -i <INPUT_FILES> -o <OUTPUT_PREFIX> --card-db-dir <DB_DIR> [OPTIONS]
```

### All Options

#### **Required Arguments**
| Flag | Description |
| :--- | :--- |
| `-i`, `--input-fastqs` | **Required.** Comma-separated list of input FASTQ files. They can be supplied in gzipped or unzipped FASTA or FASTQ format. |
| `-o`, `--output-prefix` | **Required.** A prefix for all output files and the main output directory. |
| `--card-db-dir` | **Required.** Path to the directory containing the prepared CARD database (output from `resscan_build_db`). |

#### **Database Arguments**
| Flag | Description |
| :--- | :--- |
| `--db-scg` | Path to a custom Single Copy Genes FASTA file. |
| `--db-scg-lengths` | Path to a custom SCG gene lengths TSV file. |

#### **Analysis Parameters**
| Flag | Description | Default |
| :--- | :--- | :--- |
| `--homscan-pid-cutoff` | Minimum percent identity for HomScan hits (0.0-1.0 scale). | `0.95` |
| `--varscan-pid-cutoff` | Minimum nucleotide percent identity for VarScan hits (0.0-1.0 scale). | `0.95` |
| `--homscan-pid-type` | PID type to use for HomScan filtering and WTA (`protein` or `nucleotide`). | `protein` |
| `--consensus-cutoff` | Minimum fraction of ambiguous hits that must map to the same gene family to reach a consensus. | `0.8` |
| `--homscan-gene-types` | Comma-separated list of gene types for HomScan (e.g., 'H,K'). | `H` |
| `--varscan-gene-types` | Comma-separated list of variant types for VarScan (e.g., 'V,R,O'). | `V,R` |

#### **MAP Resolver Arguments**
| Flag	| Description	| Default |
| :--- | :--- | :--- |
| `--map-priors-file` |	Path to a tab-separated file of priors. | None |
| `--map-metric-column` |	The numeric column to use for MAP abundance resolution.	| RPKG |
| `--map-base-prior` |	Baseline prior 'pseudo-count'.	| 1.0 |
| `--map-prior-strength` |	Multiplier for the influence of the priors file. | 1.0 |

#### **Performance & Control**
| Flag | Description |
| :--- | :--- |
| `-t`, `--threads` | Number of threads to use. |
| `--overwrite` | Overwrite the output directory if it already exists. |
| `--skip-mapping` | Skip the mapping steps. |
| `--skip-to-map-resolve` | Skip all steps and run only the final MAP resolution. |
| `--debug` | Enable debug mode. |
| `--version` | Show the pipeline version and exit. |

## Examples

**1. Standard Paired-End Analysis**
```bash
resscan -i sampleA_R1.fastq.gz,sampleA_R2.fastq.gz \
        -o SampleA_results \
        -t 32 \
        --card-db-dir ./resscan_CARD_v4.0.1
```

**2. Re-running only the MAP Resolver with Priors**
```bash
resscan -i sampleA_R1.fastq.gz,sampleA_R2.fastq.gz \
        -o SampleA_results \
        --card-db-dir ./resscan_CARD_v4.0.1 \
        --map-priors-file /path/to/clinical_priors.tsv \
        --skip-to-map-resolve \
        --overwrite
```

## Output Files

| File / Directory	| Description |
| :--- | :--- |
| `[prefix]_homscan_MAP.tsv` |	***(Primary Homology Result)*** The final, most accurate quantitative report for homology-based gene detection with MAP resolution. |
| `[prefix]_homscan.tsv`	| A "clean" summary grouped by family;multiple. |
| `[prefix]_varscan.tsv` |	***(Primary Variant Result)*** The final, normalised summary table for confirmed resistance variant detection. |
| `[prefix]_homscan_html/` |	Directory containing interactive HTML coverage plots. |
| `[prefix]_varscan_html/` |	Directory containing HTML alignment views. |
| `logs/` |	Run logs. |
| `tmp/` |	Intermediate files. |


## Batch Processing and Aggregation

ResScan includes high-level wrappers to manage large-scale studies and high-throughput sequencing projects efficiently.

### 1. Batch Execution (resscan_batch)
The `resscan_batch` script is a standalone Nextflow-powered orchestrator. It is designed to handle hundreds of samples with robust error recovery, parallel execution, and automated resource management, bypassing the limitations of standard Python multiprocessing.

**Key Features:**
-   **Nextflow Orchestration:** Operates as a dataflow-driven pipeline, ensuring that 100+ samples are processed without the "stalling" issues common in standard Python pools.
-   **Pre-flight Validation:** Automatically verifies every FASTQ path exists and checks for duplicate Sample IDs in your samplesheet before any analysis begins.
-   **Smart Resumption:** If a run is interrupted, simply add -resume (single dash) to the command. ResScan will skip successfully completed samples and only process the remaining ones.
-   **Variable Input Support:** Handles "ragged" CSVs where samples have a different number of FASTQ files (e.g., mixing single-end, paired-end, and technical replicates).
-   **Automatic Cleanup & Notifications:** Optional flags allow for automatic deletion of intermediate work/ directories and email notifications upon completion.


**Samplesheet Format (`samplesheet.csv`)**

The samplesheet is a CSV file. The first column must be named `sample`. Subsequent columns contain paths to FASTQ files.
```tsv
#sample,fastq_1,fastq_2,fastq_3
SampleA,/data/A_R1.fq.gz,/data/A_R2.fq.gz
SampleB,/data/B_R1.fq.gz,/data/B_R2.fq.gz,/data/B_R1_L002.fq.gz
SampleC,/data/C_SE.fq.gz
```

**Run Command:**
```bash
# Run 15 samples in parallel, using 8 threads per sample
./resscan_batch \
    --samplesheet data/samplesheet.csv \
    --card resscan_CARD_v4.0.1 \
    --out project_results \
    --parallel 15 \
    --threads 8
```

Additional Batch Options:
| Flag | Description |
| :--- | :--- |
| --dry_run | Test the pipeline logic (sleeps 2s per sample) without running ResScan. |
| --cleanup | Deletes the Nextflow work/ directory upon a successful run to save disk space. |
| --args | Pass extra flags to the core ResScan tool (e.g., --args "--min-id 0.98"). |


### 2. Result Aggregation (resscan_aggregate)
Once batch processing is complete, use the aggregator to consolidate individual sample results into wide-format matrices suitable for downstream statistical analysis in R or Python.

```bash
# Aggregate all result types (homscan, varscan, map) into a new directory
resscan_aggregate -i ./resscan_batch_results -o Project_Summary -p Project_Summary
```

- `i`: The parent directory containing all individual sample folders.
- `o`: The output directory where the aggregated tables will be saved (created automatically).
- `p`: (Optional) A prefix for the filenames (Default: aggregated).

The aggregator generates separate pivoted TSV files for each metric (e.g., Project_Summary_homscan_RPK.tsv). These tables use Sample IDs as columns and Gene metadata as rows. Missing values (genes not detected in a sample) are automatically filled with 0.

**Output Data:**

The aggregator generates pivoted TSV files for each metric (e.g., `Project_Summary_homscan_RPK.tsv`).

- Rows: Resistance Genes/Alleles.
- Columns: Sample IDs.
- Values: Normalised metrics (RPK, TPM, or Depth).
- Fill: Missing values (genes not detected) are automatically filled with 0.

### 3. Post-Run Housekeeping

After a batch run has completed successfully and you have verified your results in the output directory, you can safely remove the temporary files generated by Nextflow to reclaim disk space.

Note: Deleting these files will remove the "checkpoint" data, meaning you will not be able to use the -resume feature for that specific run.

**Manual Cleanup**

You can delete the hidden Nextflow metadata and the intermediate work/ directory using the following command:
```bash
# Safely remove Nextflow hidden files and the intermediate work directory
rm -rf .nextflow* work/
```

**Automated Cleanup**
If you want the pipeline to handle this automatically upon success, you can use the `--cleanup` flag when running the script:
```bash
# Safely remove Nextflow hidden files and the intermediate work directory
resscan_batch --samplesheet data.csv --card ./db --cleanup
```

### 4. HPC and Workflow Orchestration
Because `resscan_batch` is built on Nextflow, it is natively compatible with HPC schedulers (Slurm, SGE, LSF) and Cloud environments. While the default configuration runs locally, the script can be scaled to thousands of samples across a distributed cluster by defining a Nextflow profile.

**For Bioinformaticians & HPC Users:**
For massive-scale datasets (thousands of samples) or execution on High-Performance Computing (HPC) clusters using schedulers like SLURM or SGE, we recommend wrapping the core resscan command in a workflow manager such as Nextflow or Snakemake. This provides superior error recovery (check-pointing), containerisation (Docker/Singularity), and cloud-native scaling.


## License
This project is licensed under the MIT License.

## Contact
For questions, bug reports, or suggestions, please open an issue on this GitHub repository.
