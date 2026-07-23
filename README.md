# ResScan v.1.0.1

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

## Requirements

**Operating system:** Linux or macOS (the pipeline shells out to command-line aligners; on Windows use WSL2).

**Software** (installed automatically by the Conda environment):
-   Python ≥ 3.10
-   BWA — AMR read mapping
-   DIAMOND — single-copy-gene mapping
-   Samtools
-   pandas, numpy
-   Nextflow ≥ 26.0 — required only for `resscan_batch` (needs Java 11+)

**Hardware and resources.** ResScan's own databases are small — the prepared CARD reference and the packaged 40-gene SCG database are sub-gigabyte — so the core mapping steps (BWA, DIAMOND) are modest in memory, typically a few GB. The largest disk consumer is the per-sample intermediate SAM/BAM files written under the Nextflow `work/` directory; these scale with sequencing depth and can be reclaimed after a run with `--cleanup` (see Post-Run Housekeeping). Runtime scales with the number of reads per sample and the threads allocated (`-t` / `--threads`); the batch wrapper additionally parallelises across samples (`--parallel`). Note that upstream host-read removal (not part of ResScan) is usually the most resource-intensive step — a human `bowtie2` index alone needs several GB of RAM and disk.

## Installation

Installation is handled via Conda to ensure all dependencies are managed correctly.

**Prerequisites:**
*   You must have `mamba` (recommended) or `conda` installed. Mamba is significantly faster at solving environments and is available via [Miniforge](https://github.com/conda-forge/miniforge). If you only have `conda`, replace `mamba` with `conda` in the commands below.
*   You must have `git` installed to clone the repository.

Follow these three steps:

**1. Clone the Repository**: First, clone this repository to your local machine.

```bash
git clone https://github.com/hsgweon/resscan.git
cd resscan
```

**2. Create the Environment**: Use the provided `environment.yml` file to create a self-contained environment with all necessary software (BWA, Diamond, Samtools, etc.). This command also uses `pip` to install the `resscan` scripts.

```bash
mamba env create -f environment.yml
```

**3. Activate the Environment**: Activate the newly created environment. You must do this every time you want to use the pipeline.

```bash
mamba activate resscan-env
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

## Input Requirements

ResScan quantifies AMR genes from sequencing reads; it does **not** perform read QC or host-read removal itself. Prepare your reads before running:

1. **Quality control** — remove adapters and low-quality bases (e.g. `fastp` or `Trimmomatic`).
2. **Host-read removal** — for clinical or host-associated samples, remove host (e.g. human) reads by mapping against the host genome with `bowtie2` or `bwa` and discarding mapped reads.

Both steps can be done in a single pass with [`vanish`](https://github.com/hsgweon/vanish), which performs read QC and host-read removal together.

**On host DNA and the per-cell metric.** FPKPMC is the ratio of AMR fragment density to the density of bacterial single-copy genes. Non-bacterial (host) DNA maps to neither term, so it does not bias the per-cell estimate — this is what makes FPKPMC robust to *variable host contamination* across samples, unlike relative metrics such as RPKM. Host-read removal is nonetheless recommended upstream: it reduces compute time, avoids occasional spurious mapping, and is frequently required for data-governance reasons when human reads are present. In short, the per-cell normalisation handles *residual* host DNA gracefully, but it is not a substitute for proper host depletion.

## Quick Start

To check your installation, a small test dataset (`test.fastq.gz`) is provided in the `test_data` directory. The command below runs a basic analysis using your prepared database and 16 threads. Use the correct path (relative or absolute) to the CARD database you built:

```bash
cd test_data
resscan -i test.fastq.gz \
        -o test_run \
        -t 16 \
        --card ../resscan_CARD_v4.0.1
```
This will create a directory named `test_run` containing all the results.

**Input format.** `-i` uses commas to separate mates within a run and semicolons to separate multiple sequencing runs from the same sample. ResScan derives the fragment count from the mate structure automatically.

Single-end:
```bash
resscan -i sample.fastq.gz -o sample -t 16 --card resscan_CARD_v4.0.1
```

Paired-end (one run):
```bash
resscan -i sample_R1.fastq.gz,sample_R2.fastq.gz \
        -o sample \
        -t 16 \
        --card resscan_CARD_v4.0.1
```

Paired-end (multiple runs — mates within each run separated by commas, runs separated by semicolons):
```bash
resscan -i "run1_R1.fastq.gz,run1_R2.fastq.gz;run2_R1.fastq.gz,run2_R2.fastq.gz" \
        -o sample \
        -t 16 \
        --card resscan_CARD_v4.0.1
```

All runs must have the same number of mates; ResScan will exit with an error if they differ.

## Pipeline Workflow

The `resscan` command executes a series of scripts in a coordinated workflow:

1.  ***Mapping & QC (BWA & Diamond)***: Reads are mapped against CARD (BWA) and Universal Single-Copy Genes (DIAMOND).
2.  ***Normalisation Scaffolding (scgscan_*)***: USCG abundance is used to calculate average coverage for per-cell normalisation.
3.  ***VarScan (varscan_*)***: Processes BWA alignments to identify point mutations (SNPs) for rRNA and protein models. Requires all specified mutations for a model to be present on a single read. The pipeline applies a strict confirmation step:
    -   **Point Mutations Only:** Currently, VarScan detects single nucleotide polymorphisms (SNPs) for rRNA models and amino acid substitutions for protein models. **Frameshifts, deletions, and insertions are not yet implemented.**
    -   **Multi-Mutation Requirement:** If a resistance model defined by CARD requires multiple mutations (e.g., "Mutation A AND Mutation B"), VarScan requires *all* specified mutations to be present on a **single read**. If a read covers Mutation A but ends before reaching Mutation B, it is not counted. *Note: Using short-read data, this requirement is not always achievable if mutations are distant.*
    -   Confirmed variant hits are tabulated, normalised, and visualisations are generated showing the alignment and specific mutation sites.
4.  ***HomScan (homscan_*)***: The BWA alignments are processed again, this time to identify any read that maps to a homology-type AMR gene above a specified identity cutoff.
    -   **Tabulation, Normalisation & MAP Resolution:** All passing hits are aggregated. This step produces two reports: a `_homscan_detailed.tsv` showing each distinct ambiguous mapping pattern explicitly, and a `_homscan.tsv` where the MAP solver has been applied. The MAP solver uses an iterative algorithm to distribute the abundance from ambiguous reads among candidate genes based on the evidence from uniquely mapped reads and optional user-provided priors. The result is reported as `Top_ARO` and `Allocation_Proportions` columns directly in `_homscan.tsv`.
5.  ***Consolidation***: The final summary tables from all steps are copied from the temporary directory into the main output directory for easy access.

### Understanding the Normalisation Metrics

Raw read counts are not directly comparable between genes or samples due to variations in gene length and sequencing depth. To address this, **ResScan** calculates several normalised abundance metrics.

The pipeline uses two fundamental units for these calculations: **Reads** and **Fragments**.
* A **Read** is a single sequence from a FASTQ file.
* A **Fragment** represents the original piece of DNA sequenced. For paired-end data, the two reads (R1 and R2) from the same DNA fragment are counted as a single fragment. 

The primary metric for interpreting AMR abundance in this pipeline is **FPKPMC**, which estimates the copy number of a gene per million bacterial cells. This is the **signature metric** of the pipeline and the default used by the MAP algorithm for the most accurate abundance estimation (configurable).

| Metric | Calculation | Interpretation |
| :--- | :--- | :--- |
| **RPK / FPK** | Reads (or Fragments) / (Gene Length / 1000) | **Reads/Fragments Per Kilobase.** Normalises for gene length. |
| **RPKG / FPKG** | RPK (or FPK) / (Total Sample Bases / 1e9) | **Reads/Fragments Per Kilobase per Gigabase.** Normalises for sequencing depth. |
| **RPKM** | RPK / (Total Reads / 1e6) | **Reads Per Kilobase per Million reads.** Total reads counts every individual read including all mates. |
| **FPKM** | FPK / (Total Fragments / 1e6) | **Fragments Per Kilobase per Million fragments.** Denominator is fragment count (paired mates counted once). |
| **RPKPC / FPKPC** | RPK_amr / RPK_uscg_avg | **Reads/Fragments Per Kilobase Per Cell.** Estimated average gene copies per bacterial cell. |
| **RPKPMC / FPKPMC** | RPKPC × 1,000,000 | **Reads/Fragments Per Kilobase Per Million Cells.** FPKPMC of 1,000,000 ≈ one gene copy per cell. |
| **RPKPGC / FPKPGC** | RPKPC × 1,000,000,000 | **Reads/Fragments Per Kilobase Per Billion Cells.** FPKPGC scaled by 10⁹; useful for very low-abundance genes. |

**Single-copy-gene normalisation (the per-cell denominator).** The `RPK_uscg_avg` / `FPK_uscg_avg` term used by the per-cell metrics is the average density of a packaged set of **40 universal single-copy marker genes** — ribosomal proteins, RNA polymerase subunits, aminoacyl-tRNA synthetases, and other core housekeeping genes. Sequences are derived from the **NCBI COG2020** database and quantified at the protein level with **DIAMOND**. Because these genes occur in (almost) exactly one copy per genome across bacteria, their average density approximates the number of bacterial genomes — i.e. cells — in the sample, which is what makes the per-cell metrics biologically interpretable. A custom marker set can be supplied via `--db-scg` / `--db-scg-lengths`.

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
resscan -i <INPUT_FILES> -o <OUTPUT_PREFIX> --card <DB_DIR> [OPTIONS]
```

### All Options

#### **Required Arguments**
| Flag | Description |
| :--- | :--- |
| `-i`, `--input-fastqs` | **Required.** Input FASTQ files. Mates within a run are comma-separated; multiple runs from the same sample are semicolon-separated (e.g. `run1_R1.fq,run1_R2.fq;run2_R1.fq,run2_R2.fq`). Files may be gzipped or unzipped FASTQ/FASTA. |
| `-o`, `--output-prefix` | **Required.** A prefix for all output files and the main output directory. |
| `--card` | **Required.** Path to the directory containing the prepared CARD database (output from `resscan_build_db`). |

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
| `--homscan-pid-type` | PID type to use for HomScan filtering (`protein` or `nucleotide`). | `protein` |
| `--consensus-cutoff` | Minimum fraction of ambiguous hits that must map to the same gene family to reach a consensus. | `0.8` |
| `--homscan-gene-types` | Comma-separated list of gene types for HomScan (e.g., 'H,K'). | `H` |
| `--varscan-gene-types` | Comma-separated list of variant types for VarScan (e.g., 'V,R,O'). | `V,R` |

**CARD model types (for `--homscan-gene-types` and `--varscan-gene-types`).** Each reference in CARD carries a detection model; these single-letter codes select which references an engine considers:

| Code | CARD model type | Resistance mechanism | Engine |
| :--- | :--- | :--- | :--- |
| `H` | Protein homolog model | Gene presence confers resistance | HomScan (default) |
| `K` | Protein knockout model | Gene loss/disruption confers resistance | HomScan |
| `V` | Protein variant model | A specific amino-acid substitution confers resistance | VarScan (default) |
| `R` | rRNA gene variant model | A specific rRNA mutation confers resistance | VarScan (default) |
| `O` | Protein overexpression model | A regulatory mutation raising expression confers resistance | VarScan |

By default HomScan uses `H` and VarScan uses `V,R`. Give homology-type models (`H`, `K`) to HomScan and variant-type models (`V`, `R`, `O`) to VarScan; pairing an engine with a model type it isn't designed for (e.g. giving VarScan a homolog model that has no defined mutations) will not produce meaningful results.

#### **MAP Resolver Arguments**
| Flag	| Description	| Default |
| :--- | :--- | :--- |
| `--map-priors-file` |	Path to a tab-separated file of priors. | None |
| `--map-metric-column` |	The numeric column to use for MAP abundance resolution.	| FPKPMC |
| `--map-base-prior` |	Baseline prior 'pseudo-count'.	| 1.0 |
| `--map-prior-strength` |	Multiplier for the influence of the priors file. | 1.0 |

#### **Performance & Control**
| Flag | Description |
| :--- | :--- |
| `-t`, `--threads` | Number of threads to use. |
| `--overwrite` | Overwrite the output directory if it already exists. |
| `--skip-mapping` | Skip the mapping steps. |
| `--skip-to-tabulate` | Skip mapping and SAM processing; re-run tabulation (including MAP) and consolidation only. Useful for re-running with different MAP priors. |
| `--debug` | Enable debug mode. |
| `--version` | Show the pipeline version and exit. |

## Examples

**1. Single-end analysis**
```bash
resscan -i sample.fastq.gz \
        -o SampleA_results \
        -t 32 \
        --card ./resscan_CARD_v4.0.1
```

**2. Standard paired-end analysis**
```bash
resscan -i sampleA_R1.fastq.gz,sampleA_R2.fastq.gz \
        -o SampleA_results \
        -t 32 \
        --card ./resscan_CARD_v4.0.1
```

**3. Paired-end analysis with multiple sequencing runs from the same sample**

Mates within each run are comma-separated; runs are semicolon-separated.
```bash
resscan -i "run1_R1.fastq.gz,run1_R2.fastq.gz;run2_R1.fastq.gz,run2_R2.fastq.gz" \
        -o SampleA_results \
        -t 32 \
        --card ./resscan_CARD_v4.0.1
```

**4. Paired-end analysis with MAP priors**
```bash
resscan -i sampleA_R1.fastq.gz,sampleA_R2.fastq.gz \
        -o SampleA_results \
        -t 32 \
        --card ./resscan_CARD_v4.0.1 \
        --map-priors-file /path/to/clinical_priors.tsv
```

**5. Re-running with different MAP priors**
```bash
resscan -i sampleA_R1.fastq.gz,sampleA_R2.fastq.gz \
        -o SampleA_results \
        --card ./resscan_CARD_v4.0.1 \
        --map-priors-file /path/to/clinical_priors.tsv \
        --skip-to-tabulate \
        --overwrite
```

## Output Files

| File / Directory | Description |
| :--- | :--- |
| `[prefix]_homscan.tsv` | ***(Primary Homology Result)*** The homology-based quantitative report. One row per gene (ARO) or per gene family (`multiple` rows for ambiguous reads). Includes all normalised metrics, inline annotation (ARO_Name, Drug_Class, Resistance_Mechanism), and MAP resolution columns (`Top_ARO`, `Allocation_Proportions`). |
| `[prefix]_homscan_detailed.tsv` | ***(Detailed Homology Report)*** The same data before ambiguous rows are collapsed, showing each distinct `multiple[ARO_x\|ARO_y]` combination explicitly. Includes gene length and all normalised metrics with inline annotation. |
| `[prefix]_varscan.tsv` | ***(Primary Variant Result)*** The final, normalised summary table for confirmed resistance variant detection, with inline annotation. |
| `[prefix]_homscan_html/` | Directory containing interactive HTML coverage plots. |
| `[prefix]_varscan_html/` | Directory containing HTML alignment views. |
| `logs/` | Run logs. |
| `tmp/` | Intermediate files. |


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


**Samplesheet Format (`samplesheet.csv`):**

The batch pipeline uses a flexible, positional CSV format. Unlike many pipelines, it does not require specific header names; it identifies data based on column order.

**Format Rules:**
-   **Column 1:** Must be the Sample ID (used for naming folders and files).
-   **Column 2 onwards:** Every subsequent column is treated as a path to a FASTQ file (supports single-end, paired-end, or multiple lanes).
-   **Comments:** Any line starting with a hash (#) is completely ignored.
-   **Headers:** A header row is not required. If you include one, you must start the line with a #.

**Example `samplesheet.csv`:**
```tsv
# sample_id, read_1, read_2, read_3 (Optional Header - Ignored due to #). Spaces and tabs are ignored.
SampleA, data/A_R1.fq.gz, data/A_R2.fq.gz
SampleB, data/B_L001_R1.fq.gz, data/B_L001_R2.fq.gz, data/B_L002_R1.fq.gz, data/B_L002_R2.fq.gz
SampleC, data/C_single_end.fastq

# You can comment out specific samples to skip them
# SampleD, data/D_R1.fq.gz, data/D_R2.fq.gz
```


**Building the samplesheet automatically (optional):**

If you have a mapping file that links sample IDs to file prefixes (a two-column CSV with no header: `sample_id,file_prefix`), the `create_samplesheet_from_mapping` helper can generate the samplesheet in one step. Multiple rows sharing the same `sample_id` are treated as separate sequencing runs from that sample.

```bash
create_samplesheet_from_mapping \
    --mapping sampleid2fileprefix.csv \
    --fastq-dir rawdata/ \
    --out samplesheet.csv
```

The script scans the directory for files whose names start with each prefix, sorts them (so R1 sorts before R2), and writes one row per sample with mates comma-separated and runs semicolon-separated. All paths are written as absolute paths. A per-sample summary is printed to stdout showing how many runs and mates were found.

**Example `sampleid2fileprefix.csv`:**
```
Eff,SRR8830802
Eff,SRR8830803
Pig,SRR8830792
```

**Resulting `samplesheet.csv`:**
```
Eff, /abs/rawdata/SRR8830802_1.fastq.gz,/abs/rawdata/SRR8830802_2.fastq.gz;/abs/rawdata/SRR8830803_1.fastq.gz,/abs/rawdata/SRR8830803_2.fastq.gz
Pig, /abs/rawdata/SRR8830792_1.fastq.gz,/abs/rawdata/SRR8830792_2.fastq.gz
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
Once batch processing is complete, use the aggregator to consolidate individual sample results into wide-format, zero-filled tables. This tool uses your samplesheet as the "source of truth" to ensure all samples are accounted for.

**Key Features:**
-   **Samplesheet-Driven:** Uses the sample column from your CSV to find files and name columns.
-   **Double-Nested Path Awareness:** Automatically navigates the sample_id/sample_id_results/ structure created by the batch script.
-   **Zero-Filling:** Samples with no detected AMR genes (or missing files) are explicitly included as columns of zeros.
-   **Explicit Reporting:** Provides a terminal summary identifying exactly which samples were found, empty, or missing.

```bash
# Aggregate all result types (homscan, varscan) into a new directory
resscan_aggregate -s samplesheet.csv -i project_results -o project_summary -p Project_Name
```

The aggregator generates separate pivoted TSV files for each metric (e.g., `Project_Name_homscan_FPKPMC.tsv`). These tables use Sample IDs as columns and Gene metadata as rows. Missing values (genes not detected in a sample) are automatically filled with 0.

**Output Data:**

- Rows: Resistance Genes/Alleles (identified by AMR_Gene_Family, ARO, ARO_Name, Drug_Class, Resistance_Mechanism).
- Columns: Sample IDs.
- Values: One table per metric from: `RPK`, `FPK`, `RPKG`, `FPKG`, `RPKM`, `FPKM`, `RPKPC`, `FPKPC`, `RPKPMC`, `FPKPMC`, `RPKPGC`, `FPKPGC` (plus `Read_Count`, `Fragment_Count`). Only metrics present in the per-sample reports are emitted.

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


## Limitations

ResScan's variant detection (VarScan) is deliberately conservative. Two constraints should be borne in mind when interpreting results:

-   **Substitution variants only.** VarScan reports resistance conferred by nucleotide and amino-acid substitutions. AROs whose CARD model requires a frameshift, insertion, deletion, or duplication are skipped (and reported during the run, e.g. "Skipped N AROs containing complex mutations"). Such determinants will not appear in the output.
-   **Single-read allelic phasing.** For multi-residue models, VarScan requires *all* constituent mutations to be present on a single read alignment. This guarantees physical linkage and removes false "inferred resistance" arising when mutations are distributed across different molecules — but it also means that, with short reads, mutations separated by more than the read length cannot be co-confirmed and may be missed. VarScan favours specificity over sensitivity by design. (In paired-end data each mate is a separate alignment, so the required mutations must co-occur within a single mate.)

These are intentional trade-offs in favour of precision; the per-variant HTML evidence reports let any borderline call be inspected directly.

## License
This project is licensed under the MIT License.

## Contact
For questions, bug reports, or suggestions, please open an issue on this GitHub repository.
