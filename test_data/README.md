# ResScan test data

A small, deliberately AMR-enriched dataset for checking that an installation works
end to end, and for demonstrating every supported samplesheet layout.

## What this is

Four paired-end libraries, roughly 2,700–3,600 read pairs each (1.8 MB in total),
subsampled from a published ICU sink-drain metagenomic study:

| Sample | Read pairs | HomScan families | VarScan rows |
| :--- | ---: | ---: | ---: |
| `SRR10842841` | 2,831 | – | 5 (at full depth) |
| `SRR10842842` | 2,660 | – | 4 (at full depth) |
| `SRR10842855` | 2,831 | 55 (single-end) | 2 (single-end) |
| `SRR10842857` | 3,562 | 90 (paired-end) | 5 (paired-end) |

Reads are 150 bp, already quality-trimmed and host-depleted, so they can be passed
to `resscan` directly.

> **These are not representative abundance data — do not interpret the numbers.**
> A random subsample of this size would contain almost no AMR reads at all (a plain
> 10,000-read subsample of these libraries yields between 0 and 4 gene rows). To make
> a useful test, reads known to align to CARD were deliberately **over-represented**
> relative to the single-copy marker genes that form the per-cell denominator. Every
> abundance ResScan reports from this dataset is therefore inflated by a large and
> arbitrary factor. The data are here to prove the pipeline runs and produces
> sensible-looking output, nothing more.

Each file contains four kinds of read: pairs that align to CARD, pairs that confirm a
resistance-conferring point mutation (so VarScan has something to report), pairs that
align to the universal single-copy genes (so per-cell normalisation has a real
denominator), and a small number of background pairs. The samples were chosen because
they are the richest in variant-confirming reads; most libraries in the source study
contain none at all.

## Running the test

From the repository root, with a prepared CARD database (see the main README):

```bash
resscan -i test_data/fastq/SRR10842857_1.fastq.gz,test_data/fastq/SRR10842857_2.fastq.gz \
        -o test_results --card /path/to/resscan_CARD_v4.0.1 -t 4
```

Or run the whole samplesheet, which exercises all four input layouts at once:

```bash
resscan_batch --samplesheet test_data/samplesheet.csv \
              --card /path/to/resscan_CARD_v4.0.1 \
              --out test_batch_results --parallel 2 --threads 4
```

## What the samplesheet covers

`samplesheet.csv` defines four entries, one per supported input layout. In the file
list, `,` separates the mates of a single run and `;` separates independent runs of
the same sample.

| Entry | Layout | Exercises |
| :--- | :--- | :--- |
| `TEST_PE` | paired-end, one run | the usual case; fragment counting with two mates |
| `TEST_SE` | single-end, one run | single-mate handling — fragment counts equal read counts |
| `TEST_PE_MULTIRUN` | paired-end, two runs | merging a sample sequenced across several runs |
| `TEST_SE_MULTIRUN` | single-end, two runs | multi-run merging without pairing |

All runs of a sample must have the same number of mates; mixing paired and single-end
runs within one sample is rejected with an error.

## Expected outcome

The run should complete without errors and produce a non-empty `*_homscan.tsv`
together with per-cell (`RPKPMC`) values. Observed on this dataset with CARD v4.0.1:

| Entry | HomScan families | VarScan rows |
| :--- | ---: | ---: |
| `TEST_PE` (SRR10842857, paired-end) | 90 | 5 |
| `TEST_SE` (SRR10842855, single-end) | 55 | 2 |

The VarScan calls are real determinants — *Escherichia coli* UhpT mutations conferring
fosfomycin resistance (105 supporting reads in `TEST_PE`), plus 16S and 23S rRNA
variants. VarScan counts are small by nature: resistance-conferring point mutations are
rare, and a read is only counted when it spans and confirms the mutation itself.

Two things worth checking in the output, because they confirm the counting behaves
correctly:

-   `Effective_Length_bp` equals `Gene_Length_bp` for every gene, and the run log
    reports an effective-length offset of `+0.00 bp`. With 150 bp reads and the default
    `--homscan-min-aln-frac 0.5`, per-kilobase abundance is unbiased by construction, so
    no correction is applied. (A trailing `-0.01 bp` may appear where a few trimmed
    reads have odd lengths; this is negligible.)
-   For the single-end entries, `Fragment_Count` equals `Read_Count` for every gene. A
    single-end library has no fragment length, so the two units necessarily coincide.

## Provenance

Subsampled from public data (ENA accessions above), originally from a study of
antimicrobial resistance in hospital sink drains. Redistributed here in reduced form
purely as installation test material.
