# TRails input formats

This documents every input TRails reads, how to produce the ones you must supply, and which
inputs are public vs. licensed vs. your own data. **Only a few columns are required; almost
everything is optional — a missing optional column disables just the feature that depends on
it, never the whole run.** Column-name matching for the fixed/known columns is
**case-insensitive and underscore-insensitive** across every TSV input (`sample_id` ≡
`Sample_id` ≡ `SampleId` ≡ `SAMPLE_ID`); the per-sample column headers and `sample_id` *values*
still join exactly. Schemas were derived from the live analysis code.

## Data classes

| Class | Meaning | Examples | TRails behavior |
|-------|---------|----------|-----------------|
| **A** | Public, auto-fetchable | HPO `genes_to_phenotype.txt`; public repeat catalog; STRchive | downloaded into `reference_data/` by `install.sh` / `trails.py` |
| **B** | Licensed / controlled-access | gene–disease + phenotype-summary table | NEVER bundled; **optional** — pass `--gene-table` only if you are licensed |
| **C** | Your own cohort data | repeat-copy-numbers TSV; sample-metadata & phenotype TSVs; the generated result DB | NEVER bundled; you supply these |

TRails must start and run with class A only. Every class-B input is guarded so its absence
degrades gracefully (the dependent columns/annotations become NULL/disabled) rather than
crashing.

---

## Class C — the data you supply

### 1. Repeat copy numbers TSV (`--repeat-copy-numbers-tsv`)  *(required)*

Tab-separated (optionally gzipped). **One row per locus.** The first two columns identify the
locus; every remaining column is **one sample**, headed by that sample's `sample_id`. Each
cell holds the sample's per-allele repeat copy numbers, comma-separated — two values for a
diploid call, one for a hemizygous/haploid call; blank or `.` for a no-call.

```
trid                            motif   SAMPLE_A   SAMPLE_B   SAMPLE_C
chr1-57367043-57367119-AAAAT    AAAAT   12,12      11,13      12,40
chrX-148500631-148500691-GCC    GCC     21         20,21      45
```

| Column | Required? | Notes |
|--------|-----------|-------|
| `trid` | **required** | locus id `chrom-start-end-motif` (e.g. the TRGT trid); also used as `LocusId` |
| `motif` | **required** | repeat unit (e.g. `AAAAT`) |
| one column per sample | **≥1 required** | header = `sample_id`; cell = comma-separated per-allele repeat copy numbers |
| optional annotation columns | optional | recognized per-locus annotation columns are passed through to the loci table instead of being read as samples — see below |

**Optional per-locus annotation columns.** You may add extra columns to the matrix that
annotate the *locus* rather than a sample; they are recognized by name (case/underscore-insensitive)
and written to the loci table, not treated as samples. Recognized: `gene_id` /
`GencodeGeneId`, `gene_region` / `GencodeGeneRegion`, `ReferenceRegion`,
`NumRepeatsInReference`, and per-cohort population-stat columns (`HPRC256_*`, `AoU*`,
`TenK10K_*`, `TRExplorer*`). In particular, supplying a `gene_id` column is what lets the
optional `--gene-table` populate the gene-symbol / inheritance / pLI / disease-category columns
and enables gene-level phenotype scoring. Any header that is **not** a recognized annotation
column is treated as a sample.

This is a per-allele genotype matrix from **any** TR genotyping tool — short-read or long-read
(TRGT, ExpansionHunter, straglr, LongTR, …). TRails only needs the per-allele repeat copy
numbers, not a tool-specific VCF. (It generalizes the former TRGT-LPS table.) `trails.py` (via
`build_database.build()`) reads it **directly** into the result database — no intermediate
files are written; the per-locus records described below are built in memory.

### 2. Sample metadata TSV (`--sample-metadata-tsv`)  *(required)*

Tab-separated (optionally gzipped), **one row per sample**. **Only `sample_id` is required**;
every other column is optional and enables additional functionality when present.

```
sample_id      affected_status   analysis_status   family_id   maternal_id    paternal_id    sex      phenotype_description
PMGRC-1-1-0    Affected          Unsolved          FAM1        PMGRC-1-2-1    PMGRC-1-3-2    Male     Muscle weakness; Myopathy
PMGRC-1-2-1    Unaffected        Unaffected        FAM1                                      Female   Unaffected mother
```

| Column | Required? | Enables / values |
|--------|-----------|------------------|
| `sample_id` | **required** | joins to the sample columns of the repeat-copy-numbers TSV; may not contain `:` |
| `affected_status` | optional | affected-vs-unaffected comparison + affected-unsolved prioritization. `Affected` / `Possibly Affected` (treated as `Affected`) / `Unaffected` / `Unknown`. Absent → unknown |
| `analysis_status` | optional | solved/unsolved filtering + counts. `Solved` / `Unsolved` / `Unknown` / `Probably Solved` / `Partially Solved`. Absent → unknown |
| `family_id` | optional | distinct-family outlier counts |
| `maternal_id`, `paternal_id` | optional | Mendelian-violation QC (requires trios) |
| `sex` | optional | sex-aware handling of hemizygous loci. `Male` / `Female` / `Unknown` |
| `phenotype_description` | optional | free-text shown in the UI |

Any additional columns are accepted and ignored. (Status values are case-insensitive.)

### 3. Phenotype TSV (`--phenotypes-table`)  *(optional)*

Tab-separated (optionally gzipped), **one row per (sample, HPO term)**. Supplying it enables
phenotype-aware prioritization; omit it and TRails still runs (those scores are simply absent).

```
participant_id   term_id       hpo_description
PMGRC-1-1-0      HP:0000175    Cleft palate
PMGRC-1-1-0      HP:0001250    Seizure
```

| Column | Required? | Notes |
|--------|-----------|-------|
| `participant_id` | **required if file given** | must match `sample_id` |
| `term_id` | **required if file given** | HPO id, e.g. `HP:0000175` |
| `hpo_description` | optional | human-readable term; looked up from the HPO ontology if absent |

Extra columns (e.g. `age_of_onset`, `modifier`) are accepted and ignored.

### 4. Internal per-locus representation  *(in-memory; no file written; advanced)*

You do **not** author this, and TRails does **not** write it to disk — `build_database` builds
these per-locus records in memory directly from your TSV and populates the database. It is
documented here only as the internal data model, for contributors. Each record carries:

**Fields the build computes/reads:**

| Key | Required? | Format |
|-----|-----------|--------|
| `LocusId` | **required** | `chrom-start-end-motif` |
| `Motif` | **required** | repeat unit |
| `AllAlleleHistogram` | **required** | `allele:count` pairs, comma-separated, e.g. `9x:6,10x:4086,11x:178` |
| `ShortAlleleHistogram` | **required** | same format (shorter allele of each genotype); may be empty |
| `HemizygousAlleleHistogram` | **required** | same format; may be empty |
| `OutlierSampleIds_AllAlleles` | **required** | `allele:sample_id[:purity[:methylation]]` entries, comma-separated, **sorted by allele size DESC**, e.g. `11x:GSS225379,11x:PMGRC-111-107-2` |
| `OutlierSampleIds_ShortAlleles` | **required** | same format (short alleles) |
| `OutlierSampleIds_HemizygousAlleles` | **required** | same format (hemizygous) |
| `ReferenceRegion` | optional | `chrom:start-end`; computed from `LocusId` if absent |
| `NumRepeatsInReference` | optional | reference allele size in repeat units; computed if absent |

**Optional annotation keys** (absent ⇒ the matching output column is NULL): `CanonicalMotif`,
`Gencode*` (`GencodeGeneId`/`GencodeGeneRegion`/…), reference-quality fields, `TRExplorer*`,
RepeatMasker intervals, and per-cohort population-stat columns (`HPRC256_*`, `AoU1027_*`,
`TenK10K_*`). These come from an upstream annotation step; supply them only if you have them.
`VariantType` is **not** read by the pipeline (the server synthesizes it only for the
ExpansionHunter export) — optional/ignored.

### 5. Generated result database  *(server input; built for you)*

`trails.py` builds a single `*.with_analysis_columns.db` (SQLite) and points the server at it
(`results_server.py --db <path>`). You don't author it. It contains 7 tables:

| Table | Contents |
|-------|----------|
| `loci` | one row per locus (the analysis output columns + annotations) |
| `swim_plot` | one row per outlier allele (sample, allele size, affected status, motif, gene, …) |
| `per_locus_phenotype_scores` | gene/phenotype similarity per locus *(only if a phenotype TSV was supplied)* |
| `per_outlier_phenotype_scores` | gene/phenotype similarity per outlier sample *(only if a phenotype TSV was supplied)* |
| `sk_AllAlleles`, `sk_ShortAlleles`, `sk_HemizygousAlleles` | narrow "skinny" projections of `loci` for fast filtering/sorting (one per outlier type) |

The server also creates a small local `annotations.db` (your notes/tags) on first run.

### 6. Mendelian-violation QC  *(optional)*

If your sample-metadata TSV includes trio columns (`maternal_id` / `paternal_id`), TRails can
build a `mendelian_violations.db` (tables `mendelian_violations`,
`mendelian_violations_per_motif`) from the per-allele genotypes + family structure, enabling
the QC page. Absent → the QC page is hidden and the server still starts. Everything stays
local; nothing is uploaded.

---

## Class A — public reference data (fetched automatically)

| File | Flag | Source |
|------|------|--------|
| `genes_to_phenotype.txt` | `--genes-to-phenotype` | Human Phenotype Ontology annotation release (github.com/obophenotype/human-phenotype-ontology) |
| `variant_catalog_without_offtargets.GRCh38.json` | `--known-loci-json` | public `str-analysis` repo (github.com/broadinstitute/str-analysis) |
| STRchive disease-locus table | *(built-in default URL)* | github.com/dashnowlab/STRchive — cached locally by the installer; offline-safe |

---

## Class B — licensed / controlled-access (optional, user-supplied)

| Input | Flag | Restriction |
|-------|------|-------------|
| gene–disease + phenotype-summary table | `--gene-table` | OMIM/HPO-derived; the gene-symbol / inheritance / pLI / disease-category / phenotype-summary columns and gene-level phenotype scoring degrade without it |

If you don't pass this, TRails still runs — only the dependent columns/annotations are omitted.
