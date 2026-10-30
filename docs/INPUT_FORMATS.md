# TRails input formats

This documents every input TRails reads, how to produce the ones you must supply, and which
inputs TRails downloads for you vs. which are your own data. **Only a few columns are required; almost
everything is optional — a missing optional column disables just the feature that depends on
it, never the whole run.** Column-name matching for the fixed/known columns is
**case-insensitive and underscore-insensitive** across every TSV input (`sample_id` ≡
`Sample_id` ≡ `SampleId` ≡ `SAMPLE_ID`); the per-sample column headers and `sample_id` *values*
still join exactly. Schemas were derived from the live analysis code.

## Data classes

| Class | Meaning | Examples | TRails behavior |
|-------|---------|----------|-----------------|
| **A** | Public, auto-fetchable | TRExplorer population statistics and locus annotations; genelookup gene table; HPO `genes_to_phenotype.txt`; known-disease-locus catalogs | downloaded into `reference_data/` by `install.sh` / `trails.py` |
| **C** | Your own cohort data | repeat-copy-numbers TSV; sample-metadata & phenotype TSVs; the generated result DB | NEVER bundled; you supply these |

Your class C tables plus the class A downloads are everything TRails needs. A class A download
that fails (e.g. offline) only leaves the columns that depend on it empty; the build still runs.
All coordinates are GRCh38 (hg38).

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
| `trid` | **required** | locus id `chrom-start-end-motif` (e.g. the TRGT trid); also used as `LocusId`. Exactly four `-`-separated fields, with integer start/end; any other shape fails the build with an error naming the column, the value and this file |
| `motif` | **required** | repeat unit (e.g. `AAAAT`). Only IUPAC bases (`ACGTNYRSWMKBVDH`, case-insensitive); a decorated or placeholder motif such as `(CAG)n` fails the build with an error naming the locus, the column and the value |
| one column per sample | **≥1 required** | header = `sample_id`; cell = comma-separated per-allele repeat copy numbers |
| optional annotation columns | optional | recognized per-locus annotation columns annotate the locus instead of being read as samples; most are written to the loci table, a few are always recomputed — see below |

**Annotations from TRExplorer.** You do not need to supply any annotation: each locus gets the
population statistics (`HPRC256_*`, `AoU1027_*`, `TenK10K_*`), Gencode gene id and region,
variation cluster and reference-repeat fields of its locus in TRExplorer's repeat catalog
(class A, below). A locus matches the TRExplorer locus with the same id, ignoring a `chr`
prefix; otherwise the TRExplorer locus whose interval overlaps it with a Jaccard similarity
(intersection over union) of at least 0.66 and whose motif matches (the same canonical motif for
1-6 bp motifs, the same motif length for 7 bp or longer), preferring the largest HPRC256, then
AoU1027, standard deviation. `TRExplorerLocusId` and `TRExplorerLocusJaccardSimilarity` record
the match. A locus with no match keeps these columns empty. The build prints how many loci
matched by id, by overlap, or not at all.

**Optional per-locus annotation columns.** You may add extra columns to the matrix that
annotate the *locus* rather than a sample; they are recognized by name (case-, underscore- and
space-insensitive) and are not treated as samples. Any header that is **not** recognized by one
of the two rules below is treated as a sample, so check this list before naming a sample column.

Being recognized is not the same as being kept: some of these columns are always recomputed by
the build from `trid` and `motif`, and whatever you supply for them is silently replaced. The
"Value you supply" column below says which is which. A kept value wins over TRExplorer's:
TRExplorer fills only the columns your matrix leaves empty.

*Rule 1: recognized by exact name.* The full set, grouped by what it annotates:

| Group | Recognized headers | Value you supply |
|-------|--------------------|------------------|
| gene | `gene_id` / `GencodeGeneId`, `gene_region` / `GencodeGeneRegion` | **kept.** The `Gencode*` alias is used only when the plain name is absent |
| gene | `gene_region_rank` | **discarded.** Always recomputed from `gene_region`, so the rank can never disagree with the region stored beside it (the results table sorts gene regions by this rank while displaying `gene_region`). Supply it only to keep the header out of the sample list; the build prints one line saying how many supplied values it ignored |
| gene | `IsInMendelianGene` | **kept only without a gene table.** With one (genelookup's is downloaded by default) it is recomputed from `gene_id` |
| locus coordinates / motif | `Chrom`, `Start0Based`, `End1Based`, `MotifSize`, `CanonicalMotif` | **discarded.** Always recomputed from `trid` and `motif`. Supply them only to keep those headers out of the sample list |
| locus coordinates / motif | `ReferenceRegion`, `NumRepeatsInReference` | **kept.** Derived from the coordinates and motif only when absent |
| known-locus flags | `KnownDiseaseLocus`, `IsKnownMotif` | **kept only without a known-loci catalog** (`--known-loci-json` / STRchive). With a catalog the catalog wins |
| other annotations | `Source` | **kept.** Falls back to `--source-label` when absent |
| other annotations | `NonCodingAnnotations`, `RepeatMaskerIntervals` | **kept** (pure pass-through) |
| variation cluster | `VariationCluster`, `VariationClusterSizeDiff`, `VariationClusterFilterReason` | **kept** (pure pass-through) |

The recomputed columns are assigned unconditionally in
`locus_annotations.add_derived_locus_columns` (e.g. `record["Chrom"] = f"chr{chrom_field.replace('chr', '')}"`,
`record["MotifSize"] = len(record["Motif"])`, `record["CanonicalMotif"] = compute_canonical_motif(...)`,
`record["gene_region_rank"] = gene_region_rank(record.get("gene_region"))`), while the kept ones
are guarded by `if record.get("...") is None`. Replacing a supplied value is silent except for
`gene_region_rank`, where the build prints a single summary line naming how many loci carried a
supplied rank that was ignored.

*Rule 2: recognized by prefix.* A header whose normalized name starts with `hprc256`,
`aou1027`, `tenk10k` or `trexplorer` (e.g. `HPRC256_99thPercentile`, `AoU1027_Stdev`,
`TenK10K_MaxAllele`, `TRExplorerMotif`) is a per-cohort population-stat / TRExplorer annotation.
If the name matches a loci-table output column it is stored under that column's canonical
spelling; otherwise it is passed through verbatim under your own header. The prefixes are the
specific cohort/version tokens (not the bare cohort abbreviation), so a sample genuinely named
e.g. `AoU_0001` is still read as a sample.

*Population gate.* The `HPRC256_99thPercentile`, `AoU1027_99thPercentile` and
`TenK10K_99thPercentile` columns (TRExplorer's, unless your matrix supplies them) also gate the
precomputed `NumAffectedUnsolved*` counts (the
#Aff / #Fam columns), the `affected_unsolved_outliers_above_unaffected` table and the phenotype
scores: an affected allele counts only when it is larger than the 99th percentile of every
cohort named by `--population-gate-cohorts` that has a value at the locus. The default is
`HPRC256,AoU1027` (the long-read cohorts); the choices are `HPRC256`, `AoU1027` and `TenK10K`,
comma-separated, e.g. `HPRC256,AoU1027,TenK10K` to also gate on the short-read TenK10K cohort.
A locus where none of the chosen cohorts has a value is not gated. The same all-alleles
percentile gates every outlier type. The build records the chosen cohorts in the `metadata`
table, and `trails.py` rebuilds the database when the flag changes.

Both rules live in one place in the code, which is the source of truth if this table ever falls
behind: `ANNOTATION_COLUMN_CANONICAL_BY_NORMALIZED` and `ANNOTATION_COLUMN_PREFIXES` in
`trails/input_tables.py` (applied by `_is_annotation_column`).

In particular, a locus's `gene_id` (TRExplorer's Gencode gene id, or your own column) is what
joins it to the gene table, which populates the gene-symbol / inheritance / pLI /
phenotype-summary columns and enables gene-level phenotype scoring.

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
| `sample_id` | **required** | joins to the sample columns of the repeat-copy-numbers TSV; may not contain `:` or `,` (both are delimiters in the packed `OutlierSampleIds_*` strings, so either one aborts the build with a `ValueError`) |
| `affected_status` | optional | affected-vs-unaffected comparison + affected-unsolved prioritization. `Affected` / `Possibly Affected` (treated as `Affected`) / `Unaffected` / `Unknown`. Absent → unknown |
| `analysis_status` | optional | solved/unsolved filtering + counts. `Solved` / `Unsolved` / `Unknown` / `Probably Solved` / `Partially Solved`. Absent → unknown |
| `family_id` | optional | distinct-family outlier counts |
| `maternal_id`, `paternal_id` | optional | Mendelian-violation QC (requires trios) |
| `sex` | optional | displayed only: shown in the sample tables (in the swim plot's outlier table only for chrX/chrY loci). No filter or calculation reads it; hemizygosity comes from the genotype having a single allele. `Male` / `Female` / `Unknown` |
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
| `hpo_description` | optional | **accepted and ignored.** TRails matches phenotypes on `term_id` alone and never displays a term name; keep the column if it helps you read the file |

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
| `OutlierSampleIds_AllAlleles` | **required** | `allele:sample_id` entries, comma-separated, **sorted by allele size DESC**, e.g. `11x:GSS225379,11x:PMGRC-111-107-2` |
| `OutlierSampleIds_ShortAlleles` | **required** | same format (short alleles) |
| `OutlierSampleIds_HemizygousAlleles` | **required** | same format (hemizygous) |
| `ReferenceRegion` | optional | `chrom:start-end`; computed from `LocusId` if absent |
| `NumRepeatsInReference` | optional | reference allele size in repeat units; computed if absent |

**Optional annotation keys** (absent ⇒ the matching output column is NULL): `CanonicalMotif`,
`Gencode*` (`GencodeGeneId`/`GencodeGeneRegion`/…), reference-quality fields, `TRExplorer*`,
RepeatMasker intervals, and per-cohort population-stat columns (`HPRC256_*`, `AoU1027_*`,
`TenK10K_*`). These come from TRExplorer's catalog, or from annotation columns of your matrix.
`VariantType` is **not** read by the pipeline (the server synthesizes it only for the
ExpansionHunter export) — optional/ignored.

### 5. Generated result database  *(server input; built for you)*

`trails.py` builds a single `*.with_analysis_columns.duckdb` (DuckDB) and points the server at
it (`trails/results_server.py --db <path>`). You don't author it. It contains 6 tables, plus the two
optional Mendelian-QC tables of section 6:

| Table | Contents |
|-------|----------|
| `loci` | one row per locus (the analysis output columns + annotations) |
| `swim_plot` | one row per outlier allele (sample, allele size, affected status, motif, gene, …) |
| `affected_unsolved_outliers_above_unaffected` | one row per (locus, outlier type, sample) that the #Aff count counts: `LocusId`, `outlier_type`, `sample_id`, `family_id` (NULL when the sample has none) and the sample's largest counted `allele_size`. The Group Results By view reads it. Always present, possibly empty |
| `metadata` | build facts: the callset's sample count (`num_samples`) and the `--population-gate-cohorts` value the build used (`population_gate_cohorts`, comma-separated) |
| `per_locus_phenotype_scores` | gene/phenotype similarity per locus *(only if a phenotype TSV was supplied)* |
| `per_outlier_phenotype_scores` | gene/phenotype similarity per outlier sample *(only if a phenotype TSV was supplied)* |

The server also creates a small local annotations database (your notes/tags) on first run,
next to the result database as `*.with_analysis_columns.duckdb.annotations.duckdb`.

A *result* database built by an older, SQLite-based TRails is not read or converted. Its path
ends in `.db`, the path `trails.py` now derives ends in `.duckdb`, so the first run after the
upgrade simply rebuilds from your input TSVs. The old `.db` file is left alone; delete it
yourself once you are happy with the new one.

Notes and tags from an older, SQLite-based TRails are not carried over either. They stay in the
old `*.with_analysis_columns.db.annotations.db`, which TRails leaves alone, and the new
`*.with_analysis_columns.duckdb.annotations.duckdb` starts empty.

### 6. Mendelian-violation QC  *(optional)*

If your sample-metadata TSV includes trio columns (`maternal_id` / `paternal_id`), TRails can
add two more tables (`mendelian_violations`, `mendelian_violations_per_motif`) to the result
database above, built from the per-allele genotypes + family structure, enabling the QC page.
Absent → the QC page is hidden and the server still starts. Everything stays local; nothing
is uploaded.

---

## Class A — public reference data (fetched automatically)

| File | Flag | Source |
|------|------|--------|
| `genes_to_phenotype.txt` | `--genes-to-phenotype` | Human Phenotype Ontology annotation release (github.com/obophenotype/human-phenotype-ontology) |
| `variant_catalog_without_offtargets.GRCh38.json` | `--known-loci-json` | public `str-analysis` repo (github.com/broadinstitute/str-analysis) |
| STRchive disease-locus table (`STRchive-loci.json`) | `--strchive-loci-json` | github.com/dashnowlab/STRchive, downloaded from a built-in default URL and cached locally by the installer (offline-safe); pass the flag to pin your own copy |
| TRExplorer population statistics and locus annotations (`TRExplorer_v2.1.hg38.annotations.duckdb`) | none | TRExplorer repeat catalog v2.1, hg38 (gs://tandem-repeat-catalog). The installer downloads the annotated catalog table (about 900 MB), keeps the columns TRails uses in this file (about 500 MB) and deletes the download. See "Annotations from TRExplorer" above |
| gene-disease table (`genelookup_gene_table.tsv.gz`) | `--gene-table` | genelookup.broadinstitute.org, exported through its query service: gene id and symbol, aliases, pLI, LoF O/E upper bound, HGNC id, inheritance, LLM phenotype summary and sources. OMIM's licensed columns are not included (genelookup does not export them). Not versioned, so it is downloaded once and refreshed only with `TRAILS_FORCE=1 ./install.sh` or `trails/trails_setup.py --force`; pass `--gene-table` to use your own table instead |
