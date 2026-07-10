# TRails

<b>T</b>andem <b>R</b>epeat <b>A</b>nalysis: <b>i</b>nteractive out<b>l</b>ier <b>s</b>earch

An interactive search interface for prioritizing tandem repeat outlier expansions in rare disease callsets.

TRails loads tandem-repeat genotypes + sample metadata and provides an interactive interface for prioritizing candidate pathogenic expansions.

## Features

- **Runs on your computer** TRails is a local app - your genotypes and
  sample metadata are never sent to any remote server or third party. 
- **Supports short-read or long-read TR genotypes.** TRails only needs repeat copy-numbers and 
  sample metadata in TSV format.
- **Interactive outlier search** Filter by population allele frequencies, gene regions, motif categories, etc. Sort by allele size, outlier count, and
  more.
- **Tagging** Define custom tags to mark specific search results as Candidates, Reviewed, etc.  
- **Phenotype-aware prioritization** — rank loci based on the similarity of HPO terms among outlier sample(s) with each other as well as with the gene's known disease phenotype.
- **Sample QC** flags samples that have an excess of large expansions, or of inheritance-inconsistent genotypes.
- **Export** search results in BED, TSV, JSON, or ExpansionHunter catalog formats.

<img width="2592" height="3456" alt="ESHG Poster 2026 tandem repeats" src="https://github.com/user-attachments/assets/45b874bc-ff18-47cc-9a65-506c814e8449" />

## Install

`install.sh` requires python3 and downloads TRails + public reference data into `./TRails`:

```bash
curl -fsSL https://raw.githubusercontent.com/broadinstitute/TRails/main/install.sh | bash
```

It installs the Python dependencies and fetches the public reference data. It is
**resumable and self-updating** — re-run the same line any time to update (an up-to-date
install is a no-op; an interrupted download resumes). Optional environment variables: `TRAILS_INSTALL_DIR` (install
location; default `./TRails`), `TRAILS_FORCE=1` (force re-download).

From an existing checkout, `./install.sh` does the same minus the self-download. You can also
skip `install.sh` entirely — `trails.py` installs anything missing on first run.

## Quickstart

One command takes your two tables and starts the server:

```bash
python3 trails.py \
    --repeat-copy-numbers-tsv  your_cohort.tsv.gz \
    --sample-metadata-tsv      your_samples.tsv.gz \
    [--phenotypes-table        your_phenotypes.tsv.gz]
```

It: (1) installs any missing dependencies, (2) downloads any missing
public reference data, (3) builds the local SQLite database (or reuses it if your inputs are
unchanged — pass `--rebuild` to force), and (4) starts the local server. Open the printed
URL in your web browser to begin your analysis using TRails interface. 

Useful `trails.py` flags: 
`--rebuild` (force a rebuild), `--port` / `--host`, `-n / --n-loci N` (limit to the first N
loci, for a quick test).

## Inputs

TRails reads **two** tables you supply (a third is optional). **Only a handful of columns are
required — almost everything is optional, and any missing optional column simply disables the
feature that depends on it rather than causing an error.** Full schemas and the generated
database tables are documented in [docs/INPUT_FORMATS.md](docs/INPUT_FORMATS.md).

### 1. Repeat copy numbers TSV — `--repeat-copy-numbers-tsv`  *(required)*

Tab-separated (optionally gzipped). **One row per locus**, the first two columns identify the
locus and the remaining columns are **one per sample**. Each cell holds that sample's
per-allele repeat copy numbers, comma-separated (two values for a diploid call, one for a
hemizygous/haploid call); leave a cell blank or `.` for a no-call.

```
trid                            motif   SAMPLE_A   SAMPLE_B   SAMPLE_C
chr1-57367043-57367119-AAAAT    AAAAT   12,12      11,13      12,40
chrX-148500631-148500691-GCC    GCC     21         20,21      45
```

| Column | Required? | Notes |
|--------|-----------|-------|
| `trid` | **required** | locus id, `chrom-start-end-motif` (the TRGT repeat id) |
| `motif` | **required** | repeat unit, e.g. `AAAAT` |
| one column per sample | **≥1 required** | header = `sample_id` (must match the sample-metadata TSV); cell = comma-separated per-allele copy numbers |

> This is a per-allele genotype matrix from any TR genotyping tool. `trails.py` reads it
> directly into the local database — no intermediate files are written (see
> [docs/INPUT_FORMATS.md](docs/INPUT_FORMATS.md)).

### 2. Sample metadata TSV — `--sample-metadata-tsv`  *(required)*

Tab-separated (optionally gzipped), **one row per sample**. Only `sample_id` is required;
every other column is optional and unlocks additional functionality when present. Column names are not case-sensitive  
and the _ is optional, so `sample_id`, `Sample_id` and `SampleId` are equivalent.

```
sample_id      affected_status   analysis_status   family_id   maternal_id    paternal_id    sex      phenotype_description
PMGRC-1-1-0    Affected          Unsolved          FAM1        PMGRC-1-2-1    PMGRC-1-3-2    Male     Muscle weakness; Myopathy
PMGRC-1-2-1    Unaffected        Unaffected        FAM1                                      Female   Unaffected mother
```

| Column | Required? | Enables / notes |
|--------|-----------|-----------------|
| `sample_id` | **required** | joins to the sample columns of the repeat-copy-numbers TSV |
| `affected_status` | optional | affected-vs-unaffected outlier comparison + affected-unsolved prioritization. Values: `Affected`, `Possibly Affected` (treated as `Affected`), `Unaffected`, `Unknown`. Absent → treated as unknown |
| `analysis_status` | optional | solved/unsolved filtering + counts. Values: `Solved`, `Unsolved`, `Unknown`, … Absent → unknown |
| `family_id` | optional | distinct-family outlier counts; pairs with the parent ids for trio QC |
| `maternal_id`, `paternal_id` | optional | Mendelian-violation QC (needs trios) |
| `sex` | optional | sex-aware handling of hemizygous loci. Values: `Male`, `Female`, `Unknown` |
| `phenotype_description` | optional | free-text shown in the UI |

Any extra columns are ignored.

### 3. Phenotype TSV — `--phenotypes-table`  *(optional)*

Tab-separated (optionally gzipped), **one row per (sample, HPO term)**. Supplying it enables
phenotype-aware prioritization; omit it and TRails still runs.

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

### Optional licensed reference inputs  *(not bundled)*

This improves annotation but is license-restricted, so TRails never downloads it; pass
your own copy only if you are licensed. Absent → the dependent annotations (gene symbol,
inheritance, pLI, disease-category and phenotype-summary columns, and gene-level phenotype
scoring) are simply omitted.

| Input | Flag | Restriction |
|-------|------|-------------|
| gene–disease + phenotype-summary table | `--gene-table` | OMIM/HPO-derived |

## Dependencies

- **python3** (3.8+) and **pip**. The one-line installer additionally needs **curl** and
  **tar** (standard on macOS/Linux). No other system packages are required.
- **Python packages** (installed automatically from `requirements.txt` by `install.sh` /
  `trails.py`):
  - `flask` — the local web server
  - `pandas`, `numpy` — the analysis pipeline
  - `tqdm` — progress bars
  - `msgpack` — the server's startup cache
  - `intervaltree` — interval overlap lookups
  - `pyhpo` — HPO term similarity (phenotype scoring)
  - `requests` — fetching public reference data (e.g. STRchive)

  TRails is otherwise self-contained — it has no other third-party runtime dependency (the
  motif and locus-id utilities are built in, not pulled from an external package).

You normally never install these by hand — `install.sh` and `trails.py` do it for you. To
install manually: `python3 -m pip install -r requirements.txt`.

## License

MIT — see [LICENSE](LICENSE).
