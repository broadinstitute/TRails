**NOTE:** The first version of TRails will be made available before the end of June, 2026.  Click "Watch" to be notified of the release.

# TRails

<b>T</b>andem <b>R</b>epeat <b>A</b>nalysis: <b>i</b>nteractive out<b>l</b>ier <b>s</b>earch

An interactive search interface for prioritizing tandem repeat outlier expansions in rare disease callsets.

TRails loads tandem-repeat genotypes + sample metadata and provides a 
web-browser-based interactive interface for prioritizing candidate pathogenic expansions.
It runs locally on your computer.


<img width="2592" height="3456" alt="ESHG Poster 2026 tandem repeats" src="https://github.com/user-attachments/assets/45b874bc-ff18-47cc-9a65-506c814e8449" />


## Install

`install.sh` requires python3 to be installed already and downloads TRails + reference data into ~/TRails:

```bash
curl -fsSL https://raw.githubusercontent.com/broadinstitute/TRails/main/install.sh | bash
```

This downloads TRails into `~/TRails`, installs the Python deps, and fetches the public
reference data. It is **resumable, and self-updating** — re-run the same line
any time to update to the latest version (an up-to-date install is a no-op; an interrupted
download resumes). Knobs: `TRAILS_INSTALL_DIR` (install location: defaults to ~/TRails), 
`TRAILS_FORCE=1` (force re-download), `TRAILS_REF` (branch or tag).

From an existing checkout, `./install.sh` does the same minus the self-download.

## Quickstart

1. **Prepare your metadata TSVs.** Put your sample metadata and phenotypes into the
   sample-table and phenotype-table formats in
   [docs/INPUT_FORMATS.md](docs/INPUT_FORMATS.md). (Genotypes need no pre-conversion —
   `prepare_data.py` ingests the raw TRGT formats directly.)

2. **Build the result database(s)** from your raw genotypes in one command. TRails accepts
   a **merged-TRGT-VCF** (`--vcf`) and/or a **TRGT-LPS** table (`--lps`) — supply one or
   both. `prepare_data.py` converts each raw input to the intermediate JSON and runs the
   analysis to produce the DB:
   ```bash
   python3 prepare_data.py \
       --lps your_cohort.lps.txt.gz \
       --vcf your_cohort.trgt.vcf.gz \
       --sample-table your_samples.tsv.gz \
       --phenotypes-table your_phenotypes.tsv.gz \
       --genes-to-phenotype reference_data/genes_to_phenotype.txt
   # Drop --vcf (or --lps) to build just one DB. The gene-disease table is optional.
   # The known-loci catalog defaults to the str-analysis GitHub JSON (override with --known-loci-json).
   ```
   This writes a `*.with_analysis_columns.db` (and `.tsv.gz`) per genotype input.

3. **Start the server:**
   ```bash
   ./start.sh        # or: python3 results_server.py --lps-db <db> [--vcf-db <db>]
   ```
   Pass whichever DB(s) you built — `--lps-db`, `--vcf-db`, or both.
   Open the printed URL. Browse loci, the swim plot, filter by motif size / affected
   status / gene, and export candidates (BED / TSV / JSON / ExpansionHunter).

## License

MIT — see [LICENSE](LICENSE).
