"""Build the TRails SQLite result database directly from the user's input TSVs.

This is TRails' analysis engine and pipeline orchestrator. It reads the
repeat-copy-numbers matrix TSV + the sample-metadata TSV (plus the optional
phenotype / gene / reference inputs) and — entirely in memory — computes the
allele histograms, detects outliers, derives every analysis column, scores
phenotypes, builds the swim-plot / skinny / phenotype / Mendelian tables, and
writes them all straight into ``db_path``. Nothing is written to disk except the
database itself: there is no intermediate histogram or JSON file.

It ports the staged logic of the internal ``analyze_results.py`` pipeline (which
read a pre-built per-locus JSON); here the same per-locus records are built
directly from the parsed matrix instead of being round-tripped through a file.

``build()`` wires the modules together in this order (each module is its own
fully-tested unit):

  1. ``input_tables.read_repeat_copy_numbers``  -> locus rows + sample-id list
  2. ``input_tables.read_sample_metadata``       -> sample / affected / analysis lookups
  3. ``input_tables.read_phenotypes``            -> participant -> HPO terms (optional)
  4. ``input_tables.read_gene_table`` + ``read_gene_disease_phenotypes`` (optional)
  5. ``allele_histograms.build_histograms_and_outliers`` -> per-locus records
  6. ``locus_annotations.add_derived_locus_columns``     -> coordinate / motif columns
  7. known-disease-locus match + ``CanonicalMotif`` + ``IsKnownMotif`` + ``IsInMendelianGene``
  8. ``analysis_columns.add_gene_columns``               -> pLI / inheritance
  9. ``analysis_columns.add_all_outlier_columns``        -> the affected/unaffected analysis
 10. ``phenotype_scoring.compute_phenotype_scores``      -> phenotype-score columns + rows
 11. ``result_database``: loci table + indexes + skinny tables + swim plot +
     phenotype tables + (when complete trios exist) Mendelian tables, finalized
     with an atomic move.

CLI:
    python3 build_database.py --repeat-copy-numbers-tsv X --sample-metadata-tsv Y --db out.db
"""

import argparse

import allele_histograms
import analysis_columns
import input_tables
import locus_annotations
import mendelian_qc
import phenotype_scoring
import result_database
import swim_plot


def build(repeat_copy_numbers_tsv, sample_metadata_tsv, db_path,
          phenotypes_table=None, gene_table=None, genes_to_phenotype=None,
          known_loci_json=None, n_loci=None, skip_phenotype_scores=False,
          source_label="", fetch_strchive=False, mendelian_threshold=2,
          n_outlier_sample_ids=10, strchive_loci_json=None):
    """Read the input TSVs and write the TRails result database to db_path.

    Builds every table in memory and populates the SQLite database directly — no
    intermediate files. Optional inputs (phenotypes, gene table, known-loci
    catalog) are simply skipped when absent, and the build still runs: a missing
    gene/phenotype input leaves the corresponding columns NULL rather than
    hard-failing.

    Args:
        repeat_copy_numbers_tsv: Merged per-allele genotype matrix (trid, motif,
            one column per sample).
        sample_metadata_tsv: One row per sample (only sample_id required).
        db_path: Output SQLite database path. Written atomically (``<path>.tmp``
            then ``os.replace``).
        phenotypes_table: Optional HPO-terms-per-sample TSV.
        gene_table: Optional gene-disease table (gene_id keyed).
        genes_to_phenotype: Optional HPO genes_to_phenotype.txt (reference) for
            phenotype scoring.
        known_loci_json: Optional known-disease-locus catalog JSON (class-A
            ``variant_catalog_without_offtargets.GRCh38.json``). Read from the
            cached file; no network fetch.
        n_loci: Optional cap on the number of loci (for testing / quick builds).
        skip_phenotype_scores: Skip phenotype scoring when True.
        source_label: Value stored in the loci ``Source`` column (and swim-plot
            ``SourceDb``). Defaults to the empty string.
        fetch_strchive: When True, allow ``load_known_disease_loci`` to fetch the
            STRchive locus set over the network. Off by default for hermetic builds.
        mendelian_threshold: Alleles match within strictly fewer than this many
            repeats when checking Mendelian violations. Defaults to 2.
        n_outlier_sample_ids: Outlier-sample cap passed to the histogram builder.
            Defaults to 10.
        strchive_loci_json: Optional path to a cached STRchive-loci.json file used
            as a known-disease-locus fallback (no network access).

    Returns:
        The path of the written database (``db_path``).
    """
    # Stage 1: read the repeat-copy-numbers matrix into per-locus rows.
    print(f"Reading repeat-copy-numbers matrix: {repeat_copy_numbers_tsv}")
    locus_rows, sample_id_list = input_tables.read_repeat_copy_numbers(repeat_copy_numbers_tsv)
    if n_loci is not None:
        locus_rows = locus_rows[:n_loci]
    print(f"  {len(locus_rows):,} loci, {len(sample_id_list):,} samples")

    # Stage 2: read the sample metadata and build the lookups.
    print(f"Reading sample metadata: {sample_metadata_tsv}")
    sample_lookup, affected_lookup, analysis_lookup, sample_df = (
        input_tables.read_sample_metadata(sample_metadata_tsv))
    print(f"  {len(sample_lookup):,} samples with metadata")

    # Warn about matrix sample columns that have no metadata row (treated as
    # Unknown downstream; never a hard-fail).
    missing_metadata = [s for s in sample_id_list if s not in sample_lookup]
    if missing_metadata:
        print(f"  WARNING: {len(missing_metadata):,} matrix sample column(s) have no "
              f"metadata row and will be treated as Unknown "
              f"(e.g. {missing_metadata[:5]})")

    # Stage 3: read the optional phenotypes table.
    participant_to_hpo = input_tables.read_phenotypes(phenotypes_table)
    if participant_to_hpo:
        print(f"  {len(participant_to_hpo):,} samples with HPO phenotype terms")

    # Stage 4: read the optional gene table + gene-disease phenotype table.
    gene_lookup = input_tables.read_gene_table(gene_table)
    if gene_lookup:
        print(f"  {len(gene_lookup):,} genes in the gene table")
    gene_disease_data = input_tables.read_gene_disease_phenotypes(genes_to_phenotype)
    if gene_disease_data:
        print(f"  {len(gene_disease_data):,} genes in the gene-disease phenotype table")

    # Stage 5: build the histograms and outlier-sample lists.
    print("Building allele histograms and outlier-sample lists ...")
    records = allele_histograms.build_histograms_and_outliers(
        locus_rows, sample_id_list, n_outlier_sample_ids=n_outlier_sample_ids)

    # The histogram builder keeps the raw 'trid' / 'motif' keys; promote them to
    # the canonical LocusId / Motif columns the annotation stages expect, and
    # promote any per-locus annotation columns carried in 'extra_columns'
    # (gene_id, GencodeGeneRegion, HPRC256_* ...) onto the record top level so the
    # gene / population stages can read them.
    for record in records:
        record["LocusId"] = record["trid"]
        record["Motif"] = record["motif"]
        record.update(record.pop("extra_columns", {}))

    # Stage 6: derived coordinate / motif / source / gene-region columns.
    print("Adding derived locus columns ...")
    locus_annotations.add_derived_locus_columns(records, source_label=source_label)

    # Stage 7: known-disease-locus match + IsKnownMotif + IsInMendelianGene.
    print("Annotating known-disease loci and known motifs ...")
    if known_loci_json or strchive_loci_json or fetch_strchive:
        interval_trees, strchive_trees, _locus_lookup = locus_annotations.load_known_disease_loci(
            known_loci_json, fetch_strchive=fetch_strchive,
            strchive_filepath=strchive_loci_json)
        known_canonical_motifs = locus_annotations.collect_known_disease_canonical_motifs(
            interval_trees, strchive_trees)
    else:
        interval_trees, strchive_trees = {}, {}
        known_canonical_motifs = set()
        print("  no known-loci catalog supplied; KnownDiseaseLocus / IsKnownMotif left NULL/0")

    for record in records:
        record["KnownDiseaseLocus"] = locus_annotations.matches_disease_locus(
            record["LocusId"], interval_trees, strchive_trees) if (interval_trees or strchive_trees) else None
        record["IsKnownMotif"] = locus_annotations.compute_is_known_motif(
            record["CanonicalMotif"], known_canonical_motifs)
        record["IsInMendelianGene"] = locus_annotations.is_in_mendelian_gene(
            record.get("gene_id"), gene_lookup)

    # Stage 8: gene-derived pLI / inheritance columns.
    print("Adding gene columns ...")
    analysis_columns.add_gene_columns(records, gene_lookup)

    # Stage 9: the affected/unaffected outlier analysis (the heart).
    print("Computing outlier-analysis columns ...")
    analysis_columns.add_all_outlier_columns(
        records, sample_lookup, affected_lookup, analysis_lookup)

    # Stage 10: phenotype scoring (optional).
    per_outlier_phenotype_rows, per_locus_phenotype_rows = [], []
    if skip_phenotype_scores:
        print("Skipping phenotype scoring (--skip-phenotype-scores).")
    elif not participant_to_hpo:
        print("Skipping phenotype scoring (no phenotypes table supplied).")
    else:
        print("Scoring phenotypes ...")
        per_outlier_phenotype_rows, per_locus_phenotype_rows = (
            phenotype_scoring.compute_phenotype_scores(
                records, participant_to_hpo, gene_lookup, gene_disease_data,
                affected_lookup, analysis_lookup))
        print(f"  {len(per_outlier_phenotype_rows):,} per-outlier + "
              f"{len(per_locus_phenotype_rows):,} per-locus phenotype-score rows")

    # Stage 11: write everything into the database (atomic move at the end).
    print(f"Writing database: {db_path}")
    connection, tmp_path = result_database.open_new_database(db_path)
    try:
        _, present_columns = result_database.write_loci_table(
            connection, records, analysis_columns.OUTPUT_COLUMNS)
        result_database.create_loci_indexes(connection, present_columns)
        result_database.write_skinny_tables(connection, present_columns)

        swim_rows = swim_plot.generate_swim_plot_table(
            records, sample_lookup, affected_lookup, analysis_lookup)
        result_database.write_swim_plot(connection, swim_rows)

        result_database.write_phenotype_tables(
            connection, per_outlier_phenotype_rows, per_locus_phenotype_rows)

        per_sample_mendelian_rows, per_motif_mendelian_rows = (
            mendelian_qc.compute_mendelian_violations(
                locus_rows, sample_lookup, sample_df, threshold=mendelian_threshold))
        result_database.write_mendelian_tables(
            connection, per_sample_mendelian_rows, per_motif_mendelian_rows)

        result_database.finalize_database(connection, tmp_path, db_path)
    except Exception:
        connection.close()
        raise

    print(f"Done: {db_path}")
    return db_path


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repeat-copy-numbers-tsv", required=True,
                        help="merged repeat-copy-numbers TSV (one row per locus, one column per sample)")
    parser.add_argument("--sample-metadata-tsv", required=True,
                        help="sample metadata TSV (one row per sample; only sample_id is required)")
    parser.add_argument("--db", required=True, help="output SQLite database path")
    parser.add_argument("--phenotypes-table", help="optional phenotype TSV (HPO terms per sample)")
    parser.add_argument("--gene-table", help="optional gene-disease table (class B)")
    parser.add_argument("--genes-to-phenotype", help="HPO genes_to_phenotype.txt (reference)")
    parser.add_argument("--known-loci-json", help="known-loci catalog (reference)")
    parser.add_argument("--strchive-loci-json", help="cached STRchive-loci.json (reference, optional)")
    parser.add_argument("--source-label", default="", help="value for the loci 'Source' column")
    parser.add_argument("--fetch-strchive", action="store_true",
                        help="allow a live STRchive network fetch (off by default)")
    parser.add_argument("--mendelian-threshold", type=int, default=2,
                        help="alleles match within strictly fewer than this many repeats (default 2)")
    parser.add_argument("--n-outlier-sample-ids", type=int, default=10,
                        help="outlier-sample cap for the histogram builder (default 10)")
    parser.add_argument("-n", "--n-loci", type=int, help="limit to the first N loci (for testing)")
    parser.add_argument("--skip-phenotype-scores", action="store_true", help="skip phenotype scoring")
    args = parser.parse_args()
    build(args.repeat_copy_numbers_tsv, args.sample_metadata_tsv, args.db,
          phenotypes_table=args.phenotypes_table, gene_table=args.gene_table,
          genes_to_phenotype=args.genes_to_phenotype, known_loci_json=args.known_loci_json,
          n_loci=args.n_loci, skip_phenotype_scores=args.skip_phenotype_scores,
          source_label=args.source_label, fetch_strchive=args.fetch_strchive,
          mendelian_threshold=args.mendelian_threshold,
          n_outlier_sample_ids=args.n_outlier_sample_ids,
          strchive_loci_json=args.strchive_loci_json)


if __name__ == "__main__":
    main()
