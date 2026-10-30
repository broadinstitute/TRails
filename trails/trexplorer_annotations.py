"""Population statistics and locus annotations from TRExplorer's repeat catalog.

trails_setup downloads TRExplorer's annotated hg38 catalog table once and keeps only the columns
TRails uses, in a small DuckDB database (build_trexplorer_annotations_database). The build then
gives each locus of the repeat-copy-numbers table the annotations of its TRExplorer locus
(add_trexplorer_annotations), chosen as follows:

  * a TRExplorer locus is a candidate when it is on the same chromosome, its interval overlaps
    the locus with a Jaccard similarity (intersection over union) of at least
    MIN_JACCARD_SIMILARITY, and its motif matches: the same canonical motif when both motifs
    are 1-6 bp, the same motif length when the locus's motif is 7 bp or more;
  * the candidate with the locus's own id wins (chr prefixes ignored);
  * otherwise each candidate is ranked by the first of these it has: its HPRC256 standard
    deviation, its AoU1027 standard deviation, its reference repeat purity. A candidate ranked
    by an earlier one of the three beats a candidate ranked by a later one or by none, and
    between two candidates ranked by the same one, the larger value wins;
  * a tie that is left goes to the candidate with the last locus id in sort order.
"""

import pandas

import duckdb_compat
import locus_annotations
import result_database

MIN_JACCARD_SIMILARITY = 0.66

TREXPLORER_LOCI_TABLE = "trexplorer_loci"

# Columns copied from the matched TRExplorer locus onto the TRails locus: {TRExplorer column:
# TRails column}. A value the repeat-copy-numbers table already supplies for a column is kept.
TRAILS_COLUMN_BY_TREXPLORER_COLUMN = {
    "LocusId": "TRExplorerLocusId",
    "Motif": "TRExplorerMotif",
    "ReferenceRegion": "TRExplorerReferenceRegion",
    "ReferenceRepeatPurity": "TRExplorerReferenceRepeatPurity",
    "Source": "TRExplorerSource",
    "GencodeGeneId": "gene_id",
    "GencodeGeneRegion": "gene_region",
    "VariationCluster": "VariationCluster",
    "VariationClusterSizeDiff": "VariationClusterSizeDiff",
    "VariationClusterFilterReason": "VariationClusterFilterReason",
}
INTEGER_POPULATION_COLUMNS = [
    "HPRC256_MaxAllele", "HPRC256_StdevRankByMotif", "HPRC256_StdevRankTotalNumberByMotif",
    "AoU1027_MaxAllele", "AoU1027_StdevRankByMotif", "AoU1027_StdevRankTotalNumberByMotif",
    "TenK10K_MaxAllele",
]
FLOAT_POPULATION_COLUMNS = [
    "HPRC256_99thPercentile", "HPRC256_Stdev",
    "AoU1027_99thPercentile", "AoU1027_Stdev", "AoU1027_OE_Length", "AoU1027_OE_LengthPercentile",
    "TenK10K_99thPercentile", "TenK10K_Stdev",
]

# A chromosome's place in one genome-wide coordinate (CHROMOSOME_INDEX * GENOME_COORDINATE_STRIDE
# + position), which lets the overlap join compare intervals across all chromosomes with nothing
# but two inequalities. DuckDB answers that with its sort-based IEJoin; joining on the chromosome
# name instead makes it hash-join per chromosome and compare every pair of loci within one.
GENOME_COORDINATE_STRIDE = 10_000_000_000
CHROMOSOME_INDEX_SQL = ("CASE {chrom} WHEN 'X' THEN 23 WHEN 'Y' THEN 24 WHEN 'M' THEN 25 WHEN 'MT' THEN 25 "
                        "ELSE TRY_CAST({chrom} AS INTEGER) END")


def _sql_string_literal(text):
    """Return text as a single-quoted SQL string literal."""
    return "'" + text.replace("'", "''") + "'"


def _without_chr_prefix(text):
    """Return a chromosome or locus id without a leading "chr", as TRExplorer spells them."""
    return text[3:] if text.startswith("chr") else text


def build_trexplorer_annotations_database(annotations_tsv_path, database_path):
    """Keep the columns TRails uses from TRExplorer's annotated catalog table, in a DuckDB file.

    The database is built at a temporary path and moved into place with the result database's
    own helpers (result_database.open_new_database and finalize_database), so an interrupted
    build never leaves a partial database that a later run would take for a finished one, and
    a write-ahead log left next to the old copy is not replayed into the new one.

    Args:
        annotations_tsv_path: TRExplorer's ``*.EH.with_annotations.tsv.gz`` catalog table.
        database_path: Where to write the database.

    Returns:
        The number of TRExplorer loci written.
    """
    region_coordinates = "split_part(ReferenceRegion, ':', 2)"
    numeric_columns = (
        [f"TRY_CAST({c} AS BIGINT) AS {c}" for c in INTEGER_POPULATION_COLUMNS + ["VariationClusterSizeDiff"]]
        + [f"TRY_CAST({c} AS DOUBLE) AS {c}" for c in FLOAT_POPULATION_COLUMNS + ["ReferenceRepeatPurity"]])
    text_columns = ["Motif", "CanonicalMotif", "ReferenceRegion", "Source", "GencodeGeneId", "GencodeGeneRegion",
                    "VariationCluster", "VariationClusterFilterReason"]
    connection, tmp_path = result_database.open_new_database(database_path)
    try:
        connection.execute(
            f"CREATE TABLE {TREXPLORER_LOCI_TABLE} AS SELECT "
            f"regexp_replace(LocusId, '^chr', '') AS LocusId, "
            f"regexp_replace(split_part(ReferenceRegion, ':', 1), '^chr', '') AS Chrom, "
            f"CAST(split_part({region_coordinates}, '-', 1) AS BIGINT) AS Start0Based, "
            f"CAST(split_part({region_coordinates}, '-', 2) AS BIGINT) AS End1Based, "
            f"{', '.join(text_columns + numeric_columns)} "
            f"FROM read_csv({_sql_string_literal(annotations_tsv_path)}, delim='\t', header=true, "
            f"all_varchar=true, quote='', escape='')")
        num_loci = connection.execute(f"SELECT COUNT(*) FROM {TREXPLORER_LOCI_TABLE}").fetchone()[0]
        result_database.finalize_database(connection, tmp_path, database_path)
    except Exception:
        connection.close()
        raise
    return num_loci


def _match_loci(locus_records, annotations_database_path):
    """Find each locus's TRExplorer locus.

    Args:
        locus_records: Per-locus dicts carrying LocusId, Chrom, Start0Based, End1Based, Motif and
            CanonicalMotif (see locus_annotations.add_derived_locus_columns).
        annotations_database_path: The database build_trexplorer_annotations_database wrote.

    Returns:
        A dict mapping the index of each matched record in locus_records to a dict of its
        TRExplorer locus's columns plus ``jaccard_similarity``.
    """
    user_loci = pandas.DataFrame({
        "record_index": range(len(locus_records)),
        "LocusId": [_without_chr_prefix(r["LocusId"]) for r in locus_records],
        "Chrom": [_without_chr_prefix(r["Chrom"]) for r in locus_records],
        "Start0Based": [r["Start0Based"] for r in locus_records],
        "End1Based": [r["End1Based"] for r in locus_records],
        "MotifSize": [len(r["Motif"]) for r in locus_records],
        "CanonicalMotif": [r["CanonicalMotif"] for r in locus_records],
    })
    overlap = "(LEAST(u.End1Based, t.End1Based) - GREATEST(u.Start0Based, t.Start0Based))"
    union = f"((u.End1Based - u.Start0Based) + (t.End1Based - t.Start0Based) - {overlap})"
    # The tie-break of the module docstring: the tier says which value a candidate is ranked by
    # (HPRC256 stdev, then AoU1027 stdev, then reference purity), then that value, then the locus id.
    tier = ("CASE WHEN t.HPRC256_Stdev IS NOT NULL THEN 2 WHEN t.AoU1027_Stdev IS NOT NULL THEN 1 "
            "WHEN t.ReferenceRepeatPurity IS NOT NULL THEN 0 ELSE -1 END")
    tier_value = "COALESCE(t.HPRC256_Stdev, t.AoU1027_Stdev, t.ReferenceRepeatPurity, 0)"

    def genome_coordinate(alias, column):
        return (f"{CHROMOSOME_INDEX_SQL.format(chrom=f'{alias}.Chrom')} * {GENOME_COORDINATE_STRIDE} "
                f"+ {alias}.{column}")

    connection = duckdb_compat.connect(":memory:")
    try:
        connection.execute(f"ATTACH {_sql_string_literal(annotations_database_path)} AS trexplorer (READ_ONLY)")
        connection.register("user_loci", user_loci)
        cursor = connection.execute(
            f"WITH candidates AS ("
            f"  SELECT u.record_index, u.LocusId = t.LocusId AS same_locus_id, "
            f"         {overlap}::DOUBLE / NULLIF({union}, 0) AS jaccard_similarity, t.* "
            f"  FROM user_loci u JOIN trexplorer.{TREXPLORER_LOCI_TABLE} t "
            f"    ON {genome_coordinate('t', 'Start0Based')} < {genome_coordinate('u', 'End1Based')} "
            f"   AND {genome_coordinate('t', 'End1Based')} > {genome_coordinate('u', 'Start0Based')} "
            f"  WHERE (u.MotifSize <= 6 AND length(t.Motif) <= 6 AND u.CanonicalMotif = t.CanonicalMotif) "
            f"     OR (u.MotifSize >= 7 AND length(t.Motif) = u.MotifSize)"
            f") "
            f"SELECT * EXCLUDE (same_locus_id) FROM candidates t "
            f"WHERE jaccard_similarity >= {MIN_JACCARD_SIMILARITY} "
            f"QUALIFY ROW_NUMBER() OVER (PARTITION BY record_index ORDER BY same_locus_id DESC, "
            f"  {tier} DESC, {tier_value} DESC, LocusId DESC) = 1")
        column_names = [description[0] for description in cursor.description]
        matches = {row[0]: dict(zip(column_names, row)) for row in cursor.fetchall()}
    finally:
        connection.close()
    return matches


def add_trexplorer_annotations(locus_records, annotations_database_path):
    """Copy each locus's TRExplorer population statistics and annotations onto it.

    A value the repeat-copy-numbers table already supplies for a column is kept; TRExplorer
    fills only the columns it left empty. A locus whose gene region is filled here gets its
    gene_region_rank recomputed, since add_derived_locus_columns derived it before.

    Args:
        locus_records: Per-locus dicts after locus_annotations.add_derived_locus_columns.
            Mutated in place.
        annotations_database_path: The database build_trexplorer_annotations_database wrote.

    Returns:
        A tuple (number of loci matched by locus id, number matched by overlap).
    """
    matches = _match_loci(locus_records, annotations_database_path)
    num_matched_by_locus_id = 0
    for record_index, match in matches.items():
        record = locus_records[record_index]
        if match["LocusId"] == _without_chr_prefix(record["LocusId"]):
            num_matched_by_locus_id += 1
        if record.get("TRExplorerLocusJaccardSimilarity") is None:
            record["TRExplorerLocusJaccardSimilarity"] = round(match["jaccard_similarity"], 4)
        for trexplorer_column, trails_column in TRAILS_COLUMN_BY_TREXPLORER_COLUMN.items():
            if record.get(trails_column) is None:
                record[trails_column] = match[trexplorer_column]
        for column in INTEGER_POPULATION_COLUMNS + FLOAT_POPULATION_COLUMNS:
            if record.get(column) is None:
                record[column] = match[column]
        record["gene_region_rank"] = locus_annotations.gene_region_rank(record.get("gene_region"))
    return num_matched_by_locus_id, len(matches) - num_matched_by_locus_id
