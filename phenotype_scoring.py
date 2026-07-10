"""HPO gene-phenotype and pairwise phenotype scoring for TRails.

This module ports the *inlined* phenotype-scoring path of the original
``analyze_results.py`` pipeline (its ``compute_phenotype_scores`` plus the
low-level helpers it borrowed from ``gene_phenotype_scorer.py``). It does NOT
port the standalone ``process_all_loci`` CLI.

For each locus and outlier type (AllAlleles / ShortAlleles / HemizygousAlleles)
it selects the *qualifying* affected/unsolved outlier samples, then computes two
families of scores:

  * A gene-phenotype similarity: how well a patient's HPO terms match the known
    disease phenotypes of the gene at that locus, restricted to diseases whose
    inheritance mode is compatible with the outlier type.
  * Pairwise phenotype similarity and shared-term counts between consecutive
    qualifying samples (ordered largest allele first), which surface loci where
    multiple expanded patients share phenotypes.

The semantic similarity uses the ``pyhpo`` ontology when available; ``pyhpo`` is
a SOFT dependency, so on ``ImportError`` everything degrades to a Jaccard
overlap (and the information-content weighting collapses to a flat weight of
1.0 per shared term). The build still runs end to end without ``pyhpo``.

Internal representation follows the rest of TRails: a ``records`` list with one
dict per locus, and ``participant_to_hpo`` mapping each sample_id to a set of
HPO term-id strings (as produced by ``input_tables.read_phenotypes``).
"""

from analysis_columns import (
    is_affected_unsolved,
    is_above_population_p99,
    is_above_unaffected,
    is_missing_outlier_value,
    parse_outlier_entries,
)


# Which disease inheritance modes are compatible with each outlier type:
#   - AllAlleles (expansions on one allele) -> dominant diseases (AD, XD)
#   - ShortAlleles (the shorter of two alleles) -> recessive diseases (AR)
#   - HemizygousAlleles (single allele, X in males) -> XR, XL, XD all manifest
OUTLIER_TYPE_INHERITANCE = {
    "AllAlleles": {"AD", "XD"},
    "ShortAlleles": {"AR"},
    "HemizygousAlleles": {"XR", "XL", "XD"},
}

OUTLIER_TYPES = ["AllAlleles", "ShortAlleles", "HemizygousAlleles"]


# Whether the pyhpo Ontology() singleton has already been initialized. Building
# the ontology is expensive, so it is done at most once per process.
_ontology_initialized = False


def ensure_ontology():
    """Initialize the pyhpo Ontology singleton once, if pyhpo is installed.

    Safe to call repeatedly: the first successful call builds the ontology and
    subsequent calls are no-ops. If pyhpo is not installed the call returns
    quietly (the Jaccard fallback needs no ontology).

    Returns:
        bool: True if the ontology is available (now or already), False if pyhpo
        could not be imported.
    """
    global _ontology_initialized
    if _ontology_initialized:
        return True
    try:
        from pyhpo import Ontology
    except ImportError:
        return False
    Ontology()
    _ontology_initialized = True
    return True


def compute_similarity_pyhpo(patient_hpo_ids, disease_hpo_ids):
    """Compute semantic similarity between two HPO term sets.

    Uses pyhpo's graph-information-content similarity (kind='omim',
    method='graphic', combine='funSimAvg') when pyhpo is installed; otherwise
    falls back to the Jaccard overlap of the two term sets. Returns 0.0 when
    either set is empty (or, with pyhpo, when neither set has any HP: terms).

    Args:
        patient_hpo_ids: Iterable of the patient's HPO term ids.
        disease_hpo_ids: Iterable of a disease's HPO term ids.

    Returns:
        float: A similarity score in [0.0, 1.0].
    """
    try:
        from pyhpo import HPOSet

        if not ensure_ontology():
            raise ImportError("pyhpo ontology unavailable")

        patient_valid = [term for term in patient_hpo_ids if term.startswith("HP:")]
        disease_valid = [term for term in disease_hpo_ids if term.startswith("HP:")]
        if not patient_valid or not disease_valid:
            return 0.0

        try:
            return HPOSet.from_queries(patient_valid).similarity(
                HPOSet.from_queries(disease_valid),
                kind="omim",
                method="graphic",
                combine="funSimAvg",
            )
        except Exception:
            return 0.0

    except ImportError:
        patient_set = set(patient_hpo_ids)
        disease_set = set(disease_hpo_ids)
        union = len(patient_set | disease_set)
        return len(patient_set & disease_set) / union if union > 0 else 0.0


def compute_combined_ic(hpo_term):
    """Return the combined Information Content of a pyhpo HPO term.

    Combined IC is the mean of the OMIM-based and gene-based IC when both are
    positive, otherwise whichever is positive, otherwise 1.0 (a default weight
    for terms with no IC). Any error reading the term's IC also yields 1.0.

    Args:
        hpo_term: A pyhpo HPOTerm object.

    Returns:
        float: The combined information content.
    """
    try:
        omim_ic = hpo_term.information_content.omim
        gene_ic = hpo_term.information_content.gene
        if omim_ic > 0 and gene_ic > 0:
            return (omim_ic + gene_ic) / 2.0
        if omim_ic > 0:
            return omim_ic
        if gene_ic > 0:
            return gene_ic
        return 1.0
    except Exception:
        return 1.0


def compute_pairwise_shared_counts(hpo_set1, hpo_set2):
    """Count HPO terms shared by two patients, raw and IC-weighted.

    The raw count is the size of the intersection. The IC-weighted count sums
    each shared term's combined information content; without pyhpo (or for terms
    that cannot be resolved) each shared term contributes a flat weight of 1.0,
    so the IC-weighted count equals the raw count in the Jaccard fallback.

    Args:
        hpo_set1: Iterable of patient 1's HPO term ids.
        hpo_set2: Iterable of patient 2's HPO term ids.

    Returns:
        tuple: (raw_count, ic_weighted_count) as (int, float).
    """
    shared_terms = set(hpo_set1) & set(hpo_set2)
    raw_count = len(shared_terms)

    if not ensure_ontology():
        return raw_count, float(raw_count)

    from pyhpo import Ontology

    ic_weighted_count = 0.0
    for term_id in shared_terms:
        try:
            ic_weighted_count += compute_combined_ic(Ontology.get_hpo_object(term_id))
        except Exception:
            ic_weighted_count += 1.0
    return raw_count, ic_weighted_count


def score_patient_vs_gene_filtered(patient_hpo_ids, gene_symbol, gene_disease_data, allowed_inheritance):
    """Score a patient's HPO terms against a gene's inheritance-filtered diseases.

    Diseases of ``gene_symbol`` are considered only when they carry phenotype
    terms and their inheritance is compatible with ``allowed_inheritance``. A
    disease with NO inheritance annotation is permissively included; a disease
    with inheritance is kept only if it intersects ``allowed_inheritance``. The
    best (highest-similarity) qualifying disease is reported.

    Args:
        patient_hpo_ids: Iterable of the patient's HPO term ids.
        gene_symbol: The gene symbol to score against, or None.
        gene_disease_data: Mapping
            ``{gene_symbol: {disease_id: {"phenotypes": set, "inheritance": set}}}``.
        allowed_inheritance: Set of inheritance codes compatible with the outlier
            type (e.g. {"AD", "XD"}).

    Returns:
        dict with keys: has_annotation (bool), best_disease (str or None),
        best_similarity (float), best_inheritance (str or None),
        overlap_count (int), n_matching_diseases (int).
    """
    if gene_symbol not in gene_disease_data:
        return {
            "has_annotation": False,
            "best_disease": None,
            "best_similarity": 0.0,
            "best_inheritance": None,
            "overlap_count": 0,
            "n_matching_diseases": 0,
        }

    disease_scores = []
    for disease_id, data in gene_disease_data[gene_symbol].items():
        disease_phenotypes = data["phenotypes"]
        disease_inheritance = data["inheritance"]
        if not disease_phenotypes:
            continue
        if disease_inheritance and not (disease_inheritance & allowed_inheritance):
            continue
        disease_scores.append({
            "disease_id": disease_id,
            "similarity": compute_similarity_pyhpo(patient_hpo_ids, disease_phenotypes),
            "overlap_count": len(set(patient_hpo_ids) & disease_phenotypes),
            "inheritance": ",".join(sorted(disease_inheritance)) if disease_inheritance else None,
        })

    if not disease_scores:
        return {
            "has_annotation": True,
            "best_disease": None,
            "best_similarity": 0.0,
            "best_inheritance": None,
            "overlap_count": 0,
            "n_matching_diseases": 0,
        }

    disease_scores.sort(key=lambda score: score["similarity"], reverse=True)
    return {
        "has_annotation": True,
        "best_disease": disease_scores[0]["disease_id"],
        "best_similarity": disease_scores[0]["similarity"],
        "best_inheritance": disease_scores[0]["inheritance"],
        "overlap_count": disease_scores[0]["overlap_count"],
        "n_matching_diseases": len(disease_scores),
    }


def get_qualifying_samples(row, outlier_type, affected_lookup, analysis_lookup, participant_to_hpo):
    """Select the qualifying outlier samples for one locus and outlier type.

    Walking the outlier entries largest-allele-first, a sample qualifies only if
    it passes all four gates:
      1. Affected and unsolved (``is_affected_unsolved``).
      2. Above the largest unaffected/solved sample (``is_above_unaffected``).
      3. Above every available population 99th percentile (``is_above_population_p99``).
      4. Has at least one HPO term in ``participant_to_hpo``.

    The result is then deduplicated by sample_id keeping the LARGEST qualifying
    allele per sample (its first occurrence in the descending-allele order). A
    diploid/homozygous sample can otherwise pass at two allele sizes, which would
    inflate the qualifying count and fabricate a perfect self-similarity in the
    consecutive-pair comparison.

    Args:
        row: A locus record dict carrying OutlierSampleIds_*, the
            FirstUnaffectedAlleleSize_* and *_99thPercentile columns, and LocusId.
        outlier_type: One of AllAlleles / ShortAlleles / HemizygousAlleles.
        affected_lookup: Mapping of sample_id -> affected_status.
        analysis_lookup: Mapping of sample_id -> analysis_status.
        participant_to_hpo: Mapping of sample_id -> set of HPO term ids.

    Returns:
        list of (sample_id, allele_size, hpo_term_set) tuples sorted by allele
        size descending, at most one entry per sample_id.
    """
    column = f"OutlierSampleIds_{outlier_type}"
    outlier_value = row.get(column)
    if is_missing_outlier_value(outlier_value):
        return []

    qualifying = []
    for allele_size, sample_id, _purity, _methylation in parse_outlier_entries(
        outlier_value, locus_id=row.get("LocusId"), column=column
    ):
        if not is_affected_unsolved(sample_id, affected_lookup, analysis_lookup):
            continue
        if not is_above_unaffected(allele_size, row, outlier_type):
            continue
        if not is_above_population_p99(allele_size, row):
            continue
        hpo_term_set = participant_to_hpo.get(sample_id, set())
        if not hpo_term_set:
            continue
        qualifying.append((sample_id, allele_size, set(hpo_term_set)))

    seen_sample_ids = set()
    deduped_qualifying = []
    for entry in qualifying:
        if entry[0] in seen_sample_ids:
            continue
        seen_sample_ids.add(entry[0])
        deduped_qualifying.append(entry)
    return deduped_qualifying


def compute_phenotype_scores(records, participant_to_hpo, gene_lookup, gene_disease_data,
                             affected_lookup, analysis_lookup):
    """Compute per-outlier and per-locus phenotype scores; augment records in place.

    For each locus record and each outlier type, the qualifying samples are
    scored two ways: a gene-phenotype similarity per sample (best matching,
    inheritance-filtered disease) and a pairwise similarity / shared-term count
    against the next qualifying sample in descending-allele order. Per-locus rows
    aggregate these (sum of pairwise similarities, max gene similarity, the
    comma-joined qualifying sample ids, plus denormalized filter columns).

    As a side effect, each record gains ``MaxGenePhenoSim_{outlier_type}`` and
    ``SumPairwiseSim_{outlier_type}`` columns (the per-locus aggregates), so the
    loci table carries the same denormalized phenotype columns the original
    pipeline merged in. Outlier types with no qualifying samples leave those
    columns as None on that record.

    If ``participant_to_hpo`` is empty there are no phenotypes to score: the
    function returns ([], []) and adds no columns (the loci-table phenotype
    columns stay NULL), so the build degrades gracefully.

    Args:
        records: List of locus record dicts (one per locus), mutated in place to
            add the MaxGenePhenoSim_*/SumPairwiseSim_* columns.
        participant_to_hpo: Mapping of sample_id -> set of HPO term ids.
        gene_lookup: Mapping of gene_id -> dict carrying at least ``gene_symbol``.
        gene_disease_data: Mapping
            ``{gene_symbol: {disease_id: {"phenotypes": set, "inheritance": set}}}``.
        affected_lookup: Mapping of sample_id -> affected_status.
        analysis_lookup: Mapping of sample_id -> analysis_status.

    Returns:
        tuple: (per_outlier_rows, per_locus_rows), each a list of dicts matching
        the per_outlier_phenotype_scores / per_locus_phenotype_scores schemas.
    """
    if not participant_to_hpo:
        return [], []

    ensure_ontology()

    per_outlier_rows = []
    per_locus_rows = []

    for row in records:
        locus_id = row.get("LocusId")
        gene_id = row.get("gene_id")
        gene_symbol = None
        if gene_id and gene_id in gene_lookup:
            gene_symbol = gene_lookup[gene_id].get("gene_symbol")

        for outlier_type in OUTLIER_TYPES:
            qualifying = get_qualifying_samples(
                row, outlier_type, affected_lookup, analysis_lookup, participant_to_hpo
            )
            if not qualifying:
                continue

            gene_similarities = []
            pairwise_similarities = []
            pairwise_shared_raw = []
            pairwise_shared_ic = []

            for index, (sample_id, allele_size, hpo_terms) in enumerate(qualifying):
                if gene_symbol and gene_symbol in gene_disease_data:
                    score_result = score_patient_vs_gene_filtered(
                        hpo_terms, gene_symbol, gene_disease_data,
                        OUTLIER_TYPE_INHERITANCE[outlier_type],
                    )
                    gene_similarity = score_result["best_similarity"]
                    overlap_count = score_result["overlap_count"]
                    n_matching_diseases = score_result["n_matching_diseases"]
                    best_disease = score_result["best_disease"]
                    best_inheritance = score_result["best_inheritance"]
                    gene_similarities.append(gene_similarity)
                else:
                    gene_similarity = None
                    overlap_count = None
                    n_matching_diseases = None
                    best_disease = None
                    best_inheritance = None

                if index < len(qualifying) - 1:
                    next_sample_id, _next_allele, next_hpo = qualifying[index + 1]
                    pairwise_similarity = compute_similarity_pyhpo(hpo_terms, next_hpo)
                    pairwise_similarities.append(pairwise_similarity)
                    raw_count, ic_count = compute_pairwise_shared_counts(hpo_terms, next_hpo)
                    pairwise_shared_raw.append(raw_count)
                    pairwise_shared_ic.append(ic_count)
                else:
                    next_sample_id = None
                    pairwise_similarity = None
                    raw_count = None
                    ic_count = None

                per_outlier_rows.append({
                    "locus_id": locus_id,
                    "sample_id": sample_id,
                    "outlier_type": outlier_type,
                    "allele_size": allele_size,
                    "gene_symbol": gene_symbol,
                    "gene_phenotype_similarity": gene_similarity,
                    "gene_phenotype_overlap_count": overlap_count,
                    "n_matching_diseases": n_matching_diseases,
                    "best_matching_disease": best_disease,
                    "best_disease_inheritance": best_inheritance,
                    "pairwise_similarity_to_next": pairwise_similarity,
                    "pairwise_shared_count_raw": raw_count,
                    "pairwise_shared_count_ic": ic_count,
                    "next_sample_id": next_sample_id,
                })

            max_gene_similarity = max(gene_similarities) if gene_similarities else None
            sum_pairwise_similarity = sum(pairwise_similarities) if pairwise_similarities else None

            row[f"MaxGenePhenoSim_{outlier_type}"] = max_gene_similarity
            row[f"SumPairwiseSim_{outlier_type}"] = sum_pairwise_similarity

            per_locus_rows.append({
                "locus_id": locus_id,
                "outlier_type": outlier_type,
                "num_qualifying_samples": len(qualifying),
                "sum_pairwise_similarity": sum_pairwise_similarity,
                "sum_pairwise_shared_raw": sum(pairwise_shared_raw) if pairwise_shared_raw else None,
                "sum_pairwise_shared_ic": sum(pairwise_shared_ic) if pairwise_shared_ic else None,
                "max_gene_phenotype_similarity": max_gene_similarity,
                "qualifying_sample_ids": ",".join(entry[0] for entry in qualifying),
                "IsKnownMotif": row.get("IsKnownMotif"),
                "gene_region": row.get("gene_region"),
                "gene_region_rank": row.get("gene_region_rank"),
                "FirstAffectedAlleleSize": row.get(f"FirstAffectedAlleleSize_{outlier_type}"),
                "FirstUnaffectedAlleleSize": row.get(f"FirstUnaffectedAlleleSize_{outlier_type}"),
                "NumRepeatsInReference": row.get("NumRepeatsInReference"),
                "HPRC256_MaxAllele": row.get("HPRC256_MaxAllele"),
                "AoU1027_MaxAllele": row.get("AoU1027_MaxAllele"),
                "TenK10K_MaxAllele": row.get("TenK10K_MaxAllele"),
                "NumAffectedAboveUnaffected": row.get(f"NumAffectedUnsolvedSamplesAboveUnaffected_{outlier_type}"),
                "NumAffectedFamiliesAboveUnaffected": row.get(f"NumAffectedUnsolvedFamiliesAboveUnaffected_{outlier_type}"),
                "MotifSize": row.get("MotifSize"),
            })

    return per_outlier_rows, per_locus_rows
