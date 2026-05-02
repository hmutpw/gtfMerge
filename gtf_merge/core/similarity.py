"""
Intron-based similarity scoring for transcript merging.

Core philosophy:
  - Introns are the fundamental unit (not junctions/exons)
  - Introns are independent of TSS/TES noise
  - 5' truncation is common → default no penalty
  - 3' truncation is rare (polyA) → light penalty
  - Intron retention is severe → heavy penalty (50pts)
  - Exon skipping = different isoform → never merge (100pts)
  - Coordinate drift → penalise per bp
  - Score is continuous [0,100]; user sets threshold
"""
from typing import List, Tuple, Optional
from gtf_merge.core.models import Intron, Transcript, MergeScore, MergeType


# ---------------------------------------------------------------------------
# Step 1: find matching intron pairs
# ---------------------------------------------------------------------------

def match_introns(
    introns_a: List[Intron],
    introns_b: List[Intron],
    threshold: int,
) -> List[Tuple[int, int, int]]:
    """
    Greedy nearest-match between two intron lists.
    Returns [(idx_a, idx_b, coord_error), ...]
    Preserves order — critical for continuity check.
    """
    matched = []
    used_b  = set()
    for i, ia in enumerate(introns_a):
        best_j   = -1
        best_err = threshold + 1
        for j, ib in enumerate(introns_b):
            if j in used_b:
                continue
            if not ia.matches(ib, threshold):
                continue
            err = ia.coord_error(ib)
            if err < best_err:
                best_err = err
                best_j   = j
        if best_j >= 0:
            matched.append((i, best_j, best_err))
            used_b.add(best_j)
    return matched


# ---------------------------------------------------------------------------
# Step 2: continuity check
# ---------------------------------------------------------------------------

def is_continuous(indices: List[int]) -> bool:
    """True if indices form an unbroken sequence."""
    if len(indices) <= 1:
        return True
    return indices == list(range(indices[0], indices[-1] + 1))


def check_internal_events(
    a_indices: List[int],
    b_indices: List[int],
    len_a:     int,
    len_b:     int,
) -> Tuple[bool, bool]:
    """
    Returns (has_retention, has_skipping).
    retention: matched indices in A are non-continuous (gap in A → intron retained in B)
    skipping:  matched indices in B are non-continuous (gap in B → exon skipped in B)
    """
    has_retention = not is_continuous(a_indices)  # gap in A = retention in B
    has_skipping  = not is_continuous(b_indices)  # gap in B = skipping in B
    return has_retention, has_skipping


# ---------------------------------------------------------------------------
# Step 3: classify truncation type
# ---------------------------------------------------------------------------

def classify_truncation(
    a_indices: List[int],
    b_indices: List[int],
    len_a: int,
    len_b: int,
) -> Tuple[int, int]:
    """
    Returns (missing_5, missing_3) from B's perspective relative to A.
    missing_5 = how many introns A has before the first match
    missing_3 = how many introns A has after the last match
    """
    missing_5 = a_indices[0]
    missing_3 = (len_a - 1) - a_indices[-1]
    return missing_5, missing_3


# ---------------------------------------------------------------------------
# Step 4: compute merge score
# ---------------------------------------------------------------------------

def compute_merge_score(
    a: Transcript,
    b: Transcript,
    threshold:      int,
    per_bp_penalty: float = 1.0,
    trunc5_penalty: float = 0.0,
    trunc3_penalty: float = 5.0,
    retention_penalty: float = 50.0,
    min_shared_introns: int = 1,
) -> MergeScore:
    """
    Compute intron-based merge score between two multi-exon transcripts.

    Scoring:
      Base score = 100
      - coord errors:    max(donor_err, acceptor_err) × per_bp_penalty  per intron
      - 5' truncation:   missing_5 × trunc5_penalty  (default 0)
      - 3' truncation:   missing_3 × trunc3_penalty  (default 5)
      - intron retention: retention_count × retention_penalty (default 50)
      - exon skipping:   100 per event (= not mergeable)

    Returns MergeScore with .score property and .is_mergeable flag.
    """
    ms = MergeScore()

    ia = a.introns
    ib = b.introns

    if not ia or not ib:
        return ms  # monoexon handled separately

    # --- match ---
    matched = match_introns(ia, ib, threshold)
    ms.matched_introns = len(matched)

    if ms.matched_introns < min_shared_introns:
        return ms  # too few shared introns

    a_indices = [m[0] for m in matched]
    b_indices = [m[1] for m in matched]

    # --- internal events (retention / skipping) ---
    has_retention, has_skipping = check_internal_events(
        a_indices, b_indices, len(ia), len(ib)
    )

    if has_skipping:
        # gap in b_indices: B itself has discontinuous exon structure
        ms.skipping_penalty = 100.0
        ms.merge_type       = MergeType.SINGLETON
        ms.is_mergeable     = False
        return ms

    if has_retention:
        # gap in a_indices: A has internal introns that B doesn't have
        # → B retains those as exon (intron retention in B)
        n_retention = (a_indices[-1] - a_indices[0] + 1) - len(a_indices)
        ms.retention_penalty = n_retention * retention_penalty
        ms.merge_type        = MergeType.SINGLETON
        ms.is_mergeable      = False   # retention = different isoform, never merge
        return ms

    # --- coordinate penalties ---
    errors = [m[2] for m in matched]
    ms.per_intron_errors = errors
    ms.coord_penalty = sum(e * per_bp_penalty for e in errors)

    # --- truncation ---
    missing_5, missing_3 = classify_truncation(a_indices, b_indices, len(ia), len(ib))
    ms.missing_5 = missing_5
    ms.missing_3 = missing_3
    ms.trunc5_penalty = missing_5 * trunc5_penalty
    ms.trunc3_penalty = missing_3 * trunc3_penalty

    # --- merge type ---
    has_coord_error = any(e > 0 for e in errors)
    if missing_5 == 0 and missing_3 == 0:
        ms.merge_type = MergeType.FUZZY_INTRON if has_coord_error else MergeType.EXACT_MATCH
    elif missing_5 > 0 and missing_3 == 0:
        ms.merge_type = MergeType.FUZZY_TRUNCATION_5 if has_coord_error else MergeType.TRUNCATION_5
    elif missing_5 == 0 and missing_3 > 0:
        ms.merge_type = MergeType.FUZZY_TRUNCATION_3 if has_coord_error else MergeType.TRUNCATION_3
    else:
        ms.merge_type = MergeType.TRUNCATION_BOTH

    ms.is_mergeable = True
    return ms


# ---------------------------------------------------------------------------
# Mono-exon similarity (unchanged logic, renamed)
# ---------------------------------------------------------------------------

def monoexon_similarity(
    a: Transcript,
    b: Transcript,
    mono_overlap: float,
    mono_3end:    int,
    mono_5end:    int,
) -> Tuple[bool, str]:
    """Check if two mono-exon transcripts should be merged."""
    if a.strand != b.strand and a.strand != "." and b.strand != ".":
        return False, ""
    ea = a.exons[0] if a.exons else None
    eb = b.exons[0] if b.exons else None
    if ea is None or eb is None:
        return False, ""

    ov_start = max(ea.start, eb.start)
    ov_end   = min(ea.end,   eb.end)
    overlap  = max(0, ov_end - ov_start)
    shorter  = min(ea.length, eb.length)
    if shorter == 0:
        return False, ""

    overlap_ratio = overlap / shorter
    diff_5 = abs(ea.start - eb.start)
    diff_3 = abs(ea.end   - eb.end)

    if diff_5 == 0 and diff_3 == 0:
        return True, MergeType.MONOEXON_EXACT

    if a.strand == "+":
        three_diff, five_diff = diff_3, diff_5
    elif a.strand == "-":
        three_diff, five_diff = diff_5, diff_3
    else:
        three_diff = min(diff_3, diff_5)
        five_diff  = max(diff_3, diff_5)

    if (overlap_ratio >= mono_overlap and
            three_diff <= mono_3end and
            five_diff  <= mono_5end):
        return True, MergeType.MONOEXON_OVERLAP

    return False, ""
