"""
gffcompare-style class code assignment.

For each query transcript, finds the best-matching reference transcript
and assigns one of these codes:
  =  exact match
  c  query contained in reference (junction-wise)
  k  reference contained in query
  m  retained intron (all)
  n  partial retained intron
  j  multi-exon, ≥1 matched junction (novel isoform)
  e  single-exon overlap with multi-exon ref
  o  generic overlap, same strand
  s  intron match on opposite strand
  x  exonic overlap on opposite strand
  i  query fully within reference intron
  y  reference fully within query
  p  possible polymerase run-on
  u  unknown / intergenic

Algorithm summary (mirrors gffcompare):
  - For each query, scan all overlapping ref transcripts
  - Pick the code with highest priority across all ref hits
  - Lower priority number = better match (= is best, u is worst)
"""
import logging
from typing import List, Tuple, Optional
from collections import defaultdict
from intervaltree import IntervalTree

from gtf_merge.core.models import (
    Transcript, MergedTranscript, Junction,
    ClassCode, CLASS_CODE_PRIORITY,
)

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────
# Reference annotation index
# ─────────────────────────────────────────────────────────────────

class RefIndex:
    """Indexes reference transcripts for fast overlap queries."""

    def __init__(self):
        # {(chrom, strand): IntervalTree}
        self.trees: dict = {}
        # {transcript_id: ref_dict}
        self.by_id: dict = {}

    def add_transcript(self, t: Transcript):
        key = (t.chrom, t.strand)
        if key not in self.trees:
            self.trees[key] = IntervalTree()
        self.trees[key][t.start:t.end] = t.transcript_id
        self.by_id[t.transcript_id] = {
            "transcript": t,
            "gene_id":   t.gene_id,
            "gene_name": t.attributes.get("gene_name", "NA"),
        }

    def overlapping(self, chrom: str, start: int, end: int,
                    strand: Optional[str] = None) -> List[Transcript]:
        """Return ref transcripts overlapping the given interval.
        If strand=None, returns both + and - strand hits."""
        results = []
        if strand:
            keys = [(chrom, strand)]
        else:
            keys = [(chrom, "+"), (chrom, "-"), (chrom, ".")]
        for key in keys:
            tree = self.trees.get(key)
            if tree is None:
                continue
            for iv in tree[start:end]:
                results.append(self.by_id[iv.data]["transcript"])
        return results


# ─────────────────────────────────────────────────────────────────
# Geometric helpers
# ─────────────────────────────────────────────────────────────────

def has_genomic_overlap(a, b) -> bool:
    return a.chrom == b.chrom and a.start < b.end and a.end > b.start


def has_exon_overlap(a, b) -> bool:
    """Check if any exon of a overlaps any exon of b."""
    if a.chrom != b.chrom:
        return False
    for ea in a.exons:
        for eb in b.exons:
            if ea.start < eb.end and ea.end > eb.start:
                return True
    return False


def junctions_match(query_juncs, ref_juncs, threshold: int) -> List[Tuple[int, int]]:
    """Return list of (query_idx, ref_idx) pairs of matching junctions."""
    pairs = []
    used_r = set()
    for qi, qj in enumerate(query_juncs):
        for ri, rj in enumerate(ref_juncs):
            if ri in used_r:
                continue
            if qj.matches(rj, threshold):
                pairs.append((qi, ri))
                used_r.add(ri)
                break
    return pairs


def all_junctions_match(query_juncs, ref_juncs, threshold: int) -> bool:
    """All query junctions must find a match in ref (allows ref to have extras)."""
    if not query_juncs:
        return False
    pairs = junctions_match(query_juncs, ref_juncs, threshold)
    return len(pairs) == len(query_juncs)


def exact_chain_match(q_juncs, r_juncs, threshold: int) -> bool:
    """Exact junction chain match: same number, all aligned."""
    if len(q_juncs) != len(r_juncs):
        return False
    if not q_juncs:
        return False
    for q, r in zip(q_juncs, r_juncs):
        if not q.matches(r, threshold):
            return False
    return True


def contained_in(inner, outer) -> bool:
    """All exons of inner are contained within outer's exonic span."""
    if inner.chrom != outer.chrom:
        return False
    return inner.start >= outer.start and inner.end <= outer.end


def in_intron(query, ref) -> bool:
    """Query falls entirely within one of ref's introns."""
    if query.chrom != ref.chrom:
        return False
    for j in ref.junctions:
        # intron span = donor to acceptor (exclusive)
        if query.start >= j.donor and query.end <= j.acceptor:
            return True
    return False


# ─────────────────────────────────────────────────────────────────
# Single-pair classification
# ─────────────────────────────────────────────────────────────────

def classify_pair(query, ref, threshold: int) -> Optional[ClassCode]:
    """
    Classify the relationship between query and ref.
    Returns None if no overlap (caller will set 'u' if no any ref matches).
    """
    if not has_genomic_overlap(query, ref):
        return None

    same_strand = (query.strand == ref.strand) or query.strand == "." or ref.strand == "."

    # ── opposite strand ─────────────────────────────────────────
    if not same_strand:
        # 's' — intron match on opposite strand (strand-blind comparison)
        if query.junctions and ref.junctions:
            for qj in query.junctions:
                for rj in ref.junctions:
                    if (qj.chrom == rj.chrom
                            and abs(qj.donor - rj.donor) <= threshold
                            and abs(qj.acceptor - rj.acceptor) <= threshold):
                        return ClassCode.INTRON_MATCH_OPP
        # 'x' — exonic overlap on opposite strand
        if has_exon_overlap(query, ref):
            return ClassCode.EXONIC_OPP
        return None

    # ── same strand ─────────────────────────────────────────────

    # 'i' — query fully in ref intron (no exonic overlap)
    if in_intron(query, ref) and not has_exon_overlap(query, ref):
        return ClassCode.INTRONIC

    # if no exonic overlap on same strand, no class
    if not has_exon_overlap(query, ref):
        return None

    # mono-exon special cases
    if query.is_monoexon:
        # 'e' — single-exon overlap with multi-exon ref
        if not ref.is_monoexon:
            # check if it's truly contained in an exon
            qe = query.exons[0]
            for re in ref.exons:
                if qe.start >= re.start and qe.end <= re.end:
                    return ClassCode.CONTAINED   # 'c' if fully inside an exon
            return ClassCode.SINGLE_OVERLAP
        # both mono-exon: classify by overlap geometry
        qe = query.exons[0]
        re = ref.exons[0]
        if qe.start >= re.start and qe.end <= re.end:
            if qe.start == re.start and qe.end == re.end:
                return ClassCode.EQUAL
            return ClassCode.CONTAINED
        if re.start >= qe.start and re.end <= qe.end:
            return ClassCode.CONTAINING
        return ClassCode.GENERIC_OVERLAP

    # mono-exon ref vs multi-exon query
    if ref.is_monoexon:
        qe_exons = query.exons
        re = ref.exons[0]
        for qe in qe_exons:
            if re.start >= qe.start and re.end <= qe.end:
                return ClassCode.CONTAINING  # ref contained in one query exon
        if has_exon_overlap(query, ref):
            return ClassCode.GENERIC_OVERLAP
        return None

    # both multi-exon — check chain match
    if exact_chain_match(query.junctions, ref.junctions, threshold):
        return ClassCode.EQUAL

    # check 'c' — query contained: all query junctions match a contiguous subset of ref
    if all_junctions_match(query.junctions, ref.junctions, threshold) and \
       len(query.junctions) < len(ref.junctions):
        return ClassCode.CONTAINED

    # check 'k' — ref contained in query
    if all_junctions_match(ref.junctions, query.junctions, threshold) and \
       len(ref.junctions) < len(query.junctions):
        return ClassCode.CONTAINING

    # check 'y' — ref fully within query (incl introns)
    if contained_in(ref, query) and not contained_in(query, ref):
        return ClassCode.REF_WITHIN_QUERY

    # check retention — query has fewer junctions because some ref introns are exonic
    pairs = junctions_match(query.junctions, ref.junctions, threshold)
    if pairs and len(query.junctions) < len(ref.junctions):
        # some ref junctions are missing from query (potentially retained)
        q_idx = [p[0] for p in pairs]
        r_idx = [p[1] for p in pairs]
        # retention: ref indices have gaps
        if r_idx != list(range(r_idx[0], r_idx[-1] + 1)):
            # are the missing intron regions exonic in query?
            missing_count = (r_idx[-1] - r_idx[0] + 1) - len(r_idx)
            if missing_count >= 1:
                # if all missing introns are within a query exon → 'm', else 'n'
                return ClassCode.RETAINED_INTRON

    # 'j' — at least one matched junction (novel isoform)
    if pairs:
        return ClassCode.JUNCTION_MATCH

    # 'p' — possible run-on: adjacent but no exon overlap (already handled by overlap check)
    # if exonic overlap but no junction match
    return ClassCode.GENERIC_OVERLAP


# ─────────────────────────────────────────────────────────────────
# Top-level: classify a transcript against a ref index
# ─────────────────────────────────────────────────────────────────

def classify_transcript(query, ref_index: RefIndex, threshold: int) -> Tuple[ClassCode, Optional[Transcript]]:
    """
    Find the best class code for query among all overlapping refs.
    Returns (best_code, best_ref_transcript_or_None).
    """
    refs = ref_index.overlapping(query.chrom, query.start, query.end)
    if not refs:
        return ClassCode.UNKNOWN, None

    best_code = ClassCode.UNKNOWN
    best_ref  = None
    best_priority = CLASS_CODE_PRIORITY[ClassCode.UNKNOWN.value]

    for ref in refs:
        code = classify_pair(query, ref, threshold)
        if code is None:
            continue
        priority = CLASS_CODE_PRIORITY[code.value]
        if priority < best_priority:
            best_priority = priority
            best_code = code
            best_ref  = ref
        if best_code == ClassCode.EQUAL:
            break  # can't get better

    return best_code, best_ref
