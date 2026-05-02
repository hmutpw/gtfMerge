"""
Step2: Cluster - Group structurally similar transcripts (no forced merge).

Three-dimensional similarity:
  similarity = w1 * structural_score + w2 * junction_score + w3 * jaccard_score

Clustering: cd-hit style with reference-first / sample-count / length ordering.

Key features:
  - Only computes similarity for transcripts with genomic overlap + same strand
  - Allows retention/skipping pairs (similarity may be lower but informative)
  - Each transcript belongs to exactly one cluster
  - Reference transcripts are preferred as cluster representatives
"""
import logging
from collections import defaultdict
from typing import List, Dict, Tuple, Optional
from intervaltree import IntervalTree

from gtf_merge.core.models import (
    Transcript, MergedTranscript, Cluster, ClusterMember,
    SimilarityScore, ClassCode,
)

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────
# Similarity dimensions
# ─────────────────────────────────────────────────────────────────

def structural_score(a: Transcript, b: Transcript) -> float:
    """Similarity of basic structural features: exon count + transcript length."""
    n_a = len(a.exons)
    n_b = len(b.exons)
    if max(n_a, n_b) == 0:
        return 0.0
    exon_sim = 1.0 - abs(n_a - n_b) / max(n_a, n_b)

    len_a = a.length
    len_b = b.length
    if max(len_a, len_b) == 0:
        return 0.0
    len_sim = 1.0 - abs(len_a - len_b) / max(len_a, len_b)

    return exon_sim * len_sim


def junction_score(
    a: Transcript,
    b: Transcript,
    threshold:       int,
    fuzzy_threshold: int,
) -> float:
    """
    Junction match score over the union of junctions.
    Each junction contributes:
      1.0 if matched within threshold
      0.5 if matched within fuzzy_threshold but outside threshold
      0.0 otherwise
    """
    if not a.junctions and not b.junctions:
        return 1.0   # both mono-exon → don't penalise
    if not a.junctions or not b.junctions:
        return 0.0   # one is mono-exon, the other isn't

    # find strict matches
    used_b_strict = set()
    strict_matches = 0
    for ja in a.junctions:
        for k, jb in enumerate(b.junctions):
            if k in used_b_strict:
                continue
            if ja.matches(jb, threshold):
                used_b_strict.add(k)
                strict_matches += 1
                break

    # find fuzzy matches (in remaining b junctions)
    used_b_fuzzy = set(used_b_strict)
    fuzzy_matches = 0
    for ja in a.junctions:
        # skip those already strict-matched
        if any(ja.matches(b.junctions[k], threshold) for k in used_b_strict):
            continue
        for k, jb in enumerate(b.junctions):
            if k in used_b_fuzzy:
                continue
            if ja.matches(jb, fuzzy_threshold):
                used_b_fuzzy.add(k)
                fuzzy_matches += 1
                break

    union_size = len(a.junctions) + len(b.junctions) - strict_matches
    score = (strict_matches * 1.0 + fuzzy_matches * 0.5) / max(union_size, 1)
    return min(score, 1.0)


def jaccard_score(a: Transcript, b: Transcript) -> float:
    """Base-level Jaccard: |exonic_A ∩ exonic_B| / |exonic_A ∪ exonic_B|."""
    if a.chrom != b.chrom:
        return 0.0

    # collect exonic intervals
    intervals_a = [(e.start, e.end) for e in a.exons]
    intervals_b = [(e.start, e.end) for e in b.exons]

    overlap = 0
    for sa, ea in intervals_a:
        for sb, eb in intervals_b:
            ov = max(0, min(ea, eb) - max(sa, sb))
            overlap += ov

    len_a = sum(e - s for s, e in intervals_a)
    len_b = sum(e - s for s, e in intervals_b)
    union = len_a + len_b - overlap

    if union == 0:
        return 0.0
    return overlap / union


# ─────────────────────────────────────────────────────────────────
# Combined similarity
# ─────────────────────────────────────────────────────────────────

def compute_similarity(
    a: Transcript,
    b: Transcript,
    threshold:       int,
    fuzzy_threshold: int,
    w_structural: float,
    w_junction:   float,
    w_jaccard:    float,
) -> SimilarityScore:
    """Compute combined similarity. Both transcripts assumed same chrom+strand and overlapping."""

    s_struct = structural_score(a, b)
    s_junc   = junction_score(a, b, threshold, fuzzy_threshold)
    s_jacc   = jaccard_score(a, b)

    # mono-exon: junction is N/A → redistribute weights
    if (not a.junctions) and (not b.junctions):
        # both mono-exon: structural + jaccard only
        total_w = w_structural + w_jaccard
        if total_w == 0:
            total = 0.0
        else:
            ws = w_structural / total_w
            wj = w_jaccard / total_w
            total = ws * s_struct + wj * s_jacc
    else:
        total_w = w_structural + w_junction + w_jaccard
        if total_w == 0:
            total = 0.0
        else:
            ws = w_structural / total_w
            wjn = w_junction / total_w
            wjc = w_jaccard / total_w
            total = ws * s_struct + wjn * s_junc + wjc * s_jacc

    return SimilarityScore(
        structural = round(s_struct, 4),
        junction   = round(s_junc, 4),
        jaccard    = round(s_jacc, 4),
        total      = round(total, 4),
    )


# ─────────────────────────────────────────────────────────────────
# Index for fast overlap lookup
# ─────────────────────────────────────────────────────────────────

class TranscriptIndex:
    """Indexes transcripts by chrom+strand for overlap queries."""

    def __init__(self):
        self.trees: Dict[Tuple[str, str], IntervalTree] = {}
        self.by_id: Dict[str, Transcript] = {}

    def add(self, t: Transcript):
        key = (t.chrom, t.strand)
        if key not in self.trees:
            self.trees[key] = IntervalTree()
        self.trees[key][t.start:t.end] = t.transcript_id
        self.by_id[t.transcript_id] = t

    def overlapping(self, t: Transcript) -> List[Transcript]:
        key = (t.chrom, t.strand)
        tree = self.trees.get(key)
        if tree is None:
            return []
        return [self.by_id[iv.data] for iv in tree[t.start:t.end]
                if iv.data != t.transcript_id]


# ─────────────────────────────────────────────────────────────────
# cd-hit style clustering
# ─────────────────────────────────────────────────────────────────

def _representative_sort_key(
    t: Transcript,
    is_reference: callable,
    sample_counts: Dict[str, int],
):
    """
    Sort key for representative selection (smaller = higher priority).
    Order: reference > sample count > length > id (for reproducibility)
    """
    is_ref = is_reference(t.transcript_id)
    samples = sample_counts.get(t.transcript_id, 0)
    return (
        not is_ref,           # ref first
        -samples,             # more samples first
        -t.length,            # longer first
        t.transcript_id,      # alphabetical
    )


def cluster_transcripts(
    transcripts:     List[Transcript],
    ref_ids:         set,                      # set of transcript_ids that are reference
    sample_counts:   Dict[str, int],           # transcript_id -> sample count
    junction_threshold: int,
    fuzzy_threshold: int,
    w_structural:    float,
    w_junction:      float,
    w_jaccard:       float,
    min_similarity:  float,
) -> List[Cluster]:
    """
    cd-hit style clustering.
    Returns list of Cluster objects.
    """
    if not transcripts:
        return []

    is_ref = lambda tid: tid in ref_ids

    # group by chrom + strand for performance
    cs_groups: Dict[Tuple[str, str], List[Transcript]] = defaultdict(list)
    for t in transcripts:
        cs_groups[(t.chrom, t.strand)].append(t)

    clusters_all: List[Cluster] = []
    cluster_idx = 0

    for (chrom, strand), group in cs_groups.items():
        # sort for representative selection
        group.sort(key=lambda t: _representative_sort_key(t, is_ref, sample_counts))

        # interval tree of representatives so far
        rep_tree = IntervalTree()
        clusters_local: List[Cluster] = []

        for t in group:
            # find candidate clusters (overlapping representatives)
            assigned = False
            best_cluster = None
            best_sim     = SimilarityScore()

            for iv in rep_tree[t.start:t.end]:
                ci = iv.data
                rep = clusters_local[ci].representative
                rep_t = next(m for m in [rep] if m.transcript_id == rep.transcript_id)
                rep_transcript = clusters_local[ci]._rep_transcript  # stored on cluster

                sim = compute_similarity(
                    rep_transcript, t,
                    threshold        = junction_threshold,
                    fuzzy_threshold  = fuzzy_threshold,
                    w_structural     = w_structural,
                    w_junction       = w_junction,
                    w_jaccard        = w_jaccard,
                )
                if sim.total >= min_similarity and sim.total > best_sim.total:
                    best_sim     = sim
                    best_cluster = ci

            if best_cluster is not None:
                # add to existing cluster
                member = ClusterMember(
                    transcript_id    = t.transcript_id,
                    similarity       = best_sim,
                    is_representative= False,
                )
                clusters_local[best_cluster].members.append(member)
                assigned = True

            if not assigned:
                # create new cluster
                cluster_idx += 1
                cid = f"C{cluster_idx:06d}"
                rep_member = ClusterMember(
                    transcript_id     = t.transcript_id,
                    similarity        = SimilarityScore(structural=1.0, junction=1.0,
                                                          jaccard=1.0, total=1.0),
                    is_representative = True,
                )
                cluster = Cluster(
                    cluster_id     = cid,
                    chrom          = chrom,
                    strand         = strand,
                    representative = rep_member,
                    members        = [rep_member],
                )
                cluster._rep_transcript = t   # store for similarity computation
                clusters_local.append(cluster)
                rep_tree[t.start:t.end] = len(clusters_local) - 1

        clusters_all.extend(clusters_local)

    return clusters_all
