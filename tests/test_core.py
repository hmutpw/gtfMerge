"""Tests for new two-step architecture."""
import pytest
from gtf_merge.core.models import (
    Transcript, Exon, Junction,
    MergeStatus, ClassCode, CLASS_CODE_PRIORITY,
)
from gtf_merge.core.merger import (
    can_merge_multiexon, can_merge_monoexon,
    match_junctions, is_continuous,
)
from gtf_merge.core.class_code import (
    classify_pair, RefIndex, classify_transcript,
)
from gtf_merge.core.cluster import (
    structural_score, junction_score, jaccard_score,
    compute_similarity, cluster_transcripts,
)


# ─────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────

def make_t(juncs, strand="+", chrom="chr1", start=100, end=2000, sample="S1", tid=None):
    if not juncs:
        exons = [Exon(chrom=chrom, start=start, end=end, strand=strand)]
    else:
        exons = [Exon(chrom=chrom, start=start, end=juncs[0][0], strand=strand)]
        for k, (donor, acceptor) in enumerate(juncs):
            next_end = juncs[k+1][0] if k+1 < len(juncs) else end
            exons.append(Exon(chrom=chrom, start=acceptor, end=next_end, strand=strand))
    junctions = [Junction(chrom, d, a, strand) for d, a in juncs]
    return Transcript(
        transcript_id = tid or f"t_{sample}",
        gene_id       = "G1",
        chrom         = chrom,
        start         = start,
        end           = end,
        strand        = strand,
        source_sample = sample,
        source_file   = "test.gtf",
        exons         = exons,
        junctions     = junctions,
    )


# ─────────────────────────────────────────────────────────────────
# Step1: merge logic
# ─────────────────────────────────────────────────────────────────

class TestCanMergeMultiexon:

    def _kwargs(self, **overrides):
        defaults = dict(
            junction_threshold = 3,
            max_5trunc         = 1,
            max_3trunc         = 0,
            tss_window         = 200,
            tes_window         = 50,
            max_diff_count     = -1,
            max_diff_ratio     = 0.5,
        )
        defaults.update(overrides)
        return defaults

    def test_exact_match(self):
        A = make_t([(200,300),(500,600),(800,900)])
        B = make_t([(200,300),(500,600),(800,900)])
        ok, status = can_merge_multiexon(A, B, **self._kwargs())
        assert ok and status == "exact_match"

    def test_fuzzy_match(self):
        """Small drift in coords still merges (with relaxed diff_ratio)."""
        A = make_t([(200,300),(500,600)])
        B = make_t([(201,300),(500,602)])
        # 2/2 = 1.0 diff ratio, but small drift → allow with relaxed param
        ok, status = can_merge_multiexon(A, B, **self._kwargs(max_diff_ratio=1.0))
        assert ok and status == "fuzzy_match"

    def test_5prime_truncation_allowed(self):
        A = make_t([(200,300),(500,600),(800,900)])
        B = make_t([(500,600),(800,900)])
        ok, status = can_merge_multiexon(A, B, **self._kwargs(max_5trunc=1))
        assert ok and status == "truncation_5"

    def test_5prime_truncation_too_many(self):
        A = make_t([(200,300),(500,600),(800,900),(1100,1200)])
        B = make_t([(800,900),(1100,1200)])  # missing 2 from 5'
        ok, status = can_merge_multiexon(A, B, **self._kwargs(max_5trunc=1))
        assert not ok

    def test_3prime_truncation_disallowed(self):
        A = make_t([(200,300),(500,600),(800,900)])
        B = make_t([(200,300),(500,600)])
        ok, status = can_merge_multiexon(A, B, **self._kwargs(max_3trunc=0))
        assert not ok

    def test_3prime_truncation_allowed(self):
        A = make_t([(200,300),(500,600),(800,900)])
        B = make_t([(200,300),(500,600)])
        ok, status = can_merge_multiexon(A, B, **self._kwargs(max_3trunc=1))
        assert ok and status == "truncation_3"

    def test_internal_retention_blocked(self):
        """B is missing internal junction → must not merge."""
        A = make_t([(200,300),(500,600),(800,900)])
        B = make_t([(200,300),(800,900)])  # missing j2 internally
        ok, status = can_merge_multiexon(A, B, **self._kwargs())
        assert not ok and status == "internal_gap"

    def test_tss_window(self):
        A = make_t([(200,300),(500,600)], start=100)
        B = make_t([(200,300),(500,600)], start=500)  # TSS diff = 400
        ok, status = can_merge_multiexon(A, B, **self._kwargs(tss_window=200))
        assert not ok and status == "tss_too_far"

    def test_tes_window(self):
        A = make_t([(200,300),(500,600)], end=2000)
        B = make_t([(200,300),(500,600)], end=2200)
        ok, status = can_merge_multiexon(A, B, **self._kwargs(tes_window=50))
        assert not ok and status == "tes_too_far"

    def test_diff_ratio_too_high(self):
        """All junctions have drift → diff_ratio too high."""
        A = make_t([(200,300),(500,600),(800,900),(1100,1200)])
        B = make_t([(202,300),(502,600),(802,900),(1102,1200)])  # all drift
        # diff_ratio = 4/4 = 1.0 > 0.5
        ok, status = can_merge_multiexon(A, B, **self._kwargs(max_diff_ratio=0.5))
        assert not ok and status == "diff_ratio_exceeded"


class TestCanMergeMonoexon:
    def test_exact(self):
        A = make_t([], start=100, end=500)
        B = make_t([], start=100, end=500)
        ok, status = can_merge_monoexon(A, B, 0.5, 500, 50)
        assert ok and status == "monoexon_exact"

    def test_overlap(self):
        A = make_t([], start=100, end=500)
        B = make_t([], start=120, end=520)
        ok, status = can_merge_monoexon(A, B, 0.5, 500, 50)
        assert ok and status == "monoexon_overlap"

    def test_3end_too_far(self):
        A = make_t([], start=100, end=500)
        B = make_t([], start=100, end=700)  # 3' diff = 200
        ok, _ = can_merge_monoexon(A, B, 0.5, 500, 50)
        assert not ok


# ─────────────────────────────────────────────────────────────────
# Step2: similarity + cluster
# ─────────────────────────────────────────────────────────────────

class TestSimilarity:

    def test_identical_transcripts(self):
        A = make_t([(200,300),(500,600)], tid="t1")
        B = make_t([(200,300),(500,600)], tid="t2")
        s = compute_similarity(A, B, 3, 50, 0.4, 0.4, 0.2)
        assert s.total == pytest.approx(1.0)

    def test_completely_different(self):
        A = make_t([(200,300),(500,600)], chrom="chr1", tid="t1")
        # different chrom → not handled here; assume same chrom different positions
        B = make_t([(2000,2300),(2500,2600)], chrom="chr1", start=2000, end=3000, tid="t2")
        s = compute_similarity(A, B, 3, 50, 0.4, 0.4, 0.2)
        assert s.total < 0.5

    def test_drift_high_similarity(self):
        A = make_t([(200,300),(500,600),(800,900)], tid="t1")
        B = make_t([(201,300),(500,602),(799,900)], tid="t2")
        s = compute_similarity(A, B, 3, 50, 0.4, 0.4, 0.2)
        assert s.total > 0.9


# ─────────────────────────────────────────────────────────────────
# Class code
# ─────────────────────────────────────────────────────────────────

class TestClassCode:

    def _ref_index(self, ref_transcripts):
        idx = RefIndex()
        for t in ref_transcripts:
            idx.add_transcript(t)
        return idx

    def test_exact_match_equal(self):
        ref = make_t([(200,300),(500,600)], tid="ENST_001")
        q   = make_t([(200,300),(500,600)], tid="q1")
        code = classify_pair(q, ref, 3)
        assert code == ClassCode.EQUAL

    def test_junction_match(self):
        ref = make_t([(200,300),(500,600),(800,900)], tid="ENST_001")
        q   = make_t([(200,300),(700,900)], tid="q1")  # one matched, one new
        code = classify_pair(q, ref, 3)
        assert code == ClassCode.JUNCTION_MATCH

    def test_intronic(self):
        # ref has an intron from 300 to 500
        ref = make_t([(300,500)], start=100, end=800, tid="ENST_001")
        # query fits entirely inside that intron
        q   = make_t([], start=350, end=450, tid="q1")
        code = classify_pair(q, ref, 3)
        assert code == ClassCode.INTRONIC

    def test_opposite_strand_intron_match(self):
        ref = make_t([(200,300),(500,600)], strand="+", tid="ENST_001")
        q   = make_t([(200,300),(500,600)], strand="-", tid="q1")
        code = classify_pair(q, ref, 3)
        assert code == ClassCode.INTRON_MATCH_OPP

    def test_unknown_intergenic(self):
        ref = make_t([(200,300),(500,600)], start=100, end=800, tid="ENST_001")
        q   = make_t([(2200,2300)], start=2000, end=2500, tid="q1")
        code = classify_pair(q, ref, 3)
        assert code is None  # no overlap → caller assigns 'u'

    def test_classify_via_index(self):
        ref = make_t([(200,300),(500,600)], tid="ENST_001")
        q   = make_t([(200,300),(500,600)], tid="q1")
        idx = self._ref_index([ref])
        code, ref_t = classify_transcript(q, idx, 3)
        assert code == ClassCode.EQUAL
        assert ref_t.transcript_id == "ENST_001"

    def test_classify_no_overlap_returns_unknown(self):
        ref = make_t([(200,300),(500,600)], start=100, end=800, tid="ENST_001")
        q   = make_t([], start=5000, end=6000, tid="q1")
        idx = self._ref_index([ref])
        code, ref_t = classify_transcript(q, idx, 3)
        assert code == ClassCode.UNKNOWN
        assert ref_t is None


# ─────────────────────────────────────────────────────────────────
# Cluster
# ─────────────────────────────────────────────────────────────────

class TestCluster:

    def test_basic_clustering(self):
        # 3 similar, 1 different
        ts = [
            make_t([(200,300),(500,600)], tid="t1"),
            make_t([(201,300),(500,602)], tid="t2"),  # similar to t1
            make_t([(200,300),(500,600)], tid="t3"),  # identical to t1
            make_t([(2200,2300),(2500,2600)], start=2000, end=3000, tid="t4"),  # different
        ]
        clusters = cluster_transcripts(
            transcripts=ts,
            ref_ids=set(),
            sample_counts={},
            junction_threshold=3,
            fuzzy_threshold=50,
            w_structural=0.4,
            w_junction=0.4,
            w_jaccard=0.2,
            min_similarity=0.7,
        )
        # t1, t2, t3 → 1 cluster; t4 → another cluster
        assert len(clusters) == 2
        large = [c for c in clusters if c.size > 1]
        small = [c for c in clusters if c.size == 1]
        assert len(large) == 1
        assert large[0].size == 3
