"""
Step1: Merge - Merge transcripts likely representing the same molecule.

Strict merge rules (multi-exon):
  1. Same chrom + strand
  2. All matched junctions within --junction_threshold
  3. 5' missing junctions ≤ --max_5trunc (default 1)
  4. 3' missing junctions ≤ --max_3trunc (default 0)
  5. NO internal junction missing (retention/skipping forbidden)
  6. TSS difference ≤ --tss_window
  7. TES difference ≤ --tes_window
  8. Diff junction (non-zero error within threshold) count ≤ --max_diff_count
  9. Diff junction ratio ≤ --max_diff_ratio

Mono-exon rules:
  1. Same strand
  2. overlap_ratio ≥ --mono_overlap
  3. 5' diff ≤ --mono_5end
  4. 3' diff ≤ --mono_3end
"""
import logging
from collections import defaultdict, Counter
from typing import List, Dict, Tuple, Optional
from statistics import mean

from gtf_merge.core.models import (
    Transcript, MergedTranscript, Junction, Exon,
    MergeStatus, StrandStatus,
)

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────
# Junction matching utilities
# ─────────────────────────────────────────────────────────────────

def match_junctions(
    ja: List[Junction],
    jb: List[Junction],
    threshold: int,
) -> List[Tuple[int, int, int]]:
    """Greedy nearest-match. Returns [(idx_a, idx_b, coord_error), ...]."""
    matched = []
    used_b  = set()
    for i, a in enumerate(ja):
        best_j   = -1
        best_err = threshold + 1
        for j, b in enumerate(jb):
            if j in used_b:
                continue
            if not a.matches(b, threshold):
                continue
            err = a.coord_error(b)
            if err < best_err:
                best_err = err
                best_j   = j
        if best_j >= 0:
            matched.append((i, best_j, best_err))
            used_b.add(best_j)
    return matched


def is_continuous(indices: List[int]) -> bool:
    if len(indices) <= 1:
        return True
    return indices == list(range(indices[0], indices[-1] + 1))


# ─────────────────────────────────────────────────────────────────
# Mergeability check (multi-exon)
# ─────────────────────────────────────────────────────────────────

def can_merge_multiexon(
    a: Transcript,
    b: Transcript,
    junction_threshold: int,
    max_5trunc: int,
    max_3trunc: int,
    tss_window: int,
    tes_window: int,
    max_diff_count: int,
    max_diff_ratio: float,
) -> Tuple[bool, str]:
    """
    Decide whether two multi-exon transcripts can be merged.
    Returns (can_merge, status_string).
    """
    # 1. chrom + strand
    if a.chrom != b.chrom:
        return False, "diff_chrom"
    if a.strand != b.strand and a.strand != "." and b.strand != ".":
        return False, "diff_strand"

    # 2. find matched junctions
    matched = match_junctions(a.junctions, b.junctions, junction_threshold)
    if not matched:
        return False, "no_match"

    a_idx = [m[0] for m in matched]
    b_idx = [m[1] for m in matched]

    # 3. continuity check (no internal gap)
    if not is_continuous(a_idx) or not is_continuous(b_idx):
        return False, "internal_gap"  # retention/skipping

    # 4. truncation check
    a_miss_5 = a_idx[0]
    a_miss_3 = (len(a.junctions) - 1) - a_idx[-1]
    b_miss_5 = b_idx[0]
    b_miss_3 = (len(b.junctions) - 1) - b_idx[-1]

    # the longer transcript's missing = truncation in shorter
    miss_5 = max(a_miss_5, b_miss_5)
    miss_3 = max(a_miss_3, b_miss_3)

    if miss_5 > max_5trunc:
        return False, "5trunc_too_many"
    if max_3trunc >= 0 and miss_3 > max_3trunc:
        return False, "3trunc_too_many"

    # 5. TSS/TES window
    if abs(a.tss - b.tss) > tss_window:
        return False, "tss_too_far"
    if abs(a.tes - b.tes) > tes_window:
        return False, "tes_too_far"

    # 6. diff junction count/ratio (non-zero coord errors within threshold)
    diff_errors = [m[2] for m in matched if m[2] > 0]
    diff_count  = len(diff_errors)
    diff_ratio  = diff_count / len(matched) if matched else 0.0

    if max_diff_count >= 0 and diff_count > max_diff_count:
        return False, "diff_count_exceeded"
    if diff_ratio > max_diff_ratio:
        return False, "diff_ratio_exceeded"

    # determine status
    has_drift = diff_count > 0
    if miss_5 == 0 and miss_3 == 0:
        status = "fuzzy_match" if has_drift else "exact_match"
    elif miss_5 > 0 and miss_3 == 0:
        status = "truncation_5"
    elif miss_5 == 0 and miss_3 > 0:
        status = "truncation_3"
    else:
        status = "truncation_both"

    return True, status


# ─────────────────────────────────────────────────────────────────
# Mergeability check (mono-exon)
# ─────────────────────────────────────────────────────────────────

def can_merge_monoexon(
    a: Transcript,
    b: Transcript,
    mono_overlap: float,
    mono_5end:    int,
    mono_3end:    int,
) -> Tuple[bool, str]:
    if a.chrom != b.chrom:
        return False, "diff_chrom"
    if a.strand != b.strand and a.strand != "." and b.strand != ".":
        return False, "diff_strand"

    ea = a.exons[0] if a.exons else None
    eb = b.exons[0] if b.exons else None
    if ea is None or eb is None:
        return False, "no_exon"

    ov_start = max(ea.start, eb.start)
    ov_end   = min(ea.end,   eb.end)
    overlap  = max(0, ov_end - ov_start)
    shorter  = min(ea.length, eb.length)
    if shorter == 0:
        return False, "zero_length"

    overlap_ratio = overlap / shorter
    if overlap_ratio < mono_overlap:
        return False, "overlap_too_low"

    diff_5_genome = abs(ea.start - eb.start)
    diff_3_genome = abs(ea.end - eb.end)

    if a.strand == "+":
        diff_5, diff_3 = diff_5_genome, diff_3_genome
    elif a.strand == "-":
        diff_5, diff_3 = diff_3_genome, diff_5_genome
    else:
        diff_5 = max(diff_5_genome, diff_3_genome)
        diff_3 = min(diff_5_genome, diff_3_genome)

    if diff_5 > mono_5end:
        return False, "5end_too_far"
    if diff_3 > mono_3end:
        return False, "3end_too_far"

    if diff_5_genome == 0 and diff_3_genome == 0:
        return True, "monoexon_exact"
    return True, "monoexon_overlap"


# ─────────────────────────────────────────────────────────────────
# Representative selection helpers
# ─────────────────────────────────────────────────────────────────

def _consensus(coords: List[int]) -> int:
    return Counter(coords).most_common(1)[0][0]


def _pick_tss(group, strategy, strand):
    if strand == "+":
        vals = [(t.start, t.source_sample) for t in group]
        if strategy == "longest":
            return min(vals, key=lambda x: x[0])
        return (_consensus([v[0] for v in vals]), vals[0][1])
    else:
        vals = [(t.end, t.source_sample) for t in group]
        if strategy == "longest":
            return max(vals, key=lambda x: x[0])
        return (_consensus([v[0] for v in vals]), vals[0][1])


def _pick_tes(group, strategy, strand):
    if strand == "+":
        vals = [(t.end, t.source_sample) for t in group]
        if strategy == "longest":
            return max(vals, key=lambda x: x[0])
        return (_consensus([v[0] for v in vals]), vals[0][1])
    else:
        vals = [(t.start, t.source_sample) for t in group]
        if strategy == "longest":
            return min(vals, key=lambda x: x[0])
        return (_consensus([v[0] for v in vals]), vals[0][1])


def _build_representative_exons(group, tss, tes, strand):
    """Build representative exons using consensus coords + chosen TSS/TES."""
    scaffold = max(group, key=lambda t: len(t.exons))
    if not scaffold.exons:
        return []
    chrom   = scaffold.exons[0].chrom
    n       = len(scaffold.exons)
    rep_exons = []
    for ei, exon in enumerate(scaffold.exons):
        starts = [t.exons[ei].start for t in group if len(t.exons) > ei]
        ends   = [t.exons[ei].end   for t in group if len(t.exons) > ei]
        cs = _consensus(starts)
        ce = _consensus(ends)
        if strand == "+":
            if ei == 0:     cs = tss
            if ei == n - 1: ce = tes
        else:
            if ei == 0:     ce = tss
            if ei == n - 1: cs = tes
        if cs < ce:
            rep_exons.append(Exon(chrom=chrom, start=cs, end=ce, strand=strand))
    return rep_exons


def _build_merged_transcript(
    group:         List[Transcript],
    tss_strategy:  str,
    tes_strategy:  str,
    merge_status:  MergeStatus,
    total_samples: int,
    junction_threshold: int,
) -> MergedTranscript:
    """Construct a MergedTranscript from a group of source transcripts."""
    chrom  = group[0].chrom
    strand = group[0].strand

    # dedup by sample
    by_sample: Dict[str, Transcript] = {}
    for t in group:
        if t.source_sample not in by_sample:
            by_sample[t.source_sample] = t
    uniq = list(by_sample.values())

    tss_coord, tss_source = _pick_tss(uniq, tss_strategy, strand)
    tes_coord, tes_source = _pick_tes(uniq, tes_strategy, strand)
    g_start = min(tss_coord, tes_coord)
    g_end   = max(tss_coord, tes_coord)

    exons = _build_representative_exons(uniq, tss_coord, tes_coord, strand)
    junctions = [
        Junction(chrom=exons[i].chrom, donor=exons[i].end,
                 acceptor=exons[i+1].start, strand=strand)
        for i in range(len(exons) - 1)
    ]

    # junction support
    per_j_support: Dict[int, List[str]] = defaultdict(list)
    per_j_max_diff: Dict[int, int]      = {}
    for ji, jrep in enumerate(junctions):
        for t in uniq:
            if ji < len(t.junctions):
                err = jrep.coord_error(t.junctions[ji])
                per_j_support[ji].append(t.source_sample)
                per_j_max_diff[ji] = max(per_j_max_diff.get(ji, 0), err)

    # diff junction count
    diff_count = sum(1 for v in per_j_max_diff.values() if v > 0)
    diff_ratio = diff_count / len(junctions) if junctions else 0.0

    # TSS/TES diff ranges
    all_tss = [t.tss for t in uniq]
    all_tes = [t.tes for t in uniq]
    tss_diff = f"{min(abs(v - tss_coord) for v in all_tss)}-{max(abs(v - tss_coord) for v in all_tss)}bp"
    tes_diff = f"{min(abs(v - tes_coord) for v in all_tes)}-{max(abs(v - tes_coord) for v in all_tes)}bp"

    # truncation tracking
    rep_n_juncs = len(junctions)
    trunc5_count, trunc3_count = 0, 0
    miss5_parts, miss3_parts = [], []
    for t in uniq:
        if t.num_junctions < rep_n_juncs:
            # find which junctions are missing
            matched = match_junctions(junctions, t.junctions, junction_threshold)
            if matched:
                a_idx = sorted(set(m[0] for m in matched))
                m5 = a_idx[0]
                m3 = (rep_n_juncs - 1) - a_idx[-1]
                if m5 > 0:
                    trunc5_count += 1
                    miss5_parts.append(f"{t.source_sample}:first_{m5}")
                if m3 > 0:
                    trunc3_count += 1
                    miss3_parts.append(f"{t.source_sample}:last_{m3}")

    freq_score = len(uniq) / max(total_samples, 1)

    return MergedTranscript(
        transcript_id = "PENDING",  # will be assigned later
        gene_id       = "PENDING",
        chrom         = chrom,
        start         = g_start,
        end           = g_end,
        strand        = strand,
        exons         = exons,
        junctions     = junctions,
        source_transcripts = [f"{t.source_sample}:{t.transcript_id}" for t in uniq],
        sample_ids    = [t.source_sample for t in uniq],
        num_samples   = len(uniq),
        merge_status  = merge_status,
        per_sample_status = {t.source_sample: merge_status.value for t in uniq},
        tss_coord     = tss_coord,
        tes_coord     = tes_coord,
        tss_source    = tss_source,
        tes_source    = tes_source,
        tss_diff_range= tss_diff,
        tes_diff_range= tes_diff,
        truncation_5_count  = trunc5_count,
        truncation_3_count  = trunc3_count,
        missing_junctions_5 = ";".join(miss5_parts) if miss5_parts else "NA",
        missing_junctions_3 = ";".join(miss3_parts) if miss3_parts else "NA",
        per_junction_support  = dict(per_j_support),
        per_junction_max_diff = per_j_max_diff,
        diff_junction_count   = diff_count,
        diff_junction_ratio   = round(diff_ratio, 4),
        freq_score    = round(freq_score, 4),
        strand_status = StrandStatus.UNKNOWN if strand == "." else StrandStatus.KNOWN,
    )


# ─────────────────────────────────────────────────────────────────
# Multi-exon merge driver (cd-hit style: incremental clustering)
# ─────────────────────────────────────────────────────────────────

def merge_multiexon_transcripts(
    transcripts: List[Transcript],
    junction_threshold: int,
    max_5trunc:    int,
    max_3trunc:    int,
    tss_window:    int,
    tes_window:    int,
    max_diff_count: int,
    max_diff_ratio: float,
    tss_strategy: str,
    tes_strategy: str,
    total_samples: int,
) -> List[MergedTranscript]:
    """
    cd-hit style incremental merging.
    Sort by (num_samples desc, length desc) so high-confidence ones become representatives first.
    """
    if not transcripts:
        return []

    # group by chrom + strand
    cs_groups: Dict[Tuple, List[Transcript]] = defaultdict(list)
    for t in transcripts:
        cs_groups[(t.chrom, t.strand)].append(t)

    merged_all: List[MergedTranscript] = []

    for (chrom, strand), group in cs_groups.items():
        # sort: more junctions first (longer scaffold), then length
        group.sort(key=lambda t: (-t.num_junctions, -t.length))

        # incremental clustering
        clusters: List[List[Transcript]] = []
        cluster_status: List[str] = []

        for t in group:
            placed = False
            for ci, cluster in enumerate(clusters):
                rep = cluster[0]
                ok, status = can_merge_multiexon(
                    rep, t,
                    junction_threshold = junction_threshold,
                    max_5trunc = max_5trunc,
                    max_3trunc = max_3trunc,
                    tss_window = tss_window,
                    tes_window = tes_window,
                    max_diff_count = max_diff_count,
                    max_diff_ratio = max_diff_ratio,
                )
                if ok:
                    cluster.append(t)
                    placed = True
                    # upgrade status if needed
                    if status != cluster_status[ci]:
                        # priority: exact > fuzzy > truncation
                        priority = {"exact_match": 1, "fuzzy_match": 2,
                                    "truncation_5": 3, "truncation_3": 3,
                                    "truncation_both": 4}
                        if priority.get(status, 5) < priority.get(cluster_status[ci], 5):
                            cluster_status[ci] = status
                    break
            if not placed:
                clusters.append([t])
                cluster_status.append("singleton")

        # build MergedTranscript objects
        for cluster, status in zip(clusters, cluster_status):
            # determine final status
            if len(cluster) == 1:
                final_status = MergeStatus.SINGLETON
            else:
                final_status = MergeStatus(status if status in {s.value for s in MergeStatus}
                                           else "exact_match")

            mt = _build_merged_transcript(
                group         = cluster,
                tss_strategy  = tss_strategy,
                tes_strategy  = tes_strategy,
                merge_status  = final_status,
                total_samples = total_samples,
                junction_threshold = junction_threshold,
            )
            merged_all.append(mt)

    return merged_all


# ─────────────────────────────────────────────────────────────────
# Mono-exon merge driver
# ─────────────────────────────────────────────────────────────────

def merge_monoexon_transcripts(
    transcripts: List[Transcript],
    mono_overlap: float,
    mono_5end:    int,
    mono_3end:    int,
    tss_strategy: str,
    tes_strategy: str,
    total_samples: int,
) -> List[MergedTranscript]:
    if not transcripts:
        return []

    cs_groups: Dict[Tuple, List[Transcript]] = defaultdict(list)
    for t in transcripts:
        cs_groups[(t.chrom, t.strand)].append(t)

    merged_all: List[MergedTranscript] = []
    for (chrom, strand), group in cs_groups.items():
        group.sort(key=lambda t: t.start)

        clusters: List[List[Transcript]] = []
        cluster_status: List[str] = []

        for t in group:
            placed = False
            for ci, cluster in enumerate(clusters):
                rep = cluster[0]
                ok, status = can_merge_monoexon(
                    rep, t, mono_overlap, mono_5end, mono_3end,
                )
                if ok:
                    cluster.append(t)
                    if status == "monoexon_exact" and cluster_status[ci] != "monoexon_exact":
                        pass  # keep more strict
                    elif status != cluster_status[ci]:
                        cluster_status[ci] = "monoexon_overlap"
                    placed = True
                    break
            if not placed:
                clusters.append([t])
                cluster_status.append("monoexon_exact")

        for cluster, status in zip(clusters, cluster_status):
            if len(cluster) == 1:
                final_status = MergeStatus.SINGLETON
            else:
                final_status = MergeStatus(status)

            mt = _build_merged_transcript(
                group         = cluster,
                tss_strategy  = tss_strategy,
                tes_strategy  = tes_strategy,
                merge_status  = final_status,
                total_samples = total_samples,
                junction_threshold = 0,
            )
            merged_all.append(mt)

    return merged_all


# ─────────────────────────────────────────────────────────────────
# Top-level driver
# ─────────────────────────────────────────────────────────────────

def run_merge(
    transcripts: List[Transcript],
    cfg,
    total_samples: int,
) -> List[MergedTranscript]:
    """Top-level merge driver. Splits multi-exon vs mono-exon and merges each."""
    multi = [t for t in transcripts if not t.is_monoexon]
    mono  = [t for t in transcripts if t.is_monoexon]

    logger.debug(f"Merging {len(multi)} multi-exon, {len(mono)} mono-exon")

    merged_multi = merge_multiexon_transcripts(
        transcripts        = multi,
        junction_threshold = cfg.junction_threshold,
        max_5trunc         = cfg.max_5trunc,
        max_3trunc         = cfg.max_3trunc,
        tss_window         = cfg.tss_window,
        tes_window         = cfg.tes_window,
        max_diff_count     = cfg.max_diff_count,
        max_diff_ratio     = cfg.max_diff_ratio,
        tss_strategy       = cfg.tss_strategy,
        tes_strategy       = cfg.tes_strategy,
        total_samples      = total_samples,
    )

    merged_mono = merge_monoexon_transcripts(
        transcripts   = mono,
        mono_overlap  = cfg.mono_overlap,
        mono_5end     = cfg.mono_5end,
        mono_3end     = cfg.mono_3end,
        tss_strategy  = cfg.tss_strategy,
        tes_strategy  = cfg.tes_strategy,
        total_samples = total_samples,
    )

    return merged_multi + merged_mono
