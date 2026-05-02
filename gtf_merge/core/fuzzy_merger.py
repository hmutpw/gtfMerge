"""
Step 2: Fuzzy merge using intron-based scoring + graph connected components.

Key changes from v1:
  - Uses compute_merge_score() (intron-based) instead of junction similarity
  - Score threshold replaces similarity threshold
  - Handles 5' truncation natively (score=100 by default)
  - Retention/skipping auto-detected and heavily penalised
"""
import logging
from collections import defaultdict
from typing import List, Dict, Tuple
from statistics import mean

from gtf_merge.core.models import (
    Transcript, MergedTranscript, Exon, Intron, Junction,
    MergeType, StrandStatus, MergeScore
)
from gtf_merge.core.similarity import compute_merge_score, match_introns
from gtf_merge.core.exact_merger import _build_merged, _pick_tss, _pick_tes, _consensus_coord

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Union-Find
# ---------------------------------------------------------------------------

class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank   = [0] * n

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx == ry: return
        if self.rank[rx] < self.rank[ry]: rx, ry = ry, rx
        self.parent[ry] = rx
        if self.rank[rx] == self.rank[ry]: self.rank[rx] += 1

    def components(self):
        groups = defaultdict(list)
        for i in range(len(self.parent)):
            groups[self.find(i)].append(i)
        return dict(groups)


# ---------------------------------------------------------------------------
# Cluster merger
# ---------------------------------------------------------------------------

def _merge_cluster(
    group:          List[Transcript],
    scores:         Dict[Tuple[int,int], MergeScore],  # (i,j) → score
    threshold:      int,
    tss_strategy:   str,
    tes_strategy:   str,
    merged_id:      str,
    gene_id:        str,
    total_samples:  int,
    trunc5_penalty: float,
    trunc3_penalty: float,
) -> MergedTranscript:
    """Merge one fuzzy cluster into a representative MergedTranscript."""

    # pick representative: most introns → most samples
    representative = max(group, key=lambda t: (len(t.introns), 1))

    # deduplicate by sample
    by_sample: Dict[str, Transcript] = {}
    for t in group:
        if t.source_sample not in by_sample:
            by_sample[t.source_sample] = t
    uniq = list(by_sample.values())

    chrom  = representative.chrom
    strand = representative.strand

    # per-sample scores
    per_score:  Dict[str, float]     = {}
    per_mtype:  Dict[str, str]       = {}
    per_src:    Dict[str, str]       = {}
    trunc5_count = 0
    trunc3_count = 0
    miss5_parts, miss3_parts = [], []

    for t in uniq:
        ms = compute_merge_score(representative, t, threshold,
                                  trunc5_penalty=trunc5_penalty,
                                  trunc3_penalty=trunc3_penalty)
        per_score[t.source_sample] = ms.score
        per_mtype[t.source_sample] = ms.merge_type.value
        per_src[t.source_sample]   = t.transcript_id
        if ms.missing_5 > 0:
            trunc5_count += 1
            miss5_parts.append(f"{t.source_sample}:first_{ms.missing_5}")
        if ms.missing_3 > 0:
            trunc3_count += 1
            miss3_parts.append(f"{t.source_sample}:last_{ms.missing_3}")

    avg_score = mean(per_score.values()) if per_score else 100.0

    # dominant merge type
    type_votes = defaultdict(int)
    for mt in per_mtype.values(): type_votes[mt] += 1
    dominant_type = MergeType(max(type_votes, key=type_votes.get))

    # build representative structure
    mt_obj = _build_merged(
        group         = uniq,
        merged_id     = merged_id,
        gene_id       = gene_id,
        tss_strategy  = tss_strategy,
        tes_strategy  = tes_strategy,
        total_samples = total_samples,
        merge_type    = dominant_type,
        avg_score     = avg_score,
    )

    # override per-sample fields with detailed scores
    mt_obj.per_sample_merge_type  = per_mtype
    mt_obj.per_sample_score       = {k: round(v,2) for k,v in per_score.items()}
    mt_obj.per_sample_source_id   = per_src
    mt_obj.truncation_5_count     = trunc5_count
    mt_obj.truncation_3_count     = trunc3_count
    mt_obj.missing_introns_5      = ";".join(miss5_parts) if miss5_parts else "NA"
    mt_obj.missing_introns_3      = ";".join(miss3_parts) if miss3_parts else "NA"

    return mt_obj


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_fuzzy_merge(
    transcripts:    List[Transcript],
    threshold:      int,
    min_score:      float,
    tss_strategy:   str,
    tes_strategy:   str,
    id_prefix:      str,
    gene_counter:   Dict,
    total_samples:  int,
    trunc5_penalty: float,
    trunc3_penalty: float,
    retention_penalty: float,
    min_shared_introns: int,
) -> List[MergedTranscript]:
    """
    Fuzzy merge transcripts using intron-based scores.
    Graph edge: score >= min_score AND is_mergeable.
    Connected components form merge clusters.
    Singletons returned as SINGLETON.
    """
    if not transcripts:
        return []

    n  = len(transcripts)
    uf = UnionFind(n)

    # group by chrom+strand for efficiency
    cs_groups: Dict[Tuple, List[int]] = defaultdict(list)
    for i, t in enumerate(transcripts):
        cs_groups[(t.chrom, t.strand)].append(i)

    # score cache
    score_cache: Dict[Tuple[int,int], MergeScore] = {}

    for (chrom, strand), idxs in cs_groups.items():
        for ii in range(len(idxs)):
            for jj in range(ii+1, len(idxs)):
                i, j = idxs[ii], idxs[jj]
                ta, tb = transcripts[i], transcripts[j]

                # quick genomic overlap check before expensive scoring
                if ta.end < tb.start or tb.end < ta.start:
                    continue

                ms = compute_merge_score(
                    ta, tb, threshold,
                    trunc5_penalty     = trunc5_penalty,
                    trunc3_penalty     = trunc3_penalty,
                    retention_penalty  = retention_penalty,
                    min_shared_introns = min_shared_introns,
                )
                score_cache[(i,j)] = ms

                if ms.is_mergeable and ms.score >= min_score:
                    uf.union(i, j)

    components = uf.components()
    merged: List[MergedTranscript] = []

    for root, members in components.items():
        group  = [transcripts[i] for i in members]
        chrom  = group[0].chrom
        strand = group[0].strand
        key    = (chrom, strand)

        if key not in gene_counter:
            gene_counter[key] = {"g": 0, "t": defaultdict(int)}
        gene_counter[key]["g"] += 1
        gn = gene_counter[key]["g"]
        gene_counter[key]["t"][gn] += 1
        tn = gene_counter[key]["t"][gn]
        mid = f"{id_prefix}G{gn}.{tn}"
        gid = f"{id_prefix}G{gn}"

        if len(group) == 1:
            # singleton
            t  = group[0]
            mt = MergedTranscript(
                merged_id    = mid,
                gene_id      = gid,
                chrom        = t.chrom,
                start        = t.start,
                end          = t.end,
                strand       = t.strand,
                exons        = t.exons,
                introns      = t.introns,
                junctions    = t.junctions,
                merge_type   = MergeType.SINGLETON,
                merge_score  = 100.0,
                num_samples  = 1,
                sample_ids   = [t.source_sample],
                source_transcripts   = [f"{t.source_sample}:{t.transcript_id}"],
                per_sample_merge_type= {t.source_sample: MergeType.SINGLETON.value},
                per_sample_score     = {t.source_sample: 100.0},
                per_sample_source_id = {t.source_sample: t.transcript_id},
                tss_coord    = t.tss,
                tes_coord    = t.tes,
                tss_source   = t.source_sample,
                tes_source   = t.source_sample,
                freq_score   = 1 / max(total_samples, 1),
                strand_status= StrandStatus.UNKNOWN if t.strand == "." else StrandStatus.KNOWN,
            )
            merged.append(mt)
        else:
            mt = _merge_cluster(
                group           = group,
                scores          = score_cache,
                threshold       = threshold,
                tss_strategy    = tss_strategy,
                tes_strategy    = tes_strategy,
                merged_id       = mid,
                gene_id         = gid,
                total_samples   = total_samples,
                trunc5_penalty  = trunc5_penalty,
                trunc3_penalty  = trunc3_penalty,
            )
            merged.append(mt)

    logger.debug(f"Fuzzy merge: {n} transcripts → {len(merged)} merged/singleton")
    return merged
