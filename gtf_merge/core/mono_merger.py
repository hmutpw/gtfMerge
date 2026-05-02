"""Step 3: Mono-exon transcript merging."""
import logging
from collections import defaultdict
from typing import List, Dict, Tuple

from gtf_merge.core.models import (
    Transcript, MergedTranscript, Exon, Junction, MergeType, StrandStatus
)
from gtf_merge.core.similarity import monoexon_similarity

logger = logging.getLogger(__name__)


def _merge_mono_cluster(
    group:         List[Transcript],
    merged_id:     str,
    gene_id:       str,
    total_samples: int,
    tss_strategy:  str,
    tes_strategy:  str,
    merge_type_str: str,
) -> MergedTranscript:
    """Build a MergedTranscript from a group of mono-exon transcripts."""
    chrom  = group[0].chrom
    strand = group[0].strand

    # deduplicate by sample
    by_sample: Dict[str, Transcript] = {}
    for t in group:
        if t.source_sample not in by_sample:
            by_sample[t.source_sample] = t
    uniq = list(by_sample.values())

    # representative exon: largest span
    g_start = min(t.start for t in uniq)
    g_end   = max(t.end   for t in uniq)
    rep_exon = Exon(chrom=chrom, start=g_start, end=g_end, strand=strand)

    # TSS/TES
    if strand == "+":
        tss_coord = min(t.start for t in uniq) if tss_strategy == "longest" else _mode([t.start for t in uniq])
        tes_coord = max(t.end   for t in uniq) if tes_strategy == "longest" else _mode([t.end   for t in uniq])
        tss_source = next(t.source_sample for t in uniq if t.start == tss_coord)
        tes_source = next(t.source_sample for t in uniq if t.end   == tes_coord)
    elif strand == "-":
        tss_coord  = max(t.end   for t in uniq) if tss_strategy == "longest" else _mode([t.end   for t in uniq])
        tes_coord  = min(t.start for t in uniq) if tes_strategy == "longest" else _mode([t.start for t in uniq])
        tss_source = next(t.source_sample for t in uniq if t.end   == tss_coord)
        tes_source = next(t.source_sample for t in uniq if t.start == tes_coord)
    else:
        tss_coord  = min(t.start for t in uniq)
        tes_coord  = max(t.end   for t in uniq)
        tss_source = uniq[0].source_sample
        tes_source = uniq[0].source_sample

    all_tss = [t.tss for t in uniq]
    all_tes = [t.tes for t in uniq]
    tss_diff = f"{min(abs(v-tss_coord) for v in all_tss)}-{max(abs(v-tss_coord) for v in all_tss)}bp"
    tes_diff = f"{min(abs(v-tes_coord) for v in all_tes)}-{max(abs(v-tes_coord) for v in all_tes)}bp"

    merge_type = MergeType(merge_type_str) if merge_type_str else MergeType.MONOEXON_EXACT
    freq_score = len(uniq) / max(total_samples, 1)

    return MergedTranscript(
        merged_id    = merged_id,
        gene_id      = gene_id,
        chrom        = chrom,
        start        = g_start,
        end          = g_end,
        strand       = strand,
        exons        = [rep_exon],
        junctions    = [],
        merge_type   = merge_type,
        merge_score  = 100.0 if merge_type == MergeType.MONOEXON_EXACT else 80.0,
        num_samples  = len(uniq),
        sample_ids   = [t.source_sample for t in uniq],
        source_transcripts = [f"{t.source_sample}:{t.transcript_id}" for t in uniq],
        per_sample_merge_type = {t.source_sample: merge_type.value for t in uniq},
        per_sample_score = {t.source_sample: 100.0 for t in uniq},
        per_sample_source_id  = {t.source_sample: t.transcript_id for t in uniq},
        tss_coord    = tss_coord,
        tes_coord    = tes_coord,
        tss_source   = tss_source,
        tes_source   = tes_source,
        tss_diff_range = tss_diff,
        tes_diff_range = tes_diff,
        freq_score   = round(freq_score, 4),
        strand_status = StrandStatus.UNKNOWN if strand == "." else StrandStatus.KNOWN,
    )


def _mode(vals: List[int]) -> int:
    from collections import Counter
    return Counter(vals).most_common(1)[0][0]


def run_monoexon_merge(
    transcripts:    List[Transcript],
    mono_overlap:   float,
    mono_3end:      int,
    mono_5end:      int,
    tss_strategy:   str,
    tes_strategy:   str,
    id_prefix:      str,
    gene_counter:   Dict,
    total_samples:  int,
) -> List[MergedTranscript]:
    """
    Merge mono-exon transcripts using overlap + 3'/5' end thresholds.
    Uses greedy clustering: sort by start, extend clusters greedily.
    """
    if not transcripts:
        return []

    # separate by chrom+strand
    by_chrom_strand: Dict[Tuple, List[Transcript]] = defaultdict(list)
    for t in transcripts:
        by_chrom_strand[(t.chrom, t.strand)].append(t)

    merged_all: List[MergedTranscript] = []

    for (chrom, strand), group in by_chrom_strand.items():
        # sort by start coordinate
        group.sort(key=lambda t: t.start)

        clusters: List[List[Transcript]] = []
        cluster_types: List[str] = []

        for t in group:
            placed = False
            for ci, cluster in enumerate(clusters):
                # try to merge with existing cluster
                # use representative (first transcript added) for comparison
                rep = cluster[0]
                should_merge, mtype = monoexon_similarity(
                    rep, t, mono_overlap, mono_3end, mono_5end
                )
                if should_merge:
                    cluster.append(t)
                    # upgrade merge type if needed
                    if cluster_types[ci] == MergeType.MONOEXON_EXACT and mtype != MergeType.MONOEXON_EXACT:
                        cluster_types[ci] = mtype
                    placed = True
                    break
            if not placed:
                clusters.append([t])
                cluster_types.append(MergeType.MONOEXON_EXACT)

        for cluster, ctype in zip(clusters, cluster_types):
            key = (chrom, strand)
            if key not in gene_counter:
                gene_counter[key] = {"g": 0, "t": defaultdict(int)}
            gene_counter[key]["g"] += 1
            gn = gene_counter[key]["g"]
            gene_counter[key]["t"][gn] += 1
            tn = gene_counter[key]["t"][gn]

            mid = f"{id_prefix}G{gn}.{tn}"
            gid = f"{id_prefix}G{gn}"

            mt = _merge_mono_cluster(
                group          = cluster,
                merged_id      = mid,
                gene_id        = gid,
                total_samples  = total_samples,
                tss_strategy   = tss_strategy,
                tes_strategy   = tes_strategy,
                merge_type_str = ctype if isinstance(ctype, str) else ctype.value,
            )
            merged_all.append(mt)

    logger.debug(f"Monoexon merge: {len(transcripts)} → {len(merged_all)} merged")
    return merged_all
