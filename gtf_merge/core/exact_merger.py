"""Step 1: Exact merge - group transcripts with identical intron chains."""
import logging
from collections import defaultdict
from typing import List, Dict, Tuple
from statistics import mean

from gtf_merge.core.models import (
    Transcript, MergedTranscript, Exon, Intron, Junction,
    MergeType, StrandStatus
)

logger = logging.getLogger(__name__)


def _mode(vals: List[int]) -> int:
    from collections import Counter
    return Counter(vals).most_common(1)[0][0]


def _pick_tss(transcripts, strategy, strand):
    if strand == "+":
        vals = [(t.start, t.source_sample) for t in transcripts]
        return min(vals, key=lambda x: x[0]) if strategy == "longest" \
               else (_mode([v[0] for v in vals]), vals[0][1])
    else:
        vals = [(t.end, t.source_sample) for t in transcripts]
        return max(vals, key=lambda x: x[0]) if strategy == "longest" \
               else (_mode([v[0] for v in vals]), vals[0][1])


def _pick_tes(transcripts, strategy, strand):
    if strand == "+":
        vals = [(t.end, t.source_sample) for t in transcripts]
        return max(vals, key=lambda x: x[0]) if strategy == "longest" \
               else (_mode([v[0] for v in vals]), vals[0][1])
    else:
        vals = [(t.start, t.source_sample) for t in transcripts]
        return min(vals, key=lambda x: x[0]) if strategy == "longest" \
               else (_mode([v[0] for v in vals]), vals[0][1])


def _consensus_coord(coords: List[int]) -> int:
    from collections import Counter
    return Counter(coords).most_common(1)[0][0]


def _build_rep_exons(transcripts, tss, tes, strand):
    """Build representative exons using consensus internal coords + chosen TSS/TES."""
    scaffold = max(transcripts, key=lambda t: len(t.exons))
    if not scaffold.exons:
        return []
    chrom   = scaffold.exons[0].chrom
    n_exons = len(scaffold.exons)
    rep_exons = []
    for ei, exon in enumerate(scaffold.exons):
        starts = [t.exons[ei].start for t in transcripts if len(t.exons) > ei]
        ends   = [t.exons[ei].end   for t in transcripts if len(t.exons) > ei]
        cs = _consensus_coord(starts)
        ce = _consensus_coord(ends)
        if strand == "+":
            if ei == 0:         cs = tss
            if ei == n_exons-1: ce = tes
        else:
            if ei == 0:         ce = tss
            if ei == n_exons-1: cs = tes
        if cs < ce:
            rep_exons.append(Exon(chrom=chrom, start=cs, end=ce, strand=strand))
    return rep_exons


def _build_merged(
    group: List[Transcript],
    merged_id: str,
    gene_id: str,
    tss_strategy: str,
    tes_strategy: str,
    total_samples: int,
    merge_type: MergeType = MergeType.EXACT_MATCH,
    avg_score: float = 100.0,
) -> MergedTranscript:
    """Build a MergedTranscript from a group."""
    chrom  = group[0].chrom
    strand = group[0].strand

    # deduplicate by sample
    by_sample: Dict[str, Transcript] = {}
    for t in group:
        if t.source_sample not in by_sample:
            by_sample[t.source_sample] = t
    uniq = list(by_sample.values())

    tss_coord, tss_source = _pick_tss(uniq, tss_strategy, strand)
    tes_coord, tes_source = _pick_tes(uniq, tes_strategy, strand)
    g_start = min(tss_coord, tes_coord)
    g_end   = max(tss_coord, tes_coord)

    exons = _build_rep_exons(uniq, tss_coord, tes_coord, strand)

    # rebuild introns from representative exons
    introns = [
        Intron(chrom=exons[i].chrom, donor=exons[i].end,
               acceptor=exons[i+1].start, strand=strand)
        for i in range(len(exons)-1)
    ]
    junctions = [Junction(i.chrom,i.donor,i.acceptor,i.strand) for i in introns]

    # intron support stats
    per_i_support: Dict[int, List[str]] = defaultdict(list)
    per_i_max_diff: Dict[int, int] = {}
    anchor_introns = []
    splice_sites   = []

    for ii, irep in enumerate(introns):
        for t in uniq:
            if ii < len(t.introns):
                err = irep.coord_error(t.introns[ii])
                per_i_support[ii].append(t.source_sample)
                per_i_max_diff[ii] = max(per_i_max_diff.get(ii, 0), err)
        anchor_introns.append(per_i_max_diff.get(ii, 0) == 0)
        splice_sites.append(irep.splice_site)

    all_tss = [t.tss for t in uniq]
    all_tes = [t.tes for t in uniq]
    tss_diff = f"{min(abs(v-tss_coord) for v in all_tss)}-{max(abs(v-tss_coord) for v in all_tss)}bp"
    tes_diff = f"{min(abs(v-tes_coord) for v in all_tes)}-{max(abs(v-tes_coord) for v in all_tes)}bp"

    freq_score = len(uniq) / max(total_samples, 1)

    return MergedTranscript(
        merged_id    = merged_id,
        gene_id      = gene_id,
        chrom        = chrom,
        start        = g_start,
        end          = g_end,
        strand       = strand,
        exons        = exons,
        introns      = introns,
        junctions    = junctions,
        merge_type   = merge_type,
        merge_score  = round(avg_score, 2),
        num_samples  = len(uniq),
        sample_ids   = [t.source_sample for t in uniq],
        source_transcripts = [f"{t.source_sample}:{t.transcript_id}" for t in uniq],
        per_sample_merge_type = {t.source_sample: merge_type.value for t in uniq},
        per_sample_score      = {t.source_sample: avg_score for t in uniq},
        per_sample_source_id  = {t.source_sample: t.transcript_id for t in uniq},
        tss_coord    = tss_coord,
        tes_coord    = tes_coord,
        tss_source   = tss_source,
        tes_source   = tes_source,
        tss_diff_range = tss_diff,
        tes_diff_range = tes_diff,
        per_intron_sample_support = dict(per_i_support),
        per_intron_max_diff       = per_i_max_diff,
        anchor_introns            = anchor_introns,
        splice_sites              = splice_sites,
        freq_score   = round(freq_score, 4),
        strand_status = StrandStatus.UNKNOWN if strand == "." else StrandStatus.KNOWN,
    )


def run_exact_merge(
    transcripts:   List[Transcript],
    tss_strategy:  str,
    tes_strategy:  str,
    id_prefix:     str,
    gene_counter:  Dict,
    total_samples: int,
) -> List[MergedTranscript]:
    """
    Group multi-exon transcripts by identical intron_chain.
    Returns merged list (singletons included as SINGLETON type).
    """
    chain_groups: Dict[str, List[Transcript]] = defaultdict(list)
    for t in transcripts:
        if not t.is_monoexon:
            chain_groups[t.intron_chain].append(t)

    merged = []
    for chain, group in chain_groups.items():
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
        mt = _build_merged(
            group         = group,
            merged_id     = mid,
            gene_id       = gid,
            tss_strategy  = tss_strategy,
            tes_strategy  = tes_strategy,
            total_samples = total_samples,
            merge_type    = MergeType.EXACT_MATCH,
            avg_score     = 100.0,
        )
        merged.append(mt)

    logger.debug(f"Exact merge: {len(chain_groups)} intron chains → {len(merged)} merged")
    return merged
