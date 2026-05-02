"""Define gene boundaries from merged transcripts using exon overlap."""
import logging
from collections import defaultdict
from typing import List, Dict, Tuple

from gtf_merge.core.models import MergedTranscript, Gene

logger = logging.getLogger(__name__)


def _exons_overlap(mt_a: MergedTranscript, mt_b: MergedTranscript) -> bool:
    """Check if any exon in mt_a overlaps any exon in mt_b."""
    if mt_a.chrom != mt_b.chrom:
        return False
    if mt_a.strand != mt_b.strand and mt_a.strand != "." and mt_b.strand != ".":
        return False
    for ea in mt_a.exons:
        for eb in mt_b.exons:
            if ea.start < eb.end and ea.end > eb.start:
                return True
    return False


def define_genes(
    merged:       List[MergedTranscript],
    id_prefix:    str,
    ref_gene_map: Dict[str, str] = None,   # transcript_id -> ref_gene_id
) -> List[Gene]:
    """
    Group merged transcripts into genes based on exon overlap.
    Transcripts from the same chromosome and strand that share exon overlap
    are placed in the same gene.

    If a transcript has a ref_gene_id from GENCODE comparison, all transcripts
    with the same ref_gene_id are also forced into the same gene.
    """
    if not merged:
        return []

    # separate by chrom+strand
    by_cs: Dict[Tuple, List[MergedTranscript]] = defaultdict(list)
    for mt in merged:
        by_cs[(mt.chrom, mt.strand)].append(mt)

    all_genes: List[Gene] = []
    gene_counter = 0

    for (chrom, strand), group in by_cs.items():
        # sort by start coordinate
        group.sort(key=lambda x: x.start)
        n = len(group)

        # Union-Find for grouping
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x, y):
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[ry] = rx

        # 1. merge by exon overlap (sweep line: sorted by start)
        for i in range(n):
            for j in range(i + 1, n):
                if group[j].start > group[i].end:
                    break   # no more overlap possible for group[i]
                if _exons_overlap(group[i], group[j]):
                    union(i, j)

        # 2. merge by shared ref_gene_id
        if ref_gene_map:
            ref_groups: Dict[str, List[int]] = defaultdict(list)
            for i, mt in enumerate(group):
                if mt.ref_gene_id != "NA":
                    ref_groups[mt.ref_gene_id].append(i)
            for members in ref_groups.values():
                for m in members[1:]:
                    union(members[0], m)

        # collect components
        components: Dict[int, List[int]] = defaultdict(list)
        for i in range(n):
            components[find(i)].append(i)

        for root, members in components.items():
            gene_counter += 1
            gene_transcripts = [group[i] for i in members]

            g_start = min(mt.start for mt in gene_transcripts)
            g_end   = max(mt.end   for mt in gene_transcripts)

            # determine gene_id
            # prefer ref_gene_id if all transcripts share one
            ref_ids = set(mt.ref_gene_id for mt in gene_transcripts if mt.ref_gene_id != "NA")
            if len(ref_ids) == 1:
                gene_id = next(iter(ref_ids))
                ref_match = "known"
            elif len(ref_ids) > 1:
                gene_id   = f"{id_prefix}G{gene_counter}"
                ref_match = "multi_ref"
            else:
                gene_id   = f"{id_prefix}G{gene_counter}"
                ref_match = "novel"

            # update gene_id on all member transcripts
            for mt in gene_transcripts:
                mt.gene_id = gene_id

            # gene name from ref
            ref_gene_names = set(
                mt.attributes.get("gene_name", "NA")
                if hasattr(mt, "attributes") else "NA"
                for mt in gene_transcripts
            )
            ref_gene_name = next(iter(ref_gene_names - {"NA"}), "NA")

            gene = Gene(
                gene_id       = gene_id,
                chrom         = chrom,
                start         = g_start,
                end           = g_end,
                strand        = strand,
                transcripts   = gene_transcripts,
                ref_gene_id   = next(iter(ref_ids), "NA"),
                ref_gene_name = ref_gene_name,
                
            )
            all_genes.append(gene)

    # sort by chrom, start
    all_genes.sort(key=lambda g: (g.chrom, g.start))
    logger.debug(f"Gene definition: {len(merged)} transcripts → {len(all_genes)} genes")
    return all_genes
