"""Output writers: GTF and tracking TSV."""
import logging
from typing import List, TextIO
from gtf_merge.core.models import MergedTranscript, Gene

logger = logging.getLogger(__name__)


def _attr(key: str, val) -> str:
    if isinstance(val, float):
        return f'{key} "{val:.4f}"'
    return f'{key} "{val}"'


def write_merged_gtf(
    genes: List[Gene],
    out_fh: TextIO,
    source_tag: str = "GTFmerge",
):
    for gene in genes:
        gene_attrs = "; ".join([
            _attr("gene_id",         gene.gene_id),
            _attr("gene_name",       gene.ref_gene_name),
            _attr("ref_gene_id",     gene.ref_gene_id),
            _attr("num_transcripts", gene.num_transcripts),
            _attr("num_samples",     gene.num_samples),
        ])
        out_fh.write("\t".join([
            gene.chrom, source_tag, "gene",
            str(gene.start + 1), str(gene.end),
            ".", gene.strand, ".",
            gene_attrs,
        ]) + "\n")

        for mt in gene.transcripts:
            _write_transcript(mt, out_fh, source_tag)


def _write_transcript(mt: MergedTranscript, out_fh: TextIO, source_tag: str):
    jcoords = ",".join(f"{j.donor}-{j.acceptor}" for j in mt.junctions) or "NA"
    j_known = ",".join("known" if j.is_known else "novel" for j in mt.junctions) or "NA"
    splice  = ",".join(mt.splice_sites) if mt.splice_sites else "NA"
    per_j_diff = ",".join(
        f"j{i}:{v}bp" for i, v in sorted(mt.per_junction_max_diff.items())
    ) or "NA"
    per_j_support = ",".join(
        f"j{i}:{':'.join(s)}" for i, s in sorted(mt.per_junction_support.items())
    ) or "NA"

    trans_attrs = "; ".join([
        _attr("gene_id",              mt.gene_id),
        _attr("transcript_id",        mt.transcript_id),
        _attr("merge_status",         mt.merge_status.value),
        _attr("class_code",           mt.class_code.value),
        _attr("ref_transcript_id",    mt.ref_transcript_id),
        _attr("ref_gene_id",          mt.ref_gene_id),
        _attr("num_samples",          mt.num_samples),
        _attr("sample_ids",           ",".join(mt.sample_ids)),
        _attr("freq_score",           f"{mt.freq_score:.4f}"),
        _attr("tss_coord",            mt.tss_coord),
        _attr("tes_coord",            mt.tes_coord),
        _attr("tss_source",           mt.tss_source),
        _attr("tes_source",           mt.tes_source),
        _attr("tss_diff_range",       mt.tss_diff_range),
        _attr("tes_diff_range",       mt.tes_diff_range),
        _attr("truncation_5_count",   mt.truncation_5_count),
        _attr("truncation_3_count",   mt.truncation_3_count),
        _attr("truncation_5_ratio",   f"{mt.truncation_5_ratio:.4f}"),
        _attr("truncation_3_ratio",   f"{mt.truncation_3_ratio:.4f}"),
        _attr("missing_junctions_5",  mt.missing_junctions_5),
        _attr("missing_junctions_3",  mt.missing_junctions_3),
        _attr("num_junctions",        mt.num_junctions),
        _attr("junction_coords",      jcoords),
        _attr("junction_known",       j_known),
        _attr("per_junction_support", per_j_support),
        _attr("per_junction_max_diff",per_j_diff),
        _attr("diff_junction_count",  mt.diff_junction_count),
        _attr("diff_junction_ratio",  f"{mt.diff_junction_ratio:.4f}"),
        _attr("strand_status",        mt.strand_status.value),
    ])
    out_fh.write("\t".join([
        mt.chrom, source_tag, "transcript",
        str(mt.start + 1), str(mt.end),
        ".", mt.strand, ".",
        trans_attrs,
    ]) + "\n")

    for ei, exon in enumerate(mt.exons):
        exon_attrs = "; ".join([
            _attr("gene_id",       mt.gene_id),
            _attr("transcript_id", mt.transcript_id),
            _attr("exon_number",   ei + 1),
        ])
        out_fh.write("\t".join([
            exon.chrom, source_tag, "exon",
            str(exon.start + 1), str(exon.end),
            ".", exon.strand, ".",
            exon_attrs,
        ]) + "\n")


# ─────────────────────────────────────────────────────────────────
# Tracking writer
# ─────────────────────────────────────────────────────────────────

TRACKING_HEADER = "\t".join([
    "transcript_id", "gene_id",
    "merge_status", "class_code",
    "ref_transcript_id", "ref_gene_id",
    "num_samples", "sample_ids",
    "source_transcript_ids",
    "per_sample_status",
    "num_junctions", "junction_coords",
    "junction_known_status",
    "per_junction_sample_support",
    "per_junction_max_coord_diff",
    "diff_junction_count", "diff_junction_ratio",
    "tss_coord", "tes_coord",
    "tss_source", "tes_source",
    "tss_diff_range", "tes_diff_range",
    "truncation_5_count", "truncation_3_count",
    "truncation_5_ratio", "truncation_3_ratio",
    "missing_junctions_5", "missing_junctions_3",
    "freq_score", "strand_status",
])


def write_tracking(merged: List[MergedTranscript], out_fh: TextIO):
    out_fh.write(TRACKING_HEADER + "\n")
    for mt in merged:
        jcoords = ",".join(f"{j.donor}-{j.acceptor}" for j in mt.junctions) or "NA"
        j_known = ",".join("known" if j.is_known else "novel" for j in mt.junctions) or "NA"
        per_j_support = ",".join(
            f"j{i}:{':'.join(s)}" for i, s in sorted(mt.per_junction_support.items())
        ) or "NA"
        per_j_diff = ",".join(
            f"j{i}:{v}bp" for i, v in sorted(mt.per_junction_max_diff.items())
        ) or "NA"
        per_sample_status = ",".join(f"{s}:{v}" for s, v in mt.per_sample_status.items())
        src_ids = ",".join(mt.source_transcripts)

        row = [
            mt.transcript_id, mt.gene_id,
            mt.merge_status.value, mt.class_code.value,
            mt.ref_transcript_id, mt.ref_gene_id,
            str(mt.num_samples), ",".join(mt.sample_ids),
            src_ids, per_sample_status,
            str(mt.num_junctions), jcoords, j_known,
            per_j_support, per_j_diff,
            str(mt.diff_junction_count), f"{mt.diff_junction_ratio:.4f}",
            str(mt.tss_coord), str(mt.tes_coord),
            mt.tss_source, mt.tes_source,
            mt.tss_diff_range, mt.tes_diff_range,
            str(mt.truncation_5_count), str(mt.truncation_3_count),
            f"{mt.truncation_5_ratio:.4f}", f"{mt.truncation_3_ratio:.4f}",
            mt.missing_junctions_5, mt.missing_junctions_3,
            f"{mt.freq_score:.4f}", mt.strand_status.value,
        ]
        out_fh.write("\t".join(row) + "\n")


# ─────────────────────────────────────────────────────────────────
# Cluster writers
# ─────────────────────────────────────────────────────────────────

def write_clusters(clusters, out_fh: TextIO):
    """Write cluster membership in cd-hit-like format."""
    header = "\t".join([
        "cluster_id", "transcript_id", "is_representative",
        "similarity", "structural", "junction", "jaccard",
        "class_code", "ref_id", "chrom", "strand",
    ])
    out_fh.write(header + "\n")

    for cluster in clusters:
        for member in cluster.members:
            row = [
                cluster.cluster_id,
                member.transcript_id,
                "yes" if member.is_representative else "no",
                f"{member.similarity.total:.4f}",
                f"{member.similarity.structural:.4f}",
                f"{member.similarity.junction:.4f}",
                f"{member.similarity.jaccard:.4f}",
                member.class_code.value,
                member.ref_id,
                cluster.chrom,
                cluster.strand,
            ]
            out_fh.write("\t".join(row) + "\n")


def write_similarity_pairs(pairs, out_fh: TextIO):
    """Write all-pairs similarity (optional, can be very large)."""
    header = "\t".join([
        "transcript_a", "transcript_b",
        "similarity", "structural", "junction", "jaccard",
    ])
    out_fh.write(header + "\n")
    for tid_a, tid_b, sim in pairs:
        out_fh.write("\t".join([
            tid_a, tid_b,
            f"{sim.total:.4f}",
            f"{sim.structural:.4f}",
            f"{sim.junction:.4f}",
            f"{sim.jaccard:.4f}",
        ]) + "\n")
