"""Step 3 (annotation): Compare merged transcripts against reference annotation."""
import logging
from collections import defaultdict
from typing import List, Dict, Optional, Tuple

from gtf_merge.core.models import MergedTranscript, RefMatchType
from gtf_merge.core.parser import parse_gtf_by_chrom
from gtf_merge.core.similarity import match_introns as match_junctions

logger = logging.getLogger(__name__)


class RefAnnotation:
    """Holds reference transcripts indexed by chromosome for fast lookup."""

    def __init__(self):
        # chrom -> list of ref transcript dicts
        self._by_chrom: Dict[str, List[dict]] = defaultdict(list)

    def load(self, ref_gtf: str):
        """Load all reference transcripts from GTF."""
        logger.info(f"Loading reference annotation: {ref_gtf}")
        count = 0
        from gtf_merge.core.parser import parse_gtf_by_chrom
        # we read without chrom filter to build full index
        # use a dummy sample name
        for t in _iter_all_chroms(ref_gtf):
            self._by_chrom[t.chrom].append({
                "transcript_id": t.transcript_id,
                "gene_id":       t.gene_id,
                "gene_name":     t.attributes.get("gene_name", "NA"),
                "chrom":         t.chrom,
                "start":         t.start,
                "end":           t.end,
                "strand":        t.strand,
                "junctions":     t.junctions,
                "exons":         t.exons,
                "junction_chain": t.junction_chain,
            })
            count += 1
        logger.info(f"Loaded {count} reference transcripts")

    def get_chrom(self, chrom: str) -> List[dict]:
        return self._by_chrom.get(chrom, [])


def _iter_all_chroms(ref_gtf: str):
    """Yield all transcripts from reference GTF using the main parser."""
    import gzip
    from pathlib import Path
    from gtf_merge.core.parser import _open_gtf, _parse_attributes
    from gtf_merge.core.models import Exon, Junction, Transcript
    from collections import defaultdict

    trans_buf: Dict[str, dict] = {}
    exon_buf:  Dict[str, List] = defaultdict(list)
    prev_tid  = None

    def flush(tid):
        if tid not in trans_buf:
            return None
        td = trans_buf.pop(tid)
        raw = exon_buf.pop(tid, [])
        if not raw:
            return None
        exons = sorted([Exon(chrom=c, start=s, end=e, strand=st) for c,s,e,st in raw], key=lambda x: x.start)
        junctions = []
        for i in range(len(exons)-1):
            junctions.append(Junction(
                chrom=exons[i].chrom, donor=exons[i].end,
                acceptor=exons[i+1].start, strand=exons[i].strand,
            ))
        return Transcript(
            transcript_id=tid, gene_id=td["gene_id"],
            chrom=td["chrom"], start=td["start"], end=td["end"],
            strand=td["strand"], source_sample="ref", source_file=ref_gtf,
            exons=exons, junctions=junctions, attributes=td.get("attrs", {}),
        )

    with _open_gtf(ref_gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs_str = fields
            start = int(start) - 1
            end   = int(end)
            attrs = _parse_attributes(attrs_str)
            tid   = attrs.get("transcript_id", "")
            gid   = attrs.get("gene_id", "")

            if feature == "transcript":
                if prev_tid and prev_tid != tid:
                    t = flush(prev_tid)
                    if t:
                        yield t
                trans_buf[tid] = {
                    "gene_id": gid, "chrom": chrom,
                    "start": start, "end": end,
                    "strand": strand if strand in ("+", "-") else ".",
                    "attrs": attrs,
                }
                prev_tid = tid
            elif feature == "exon" and tid:
                exon_buf[tid].append((chrom, start, end, strand))
                if tid not in trans_buf:
                    trans_buf[tid] = {
                        "gene_id": gid, "chrom": chrom,
                        "start": start, "end": end,
                        "strand": strand if strand in ("+", "-") else ".",
                        "attrs": attrs,
                    }
                    prev_tid = tid

    if prev_tid:
        t = flush(prev_tid)
        if t:
            yield t


def compare_to_ref(
    merged:    List[MergedTranscript],
    ref_annot: RefAnnotation,
    threshold: int = 0,
) -> List[MergedTranscript]:
    """
    Compare merged transcripts to reference annotation.
    Annotates ref_match_type, ref_gene_id, ref_transcript_id in-place.
    """
    for mt in merged:
        ref_list = ref_annot.get_chrom(mt.chrom)
        if not ref_list:
            mt.ref_match_type = RefMatchType.NOVEL
            continue

        best_match_type = RefMatchType.NOVEL
        best_ref_gene   = "NA"
        best_ref_trans  = "NA"

        for ref in ref_list:
            # quick overlap check
            if ref["end"] < mt.start or ref["start"] > mt.end:
                continue
            if ref["strand"] != mt.strand and mt.strand != ".":
                continue

            # full match: junction chain identical
            if not mt.is_monoexon:
                if mt.junctions and ref["junctions"]:
                    ref_chain = (
                        f"{ref['chrom']}:{ref['strand']}:" +
                        ":".join(f"{j.donor}-{j.acceptor}" for j in ref["junctions"])
                    )
                    my_chain = (
                        f"{mt.chrom}:{mt.strand}:" +
                        ":".join(f"{j.donor}-{j.acceptor}" for j in mt.junctions)
                    )
                    if my_chain == ref_chain:
                        best_match_type = RefMatchType.FULL_MATCH
                        best_ref_gene   = ref["gene_id"]
                        best_ref_trans  = ref["transcript_id"]
                        break  # can't do better

                    # junction match (all query junctions found in ref)
                    matched = match_junctions(ref["junctions"], mt.junctions, threshold)
                    if len(matched) == len(mt.junctions):
                        if best_match_type.value not in ("full_match",):
                            best_match_type = RefMatchType.JUNCTION_MATCH
                            best_ref_gene   = ref["gene_id"]
                            best_ref_trans  = ref["transcript_id"]
                    elif matched:
                        if best_match_type.value not in ("full_match", "junction_match"):
                            best_match_type = RefMatchType.PARTIAL_MATCH
                            best_ref_gene   = ref["gene_id"]
                            best_ref_trans  = ref["transcript_id"]
                    else:
                        if best_match_type == RefMatchType.NOVEL:
                            best_match_type = RefMatchType.GENE_MATCH
                            best_ref_gene   = ref["gene_id"]
            else:
                # monoexon: overlap
                if best_match_type == RefMatchType.NOVEL:
                    best_match_type = RefMatchType.GENE_MATCH
                    best_ref_gene   = ref["gene_id"]

        mt.ref_match_type   = best_match_type
        mt.ref_gene_id      = best_ref_gene
        mt.ref_transcript_id = best_ref_trans

    return merged
