"""Per-chromosome processing for merge command."""
import logging
from collections import defaultdict
from typing import List, Dict, Optional, Tuple

from gtf_merge.core.models import (
    Transcript, MergedTranscript, ClassCode, MergeStatus, CLASS_CODE_PRIORITY,
)
from gtf_merge.core.parser import parse_gtf_by_chrom
from gtf_merge.core.merger import run_merge
from gtf_merge.core.class_code import RefIndex, classify_transcript
from gtf_merge.core.gene_definer import define_genes
from gtf_merge.core.writer import write_merged_gtf, write_tracking

logger = logging.getLogger(__name__)


def _generate_id(novel_prefix: str, counter: int, chrom: str, code: ClassCode) -> str:
    """Generate Novel.Tx000001.chr1.j format."""
    return f"{novel_prefix}.Tx{counter:06d}.{chrom}.{code.value}"


def assign_ids_and_classify(
    merged_transcripts: List[MergedTranscript],
    ref_index: Optional[RefIndex],
    novel_prefix: str,
    junction_threshold: int,
    id_counter_start: int = 1,
) -> int:
    """
    Assign final transcript_id and gene_id to each merged transcript.
    Compares against reference (if provided) to get class_code.

    Returns next available counter.
    """
    counter = id_counter_start
    novel_gene_counter = 1

    # gene_id mapping for novel ones (chrom -> next idx)
    novel_genes_by_chrom: Dict[str, Dict[str, str]] = defaultdict(dict)

    for mt in merged_transcripts:
        # classify against reference
        if ref_index:
            class_code, ref_t = classify_transcript(mt, ref_index, junction_threshold)
        else:
            class_code, ref_t = ClassCode.UNKNOWN, None

        mt.class_code = class_code

        if ref_t:
            mt.ref_transcript_id = ref_t.transcript_id
            mt.ref_gene_id       = ref_t.gene_id
            mt.ref_gene_name     = ref_t.attributes.get("gene_name", "NA")

        # assign final transcript_id
        if class_code == ClassCode.EQUAL and ref_t:
            mt.transcript_id = ref_t.transcript_id
            mt.gene_id       = ref_t.gene_id
        else:
            mt.transcript_id = _generate_id(novel_prefix, counter, mt.chrom, class_code)
            counter += 1

            # gene_id assignment
            if class_code in (ClassCode.CONTAINED, ClassCode.CONTAINING,
                              ClassCode.JUNCTION_MATCH, ClassCode.RETAINED_INTRON,
                              ClassCode.PARTIAL_RETAINED) and ref_t:
                # novel isoform of known gene
                mt.gene_id = ref_t.gene_id
            else:
                # need a novel gene
                # use rough overlap to group: if same novel transcripts overlap → same gene
                # for now use simple per-chrom counter
                novel_gid = f"{novel_prefix}.G{novel_gene_counter:06d}.{mt.chrom}"
                novel_gene_counter += 1
                mt.gene_id = novel_gid

    return counter


def process_chrom(
    chrom:            str,
    samples:          List[Dict],
    cfg,
    ref_index:        Optional[RefIndex],
    output_gtf:       str,
    output_tracking:  str,
    id_counter_start: int = 1,
) -> Dict:
    """Process all transcripts on one chromosome."""
    logger.info(f"[{chrom}] Loading transcripts from {len(samples)} samples...")

    all_transcripts: List[Transcript] = []
    for s in samples:
        for t in parse_gtf_by_chrom(s["gtf_path"], target_chrom=chrom):
            t.source_sample = s["sample_id"]
            all_transcripts.append(t)

    if not all_transcripts:
        logger.info(f"[{chrom}] No transcripts found")
        open(output_gtf, "w").close()
        open(output_tracking, "w").close()
        return {"chrom": chrom, "input": 0, "output": 0, "next_counter": id_counter_start}

    total_input = len(all_transcripts)
    logger.info(f"[{chrom}] {total_input} transcripts loaded")

    total_samples = len(samples)

    # ── Step1: Merge ─────────────────────────────────────────────
    merged = run_merge(all_transcripts, cfg, total_samples)
    logger.info(f"[{chrom}] Merged {total_input} → {len(merged)}")

    # ── Step2: Classify against reference ────────────────────────
    next_counter = assign_ids_and_classify(
        merged, ref_index, cfg.novel_prefix,
        cfg.junction_threshold,
        id_counter_start = id_counter_start,
    )

    # ── Step3: Define genes ──────────────────────────────────────
    genes = define_genes(merged, id_prefix=cfg.novel_prefix)
    logger.info(f"[{chrom}] {len(genes)} genes")

    # ── Step4: Write outputs ─────────────────────────────────────
    with open(output_gtf, "w") as gfh, open(output_tracking, "w") as tfh:
        write_merged_gtf(genes, gfh)
        write_tracking(merged, tfh)

    return {
        "chrom":         chrom,
        "input":         total_input,
        "merged":        len(merged),
        "genes":         len(genes),
        "next_counter":  next_counter,
    }
