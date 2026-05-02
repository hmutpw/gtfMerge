"""Streaming GTF parser - memory efficient, handles gzip."""
import gzip
import logging
from pathlib import Path
from typing import Iterator, Dict, List, Optional, Tuple
from collections import defaultdict

from gtf_merge.core.models import Transcript, Exon, Junction

logger = logging.getLogger(__name__)


def _open_gtf(path: str):
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")


def _parse_attributes(attr_str: str) -> Dict[str, str]:
    """Parse GTF attribute string into dict."""
    attrs = {}
    for part in attr_str.strip().rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        if " " in part:
            key, _, val = part.partition(" ")
            attrs[key] = val.strip('"')
    return attrs


def parse_gtf_by_chrom(
    gtf_path: str,
    target_chrom: Optional[str] = None,
) -> Iterator[Transcript]:
    """
    Stream-parse a GTF file, yielding complete Transcript objects.
    If target_chrom is set, only that chromosome is parsed.
    """
    # Buffer: transcript_id -> raw data
    trans_buf: Dict[str, dict] = {}
    exon_buf:  Dict[str, List[Tuple]] = defaultdict(list)

    def _flush(tid: str, sample_id: str, source_file: str) -> Optional[Transcript]:
        if tid not in trans_buf:
            return None
        td = trans_buf.pop(tid)
        raw_exons = exon_buf.pop(tid, [])
        if not raw_exons:
            return None

        exons = []
        for (chrom, s, e, strand) in raw_exons:
            exons.append(Exon(chrom=chrom, start=s, end=e, strand=strand))
        exons.sort(key=lambda x: x.start)

        # build junctions from consecutive exons
        junctions = []
        for i in range(len(exons) - 1):
            junctions.append(Junction(
                chrom    = exons[i].chrom,
                donor    = exons[i].end,       # 0-based end of upstream exon
                acceptor = exons[i+1].start,   # 0-based start of downstream exon
                strand   = exons[i].strand,
            ))

        t = Transcript(
            transcript_id  = tid,
            gene_id        = td["gene_id"],
            chrom          = td["chrom"],
            start          = td["start"],
            end            = td["end"],
            strand         = td["strand"],
            source_sample  = sample_id,
            source_file    = source_file,
            exons          = exons,
            junctions      = junctions,
            attributes     = td.get("attrs", {}),
        )
        return t

    sample_id   = Path(gtf_path).stem
    source_file = gtf_path
    prev_tid    = None

    try:
        with _open_gtf(gtf_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                line = line.rstrip("\n")
                if not line:
                    continue

                fields = line.split("\t")
                if len(fields) < 9:
                    continue

                chrom, source, feature, start, end, score, strand, frame, attrs_str = fields
                start = int(start) - 1   # convert to 0-based
                end   = int(end)         # GTF end is inclusive, we store exclusive

                if target_chrom and chrom != target_chrom:
                    continue

                attrs = _parse_attributes(attrs_str)
                gene_id  = attrs.get("gene_id", "")
                trans_id = attrs.get("transcript_id", "")

                if feature == "transcript":
                    # flush previous if changed
                    if prev_tid and prev_tid != trans_id:
                        t = _flush(prev_tid, sample_id, source_file)
                        if t:
                            yield t

                    trans_buf[trans_id] = {
                        "gene_id": gene_id,
                        "chrom":   chrom,
                        "start":   start,
                        "end":     end,
                        "strand":  strand if strand in ("+", "-") else ".",
                        "attrs":   attrs,
                    }
                    prev_tid = trans_id

                elif feature == "exon":
                    if trans_id:
                        exon_buf[trans_id].append((chrom, start, end, strand))
                    if trans_id not in trans_buf and gene_id:
                        # GTF without explicit transcript lines
                        trans_buf[trans_id] = {
                            "gene_id": gene_id,
                            "chrom":   chrom,
                            "start":   start,
                            "end":     end,
                            "strand":  strand if strand in ("+", "-") else ".",
                            "attrs":   attrs,
                        }
                        prev_tid = trans_id

            # flush last
            if prev_tid:
                t = _flush(prev_tid, sample_id, source_file)
                if t:
                    yield t

    except Exception as e:
        logger.error(f"Error parsing {gtf_path}: {e}")
        raise


def get_chromosomes(gtf_path: str) -> List[str]:
    """Extract all chromosome names from a GTF file."""
    chroms = []
    seen   = set()
    try:
        with _open_gtf(gtf_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 1:
                    continue
                c = fields[0]
                if c not in seen:
                    seen.add(c)
                    chroms.append(c)
    except Exception as e:
        logger.warning(f"Could not read chromosomes from {gtf_path}: {e}")
    return chroms


def parse_filelist(filelist_path: str) -> List[Dict]:
    """
    Parse the input file list.
    Format (tab-separated, with header):
      gtf_path  sample_id  cap_status  priority  metadata
    """
    samples = []
    with open(filelist_path) as fh:
        for i, line in enumerate(fh):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[0].lower() in ("gtf_path", "gtf", "path"):
                continue  # skip header

            if len(fields) < 1:
                continue

            gtf_path   = fields[0].strip()
            sample_id  = fields[1].strip() if len(fields) > 1 else Path(gtf_path).stem
            cap_status = fields[2].strip() if len(fields) > 2 else "no_cap"
            priority   = int(fields[3].strip()) if len(fields) > 3 else 1
            metadata   = fields[4].strip() if len(fields) > 4 else ""

            if cap_status not in ("capped", "no_cap"):
                logger.warning(f"Line {i+1}: unknown cap_status '{cap_status}', defaulting to no_cap")
                cap_status = "no_cap"

            samples.append({
                "gtf_path":   gtf_path,
                "sample_id":  sample_id,
                "cap_status": cap_status,
                "priority":   priority,
                "metadata":   metadata,
            })

    return samples
