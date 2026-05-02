"""Input validation utilities."""
import logging
import os
from pathlib import Path
from typing import List, Dict, Tuple

logger = logging.getLogger(__name__)


def validate_inputs(samples: List[Dict], ref_gtf: str = None) -> Tuple[bool, List[str]]:
    """
    Validate all input files.
    Returns (all_ok, list_of_errors).
    """
    errors = []

    # check for unique sample IDs
    ids = [s["sample_id"] for s in samples]
    seen = set()
    dupes = []
    for sid in ids:
        if sid in seen:
            dupes.append(sid)
        seen.add(sid)
    if dupes:
        errors.append(f"Duplicate sample IDs: {', '.join(dupes)}")

    # check GTF files exist
    for s in samples:
        p = s["gtf_path"]
        if not Path(p).exists():
            errors.append(f"GTF file not found: {p}")
        elif Path(p).stat().st_size == 0:
            errors.append(f"GTF file is empty: {p}")

    # check ref_gtf if provided
    if ref_gtf:
        if not Path(ref_gtf).exists():
            errors.append(f"Reference GTF not found: {ref_gtf}")

    # check cap_status values
    for s in samples:
        if s["cap_status"] not in ("capped", "no_cap"):
            errors.append(f"Invalid cap_status '{s['cap_status']}' for sample {s['sample_id']}")

    # check priorities are integers
    for s in samples:
        if not isinstance(s["priority"], int):
            errors.append(f"Priority must be integer for sample {s['sample_id']}")

    all_ok = len(errors) == 0
    return all_ok, errors


def detect_chr_style(gtf_path: str) -> str:
    """Auto-detect chromosome naming style from GTF."""
    import gzip
    from pathlib import Path

    opener = gzip.open if gtf_path.endswith(".gz") else open
    mode   = "rt" if gtf_path.endswith(".gz") else "r"

    try:
        with opener(gtf_path, mode) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                chrom = line.split("\t")[0]
                if chrom.startswith("chr"):
                    return "ucsc"
                else:
                    return "ensembl"
    except Exception:
        pass
    return "ensembl"


def check_chr_consistency(samples: List[Dict]) -> Tuple[bool, str]:
    """
    Check that all samples use the same chromosome naming convention.
    Returns (consistent, detected_style).
    """
    styles = set()
    for s in samples[:min(5, len(samples))]:  # check first 5 samples
        style = detect_chr_style(s["gtf_path"])
        styles.add(style)

    if len(styles) > 1:
        return False, "mixed"
    return True, styles.pop() if styles else "unknown"
