"""Configuration dataclasses for merge and cluster commands."""
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class MergeConfig:
    """Config for `gtf-merge merge` subcommand."""
    # I/O
    input:           str = ""
    output:          str = ""
    ref_gtf:         Optional[str] = None
    novel_prefix:    str = "Novel"
    output_extended: bool = False
    chr_style:       str = "auto"

    # platform
    platform:           str = "pacbio"
    junction_threshold: int = 3

    # multi-exon merge
    max_5trunc:     int   = 1
    max_3trunc:     int   = 0
    tss_window:     int   = 200
    tes_window:     int   = 50
    max_diff_count: int   = -1     # -1 = unlimited
    max_diff_ratio: float = 0.5

    # mono-exon merge
    mono_overlap: float = 0.5
    mono_5end:    int   = 500
    mono_3end:    int   = 50

    # representative selection
    tss_strategy: str = "longest"
    tes_strategy: str = "freq"

    # filtering
    min_samples:    int = 1
    min_transcript_len: int = 100

    # output
    compress:       bool = False
    output_tracking: bool = True

    # performance
    threads:  int = 4
    tmp_dir:  str = "/tmp"
    keep_tmp: bool = False

    # checkpoint
    checkpoint_dir: str = ""
    resume:         bool = False
    force:          bool = False

    # logging
    verbose:  bool = False
    log_file: str  = ""

    def __post_init__(self):
        if self.platform == "ont" and self.junction_threshold == 3:
            self.junction_threshold = 5
        if not self.checkpoint_dir and self.output:
            self.checkpoint_dir = self.output + "_checkpoint"


@dataclass
class ClusterConfig:
    """Config for `gtf-merge cluster` subcommand."""
    # I/O
    input:        str = ""
    output:       str = ""
    ref_gtf:      Optional[str] = None
    chr_style:    str = "auto"

    # similarity
    junction_threshold: int = 3
    fuzzy_threshold:    int = 50
    weight_structural:  float = 0.4
    weight_junction:    float = 0.4
    weight_jaccard:     float = 0.2

    # cluster
    min_similarity:        float = 0.7
    representative_order:  str   = "ref_first"   # ref_first / length / samples / random

    # output
    output_pairs:  bool = True

    # performance
    threads:  int = 4
    tmp_dir:  str = "/tmp"

    # logging
    verbose:  bool = False
    log_file: str  = ""
