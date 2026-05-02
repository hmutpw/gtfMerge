"""Core data models for GTF Merge.

Two-step architecture:
  Step1 (merge):   Merge transcripts likely representing the same molecule
  Step2 (cluster): Group structurally similar transcripts (no forced merge)
"""
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Optional
from enum import Enum


# ─────────────────────────────────────────────────────────────────
# Enums
# ─────────────────────────────────────────────────────────────────

class MergeStatus(str, Enum):
    """How a transcript ended up in the merged set."""
    EXACT_MATCH        = "exact_match"        # all junctions match (within threshold)
    TRUNCATION_5       = "truncation_5"
    TRUNCATION_3       = "truncation_3"
    TRUNCATION_BOTH    = "truncation_both"
    FUZZY_MATCH        = "fuzzy_match"        # within threshold but coords drift
    MONOEXON_EXACT     = "monoexon_exact"
    MONOEXON_OVERLAP   = "monoexon_overlap"
    SINGLETON          = "singleton"
    REF_ONLY           = "ref_only"           # only in reference, not in samples


class ClassCode(str, Enum):
    """gffcompare-style class codes for query vs reference."""
    EQUAL                = "="   # exact match
    CONTAINED            = "c"   # query contained in ref
    CONTAINING           = "k"   # query contains ref
    RETAINED_INTRON      = "m"   # all introns of ref retained in query
    PARTIAL_RETAINED     = "n"   # partial intron retained
    JUNCTION_MATCH       = "j"   # multi-exon, ≥1 matched junction
    SINGLE_OVERLAP       = "e"   # single-exon overlap with multi-exon ref
    GENERIC_OVERLAP      = "o"   # generic overlap, same strand
    INTRON_MATCH_OPP     = "s"   # intron match on opposite strand
    EXONIC_OPP           = "x"   # exonic overlap on opposite strand
    INTRONIC             = "i"   # query fully within ref intron
    REF_WITHIN_QUERY     = "y"   # ref fully within query (incl introns)
    POLYMERASE_RUNON     = "p"   # possible polymerase run-on
    UNKNOWN              = "u"   # unknown/intergenic


# Class code priority (lower = higher priority, gffcompare convention)
CLASS_CODE_PRIORITY = {
    "=": 1,  "c": 2,  "k": 3,  "m": 4,  "n": 5,
    "j": 6,  "e": 7,  "o": 8,  "s": 9,  "x": 10,
    "i": 11, "y": 12, "p": 13, "u": 14,
}


class StrandStatus(str, Enum):
    KNOWN    = "known"
    UNKNOWN  = "unknown"
    RESOLVED = "resolved"


# ─────────────────────────────────────────────────────────────────
# Junction / Intron / Exon
# ─────────────────────────────────────────────────────────────────

@dataclass
class Junction:
    """A splice junction (donor-acceptor pair, 0-based).
    Equivalent to an intron coordinate-wise.
    """
    chrom:    str
    donor:    int   # end of upstream exon (exclusive, 0-based)
    acceptor: int   # start of downstream exon (0-based)
    strand:   str

    is_known:    bool  = False
    splice_site: str   = "NA"   # GT-AG / GC-AG / AT-AC / non-canonical

    def key(self) -> Tuple:
        return (self.chrom, self.strand, self.donor, self.acceptor)

    def coord_error(self, other: "Junction") -> int:
        return max(abs(self.donor - other.donor),
                   abs(self.acceptor - other.acceptor))

    def matches(self, other: "Junction", threshold: int) -> bool:
        if self.chrom != other.chrom:
            return False
        if self.strand != "." and other.strand != "." and self.strand != other.strand:
            return False
        return self.coord_error(other) <= threshold


# Aliases for clarity
Intron = Junction


@dataclass
class Exon:
    chrom:  str
    start:  int   # 0-based
    end:    int   # 0-based exclusive
    strand: str

    @property
    def length(self) -> int:
        return self.end - self.start


# ─────────────────────────────────────────────────────────────────
# Transcript (input from one sample)
# ─────────────────────────────────────────────────────────────────

@dataclass
class Transcript:
    transcript_id: str
    gene_id:       str
    chrom:         str
    start:         int    # 0-based, leftmost
    end:           int    # 0-based exclusive, rightmost
    strand:        str    # +, -, .
    source_sample: str
    source_file:   str

    exons:      List[Exon]     = field(default_factory=list)
    junctions:  List[Junction] = field(default_factory=list)
    attributes: Dict[str, str] = field(default_factory=dict)

    is_monoexon:    bool = False
    junction_chain: str  = ""    # canonical key for exact matching

    def __post_init__(self):
        self.exons = sorted(self.exons, key=lambda e: e.start)
        self.is_monoexon = (len(self.exons) <= 1)
        if not self.junctions:
            self.junctions = self._build_junctions()
        self.junction_chain = self._make_chain()

    def _build_junctions(self) -> List[Junction]:
        return [
            Junction(
                chrom    = self.exons[i].chrom,
                donor    = self.exons[i].end,
                acceptor = self.exons[i+1].start,
                strand   = self.strand,
            )
            for i in range(len(self.exons) - 1)
        ]

    def _make_chain(self) -> str:
        if not self.junctions:
            return ""
        return f"{self.chrom}:{self.strand}:" + ":".join(
            f"{j.donor}-{j.acceptor}" for j in self.junctions
        )

    @property
    def tss(self) -> int:
        return self.start if self.strand != "-" else self.end

    @property
    def tes(self) -> int:
        return self.end if self.strand != "-" else self.start

    @property
    def num_junctions(self) -> int:
        return len(self.junctions)

    @property
    def length(self) -> int:
        return sum(e.length for e in self.exons)

    @property
    def genomic_span(self) -> int:
        return self.end - self.start


# ─────────────────────────────────────────────────────────────────
# Merged transcript (output of Step1)
# ─────────────────────────────────────────────────────────────────

@dataclass
class MergedTranscript:
    """Output of Step1 merge.
    Represents a consensus transcript from one or more sample transcripts."""
    transcript_id: str        # final ID (Novel.Tx... or reference ID)
    gene_id:       str        # gene ID (novel or reference)
    chrom:         str
    start:         int
    end:           int
    strand:        str

    # structure
    exons:     List[Exon]     = field(default_factory=list)
    junctions: List[Junction] = field(default_factory=list)

    # source tracking
    source_transcripts: List[str] = field(default_factory=list)  # "sample:original_id"
    sample_ids:         List[str] = field(default_factory=list)
    num_samples:        int = 0

    # merge metadata
    merge_status:    MergeStatus = MergeStatus.SINGLETON
    per_sample_status: Dict[str, str] = field(default_factory=dict)

    # TSS/TES consensus
    tss_coord:        int = 0
    tes_coord:        int = 0
    tss_source:       str = ""
    tes_source:       str = ""
    tss_diff_range:   str = "0-0bp"
    tes_diff_range:   str = "0-0bp"

    # truncation tracking
    truncation_5_count: int = 0
    truncation_3_count: int = 0
    missing_junctions_5: str = "NA"
    missing_junctions_3: str = "NA"

    # junction support
    per_junction_support:  Dict[int, List[str]] = field(default_factory=dict)
    per_junction_max_diff: Dict[int, int]       = field(default_factory=dict)
    junction_known:        List[bool]           = field(default_factory=list)
    splice_sites:          List[str]            = field(default_factory=list)

    # micro-difference tracking
    diff_junction_count:   int   = 0   # junctions with non-zero coord error
    diff_junction_ratio:   float = 0.0

    # reference comparison
    class_code:        ClassCode = ClassCode.UNKNOWN
    ref_transcript_id: str       = "NA"
    ref_gene_id:       str       = "NA"
    ref_gene_name:     str       = "NA"

    # strand resolution
    strand_status: StrandStatus = StrandStatus.KNOWN

    # frequency
    freq_score: float = 0.0

    @property
    def truncation_5_ratio(self) -> float:
        return self.truncation_5_count / max(self.num_samples, 1)

    @property
    def truncation_3_ratio(self) -> float:
        return self.truncation_3_count / max(self.num_samples, 1)

    @property
    def num_junctions(self) -> int:
        return len(self.junctions)

    @property
    def is_monoexon(self) -> bool:
        return len(self.exons) <= 1


# ─────────────────────────────────────────────────────────────────
# Cluster (output of Step2)
# ─────────────────────────────────────────────────────────────────

@dataclass
class SimilarityScore:
    """Similarity score between two transcripts (Step2)."""
    structural: float = 0.0   # exon count + length similarity
    junction:   float = 0.0   # junction matching
    jaccard:    float = 0.0   # base-level Jaccard
    total:      float = 0.0


@dataclass
class ClusterMember:
    """A transcript belonging to a cluster."""
    transcript_id: str
    similarity:    SimilarityScore = field(default_factory=SimilarityScore)
    class_code:    ClassCode = ClassCode.UNKNOWN
    ref_id:        str = "NA"
    is_representative: bool = False


@dataclass
class Cluster:
    cluster_id:     str
    chrom:          str
    strand:         str
    representative: ClusterMember
    members:        List[ClusterMember] = field(default_factory=list)

    @property
    def size(self) -> int:
        return len(self.members)


# ─────────────────────────────────────────────────────────────────
# Gene
# ─────────────────────────────────────────────────────────────────

@dataclass
class Gene:
    gene_id:       str
    chrom:         str
    start:         int
    end:           int
    strand:        str
    transcripts:   List[MergedTranscript] = field(default_factory=list)
    ref_gene_id:   str = "NA"
    ref_gene_name: str = "NA"

    @property
    def num_transcripts(self) -> int:
        return len(self.transcripts)

    @property
    def num_samples(self) -> int:
        s = set()
        for t in self.transcripts:
            s.update(t.sample_ids)
        return len(s)
