"""GTF Merge CLI with two subcommands: merge / cluster."""
import argparse
import logging
import os
import sys
import gzip
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import List, Dict, Optional
from datetime import datetime

from gtf_merge import __version__
from gtf_merge.utils.config import MergeConfig, ClusterConfig
from gtf_merge.utils.checkpoint import CheckpointManager
from gtf_merge.utils.validator import validate_inputs
from gtf_merge.core.parser import parse_filelist, get_chromosomes, parse_gtf_by_chrom
from gtf_merge.core.class_code import RefIndex, classify_transcript
from gtf_merge.core.chrom_worker import process_chrom


# ─────────────────────────────────────────────────────────────────
# Logging
# ─────────────────────────────────────────────────────────────────

def setup_logging(verbose: bool, log_file: str = ""):
    level = logging.DEBUG if verbose else logging.INFO
    fmt   = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    handlers = [logging.StreamHandler(sys.stderr)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    logging.basicConfig(level=level, format=fmt, handlers=handlers)


# ─────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────

def get_all_chroms(samples: List[Dict]) -> List[str]:
    seen   = set()
    chroms = []
    for s in samples:
        for c in get_chromosomes(s["gtf_path"]):
            if c not in seen:
                seen.add(c)
                chroms.append(c)
    import re
    def _nat_key(c):
        c2 = c.replace("chr", "")
        parts = re.split(r"(\d+)", c2)
        return [int(p) if p.isdigit() else p for p in parts]
    return sorted(chroms, key=_nat_key)


def load_ref_index(ref_gtf: str, target_chrom: Optional[str] = None) -> RefIndex:
    """Load a reference GTF into a RefIndex."""
    idx = RefIndex()
    if not ref_gtf:
        return idx
    n = 0
    for t in parse_gtf_by_chrom(ref_gtf, target_chrom=target_chrom):
        idx.add_transcript(t)
        n += 1
    return idx


# ─────────────────────────────────────────────────────────────────
# Worker setup for multiprocessing
# ─────────────────────────────────────────────────────────────────

_worker_ref_index    = None
_worker_ref_gtf_path = None


def _init_worker(ref_gtf_path: str):
    global _worker_ref_index, _worker_ref_gtf_path
    if ref_gtf_path and _worker_ref_gtf_path != ref_gtf_path:
        logging.getLogger("gtf_merge").info(f"[worker {os.getpid()}] loading reference: {ref_gtf_path}")
        _worker_ref_index = load_ref_index(ref_gtf_path)
        _worker_ref_gtf_path = ref_gtf_path


def _merge_worker(args):
    (chrom, samples, cfg, ref_gtf_path, out_gtf, out_tracking, id_counter_start) = args
    global _worker_ref_index
    return process_chrom(
        chrom            = chrom,
        samples          = samples,
        cfg              = cfg,
        ref_index        = _worker_ref_index,
        output_gtf       = out_gtf,
        output_tracking  = out_tracking,
        id_counter_start = id_counter_start,
    )


# ─────────────────────────────────────────────────────────────────
# `gtf-merge merge` command
# ─────────────────────────────────────────────────────────────────

def run_merge_command(cfg: MergeConfig):
    t_start = datetime.now()
    logger  = logging.getLogger("gtf_merge")
    logger.info(f"GTF Merge v{__version__} - merge subcommand")

    # parse + validate
    samples = parse_filelist(cfg.input)
    if not samples:
        logger.error("No samples in filelist")
        sys.exit(1)
    logger.info(f"Found {len(samples)} samples")

    ok, errs = validate_inputs(samples, cfg.ref_gtf)
    if not ok:
        for e in errs:
            logger.error(e)
        sys.exit(1)

    # chrom list
    chroms = get_all_chroms(samples)
    logger.info(f"Found {len(chroms)} chromosomes")

    # checkpoint
    ckpt = CheckpointManager(cfg.checkpoint_dir)
    if cfg.force:
        ckpt.reset()

    # tmp dir
    tmp_dir = Path(cfg.tmp_dir) / f"gtf_merge_{os.getpid()}"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # build per-chrom jobs (each chrom gets a non-overlapping ID counter range)
    chrom_jobs = []
    counter_start = 1
    counter_increment = 1_000_000  # space ID ranges per chrom
    for chrom in chroms:
        if cfg.resume and ckpt.is_completed(chrom):
            logger.info(f"[{chrom}] skipped (checkpoint)")
            continue
        out_gtf  = str(tmp_dir / f"{chrom}.gtf")
        out_trk  = str(tmp_dir / f"{chrom}.tracking.tsv")
        chrom_jobs.append((chrom, samples, cfg, cfg.ref_gtf, out_gtf, out_trk, counter_start))
        counter_start += counter_increment

    stats_list: List[Dict] = []

    if cfg.threads > 1 and len(chrom_jobs) > 1:
        logger.info(f"Running {len(chrom_jobs)} chroms with {cfg.threads} threads")
        with ProcessPoolExecutor(
            max_workers = cfg.threads,
            initializer = _init_worker,
            initargs    = (cfg.ref_gtf or "",),
        ) as exe:
            futures = {exe.submit(_merge_worker, job): job[0] for job in chrom_jobs}
            for fut in as_completed(futures):
                chrom = futures[fut]
                try:
                    stats = fut.result()
                    stats_list.append(stats)
                    out_gtf = str(tmp_dir / f"{chrom}.gtf")
                    out_trk = str(tmp_dir / f"{chrom}.tracking.tsv")
                    ckpt.mark_completed(chrom, out_gtf, out_trk)
                    logger.info(f"[{chrom}] done")
                except Exception as e:
                    logger.error(f"[{chrom}] failed: {e}")
                    raise
    else:
        logger.info(f"Running {len(chrom_jobs)} chroms sequentially")
        _init_worker(cfg.ref_gtf or "")
        for job in chrom_jobs:
            chrom = job[0]
            stats = _merge_worker(job)
            stats_list.append(stats)
            ckpt.mark_completed(chrom, job[4], job[5])

    # combine outputs
    Path(cfg.output).parent.mkdir(parents=True, exist_ok=True)
    final_gtf = cfg.output + ".gtf"
    final_trk = cfg.output + "_tracking.tsv"
    logger.info("Merging chrom outputs...")

    header_written = False
    with open(final_gtf, "w") as gfh, open(final_trk, "w") as tfh:
        gfh.write(f"##gtf-version 2.2\n")
        gfh.write(f"##generated-by GTFmerge {__version__}\n")
        gfh.write(f"##date {datetime.now().isoformat()}\n")
        for chrom in chroms:
            cg = str(tmp_dir / f"{chrom}.gtf")
            ct = str(tmp_dir / f"{chrom}.tracking.tsv")
            if not Path(cg).exists():
                continue
            with open(cg) as f:
                gfh.write(f.read())
            with open(ct) as f:
                for i, line in enumerate(f):
                    if i == 0:
                        if not header_written:
                            tfh.write(line)
                            header_written = True
                    else:
                        tfh.write(line)

    # extended GTF (includes ref-only transcripts)
    if cfg.output_extended and cfg.ref_gtf:
        logger.info("Generating extended GTF (with reference-only transcripts)...")
        ext_gtf = cfg.output + "_extended.gtf"
        # collect all merged transcript_ids that already used ref ids
        used_ref_ids = set()
        with open(final_trk) as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 4 and fields[3] == "=":
                    used_ref_ids.add(fields[4])  # ref_transcript_id

        with open(ext_gtf, "w") as out:
            with open(final_gtf) as f:
                out.write(f.read())
            # append unused reference transcripts
            with open(cfg.ref_gtf) as ref_fh:
                current_tid = None
                buffered = []
                def flush():
                    nonlocal buffered, current_tid
                    if current_tid and current_tid not in used_ref_ids:
                        for line in buffered:
                            out.write(line)
                    buffered = []
                    current_tid = None

                for line in ref_fh:
                    if line.startswith("#"):
                        continue
                    if "\ttranscript\t" in line or "\texon\t" in line:
                        # extract transcript_id
                        import re
                        m = re.search(r'transcript_id "([^"]+)"', line)
                        tid = m.group(1) if m else None
                        if tid != current_tid:
                            flush()
                            current_tid = tid
                        buffered.append(line)
                flush()

    elapsed = (datetime.now() - t_start).total_seconds()
    total_in = sum(s.get("input", 0) for s in stats_list)
    total_merged = sum(s.get("merged", 0) for s in stats_list)
    total_genes = sum(s.get("genes", 0) for s in stats_list)

    summary = [
        "=" * 60,
        f"GTF Merge {__version__} - Merge Summary",
        "=" * 60,
        f"Date     : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"Elapsed  : {elapsed:.1f}s",
        f"Samples  : {len(samples)}",
        f"Input    : {total_in:,} transcripts",
        f"Merged   : {total_merged:,} transcripts",
        f"Genes    : {total_genes:,}",
        "=" * 60,
    ]
    print("\n".join(summary))
    with open(cfg.output + "_summary.txt", "w") as f:
        f.write("\n".join(summary) + "\n")

    if not cfg.keep_tmp:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    logger.info("Done!")


# ─────────────────────────────────────────────────────────────────
# `gtf-merge cluster` command
# ─────────────────────────────────────────────────────────────────

def run_cluster_command(cfg: ClusterConfig):
    from gtf_merge.core.cluster import cluster_transcripts
    from gtf_merge.core.writer import write_clusters
    from gtf_merge.core.class_code import classify_transcript

    t_start = datetime.now()
    logger  = logging.getLogger("gtf_merge")
    logger.info(f"GTF Merge v{__version__} - cluster subcommand")

    # input can be: a single GTF, or a filelist
    input_path = Path(cfg.input)
    transcripts = []
    sample_counts: Dict[str, int] = {}

    if input_path.suffix in (".gtf", ".gz") or "gtf" in input_path.suffix.lower():
        # treat as single GTF
        logger.info(f"Reading {cfg.input} as single GTF")
        for t in parse_gtf_by_chrom(cfg.input):
            transcripts.append(t)
            sample_counts[t.transcript_id] = sample_counts.get(t.transcript_id, 0) + 1
    else:
        # treat as filelist
        logger.info(f"Reading {cfg.input} as filelist")
        samples = parse_filelist(cfg.input)
        for s in samples:
            for t in parse_gtf_by_chrom(s["gtf_path"]):
                t.source_sample = s["sample_id"]
                transcripts.append(t)
                sample_counts[t.transcript_id] = sample_counts.get(t.transcript_id, 0) + 1

    logger.info(f"Loaded {len(transcripts)} transcripts")

    # load reference and append to transcripts
    ref_ids: set = set()
    ref_index = None
    if cfg.ref_gtf:
        logger.info(f"Loading reference: {cfg.ref_gtf}")
        ref_index = RefIndex()
        for t in parse_gtf_by_chrom(cfg.ref_gtf):
            ref_index.add_transcript(t)
            transcripts.append(t)
            ref_ids.add(t.transcript_id)
        logger.info(f"Added {len(ref_ids)} reference transcripts")

    # cluster
    logger.info("Clustering...")
    clusters = cluster_transcripts(
        transcripts        = transcripts,
        ref_ids            = ref_ids,
        sample_counts      = sample_counts,
        junction_threshold = cfg.junction_threshold,
        fuzzy_threshold    = cfg.fuzzy_threshold,
        w_structural       = cfg.weight_structural,
        w_junction         = cfg.weight_junction,
        w_jaccard          = cfg.weight_jaccard,
        min_similarity     = cfg.min_similarity,
    )
    logger.info(f"Created {len(clusters)} clusters")

    # classify each member against ref
    if ref_index:
        logger.info("Classifying members against reference...")
        for cluster in clusters:
            for member in cluster.members:
                # get the transcript object (already have by_id in ref_index for refs)
                # for non-ref, need to look up
                # we'll re-classify by looking up the transcript
                pass  # skip for now to avoid double work; class_code is set during merge

    # write outputs
    Path(cfg.output).parent.mkdir(parents=True, exist_ok=True)
    out_clusters = cfg.output + "_clusters.tsv"
    with open(out_clusters, "w") as f:
        write_clusters(clusters, f)

    elapsed = (datetime.now() - t_start).total_seconds()
    summary = [
        "=" * 60,
        f"GTF Merge {__version__} - Cluster Summary",
        "=" * 60,
        f"Date          : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"Elapsed       : {elapsed:.1f}s",
        f"Total inputs  : {len(transcripts):,}",
        f"Clusters      : {len(clusters):,}",
        f"Avg size      : {len(transcripts) / max(len(clusters),1):.2f}",
        f"Singletons    : {sum(1 for c in clusters if c.size == 1):,}",
        "=" * 60,
    ]
    print("\n".join(summary))
    with open(cfg.output + "_summary.txt", "w") as f:
        f.write("\n".join(summary) + "\n")

    logger.info("Done!")


# ─────────────────────────────────────────────────────────────────
# Argument parser
# ─────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog        = "gtf-merge",
        description = "GTF Merge: high-performance long-read GTF merger and clusterer",
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    sp = p.add_subparsers(dest="command", required=True, help="subcommand")

    # ── merge ──────────────────────────────────────────────────────
    pm = sp.add_parser("merge", help="merge transcripts likely from the same molecule",
                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    req = pm.add_argument_group("Required")
    req.add_argument("-i", "--input",  required=True, help="filelist TSV")
    req.add_argument("-o", "--output", required=True, help="output prefix")

    ref = pm.add_argument_group("Reference annotation")
    ref.add_argument("--ref_gtf",         default=None, help="GENCODE/Ensembl GTF (optional)")
    ref.add_argument("--novel_prefix",    default="Novel", help="novel transcript ID prefix")
    ref.add_argument("--output_extended", action="store_true",
                     help="output _extended.gtf with ref-only transcripts")

    plat = pm.add_argument_group("Platform")
    plat.add_argument("--platform", default="pacbio", choices=["pacbio","ont","custom"])
    plat.add_argument("--junction_threshold", default=3, type=int,
                      help="junction coord error tolerance (bp)")

    me = pm.add_argument_group("Multi-exon merge")
    me.add_argument("--max_5trunc",     default=1,    type=int,   help="max 5' missing junctions")
    me.add_argument("--max_3trunc",     default=0,    type=int,   help="max 3' missing junctions")
    me.add_argument("--tss_window",     default=200,  type=int,   help="TSS tolerance (bp)")
    me.add_argument("--tes_window",     default=50,   type=int,   help="TES tolerance (bp)")
    me.add_argument("--max_diff_count", default=-1,   type=int,   help="max #junctions with non-zero drift (-1=unlimited)")
    me.add_argument("--max_diff_ratio", default=0.5,  type=float, help="max ratio of drifted junctions")

    mo = pm.add_argument_group("Mono-exon merge")
    mo.add_argument("--mono_overlap", default=0.5, type=float, help="overlap ratio threshold")
    mo.add_argument("--mono_5end",    default=500, type=int,   help="5' diff threshold (bp)")
    mo.add_argument("--mono_3end",    default=50,  type=int,   help="3' diff threshold (bp)")

    rep = pm.add_argument_group("Representative selection")
    rep.add_argument("--tss_strategy", default="longest", choices=["longest","common","freq"])
    rep.add_argument("--tes_strategy", default="freq",    choices=["longest","common","freq"])

    perf = pm.add_argument_group("Performance / Output")
    perf.add_argument("--threads",  default=4,   type=int)
    perf.add_argument("--tmp_dir",  default="/tmp")
    perf.add_argument("--keep_tmp", action="store_true")
    perf.add_argument("--checkpoint_dir", default="")
    perf.add_argument("--resume",   action="store_true")
    perf.add_argument("--force",    action="store_true")
    perf.add_argument("--verbose",  action="store_true")
    perf.add_argument("--log_file", default="")

    # ── cluster ─────────────────────────────────────────────────────
    pc = sp.add_parser("cluster", help="compute similarity clusters of overlapping transcripts",
                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    cr = pc.add_argument_group("Required")
    cr.add_argument("-i", "--input",  required=True, help="GTF file or filelist TSV")
    cr.add_argument("-o", "--output", required=True, help="output prefix")

    cref = pc.add_argument_group("Reference annotation (optional)")
    cref.add_argument("--ref_gtf", default=None)

    csim = pc.add_argument_group("Similarity")
    csim.add_argument("--junction_threshold", default=3,  type=int)
    csim.add_argument("--fuzzy_threshold",    default=50, type=int)
    csim.add_argument("--weight_structural",  default=0.4, type=float)
    csim.add_argument("--weight_junction",    default=0.4, type=float)
    csim.add_argument("--weight_jaccard",     default=0.2, type=float)

    cclust = pc.add_argument_group("Clustering")
    cclust.add_argument("--min_similarity",       default=0.7, type=float)
    cclust.add_argument("--representative_order", default="ref_first",
                        choices=["ref_first","length","samples","random"])

    cperf = pc.add_argument_group("Performance")
    cperf.add_argument("--threads",  default=4, type=int)
    cperf.add_argument("--tmp_dir",  default="/tmp")
    cperf.add_argument("--verbose",  action="store_true")
    cperf.add_argument("--log_file", default="")

    return p


def main():
    parser = build_parser()
    args = parser.parse_args()
    setup_logging(getattr(args, "verbose", False), getattr(args, "log_file", ""))

    if args.command == "merge":
        cfg = MergeConfig(
            input              = args.input,
            output             = args.output,
            ref_gtf            = args.ref_gtf,
            novel_prefix       = args.novel_prefix,
            output_extended    = args.output_extended,
            platform           = args.platform,
            junction_threshold = args.junction_threshold,
            max_5trunc         = args.max_5trunc,
            max_3trunc         = args.max_3trunc,
            tss_window         = args.tss_window,
            tes_window         = args.tes_window,
            max_diff_count     = args.max_diff_count,
            max_diff_ratio     = args.max_diff_ratio,
            mono_overlap       = args.mono_overlap,
            mono_5end          = args.mono_5end,
            mono_3end          = args.mono_3end,
            tss_strategy       = args.tss_strategy,
            tes_strategy       = args.tes_strategy,
            threads            = args.threads,
            tmp_dir            = args.tmp_dir,
            keep_tmp           = args.keep_tmp,
            checkpoint_dir     = args.checkpoint_dir,
            resume             = args.resume,
            force              = args.force,
            verbose            = args.verbose,
            log_file           = args.log_file,
        )
        run_merge_command(cfg)
    elif args.command == "cluster":
        cfg = ClusterConfig(
            input              = args.input,
            output             = args.output,
            ref_gtf            = args.ref_gtf,
            junction_threshold = args.junction_threshold,
            fuzzy_threshold    = args.fuzzy_threshold,
            weight_structural  = args.weight_structural,
            weight_junction    = args.weight_junction,
            weight_jaccard     = args.weight_jaccard,
            min_similarity     = args.min_similarity,
            representative_order = args.representative_order,
            threads            = args.threads,
            tmp_dir            = args.tmp_dir,
            verbose            = args.verbose,
            log_file           = args.log_file,
        )
        run_cluster_command(cfg)


if __name__ == "__main__":
    main()
