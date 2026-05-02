# GTF Merge

**High-performance GTF merger for long-read RNA-seq data (PacBio Iso-seq / Oxford Nanopore)**

Designed for merging hundreds of samples while preserving complete provenance tracking of every transcript.

---

## Features

- **Two-step merge strategy**: exact junction match → fuzzy graph-based clustering
- **Strand-aware truncation handling**: asymmetric 5'/3' truncation thresholds
- **Mono-exon aware**: separate merging logic with 3'-end anchoring (polyA-enriched libraries)
- **Full provenance tracking**: per-transcript, per-junction, per-sample detail
- **Reference annotation integration**: compare against GENCODE/Ensembl, inherit known IDs
- **Checkpoint/resume**: safely interrupt and resume 500-sample runs
- **Multi-threaded**: parallel processing by chromosome
- **Memory efficient**: streaming GTF parser, no full-genome load

---

## Installation

```bash
# From PyPI (after release)
pip install gtf-merge

# From GitHub (recommended)
git clone https://github.com/hmutpw/gtfMerge
cd gtf-merge
pip install -e .

# With conda
conda create -n gtf-merge python=3.9
conda activate gtf-merge
pip install -e .
```

**Requirements**: Python ≥ 3.8, `intervaltree`

---

## Quick Start

```bash
# 1. Prepare input file list (filelist.tsv, tab-separated, header required)
# gtf_path       sample_id   cap_status  priority  metadata
# /path/s1.gtf   SAMPLE1     no_cap      1         tissue:brain
# /path/s2.gtf   SAMPLE2     no_cap      1         tissue:liver

# 2. Run merge
gtf-merge \
  -i filelist.tsv \
  -o results/merged \
  --ref_gtf gencode.vM35.annotation.gtf \
  --platform pacbio \
  --threads 8

# 3. Outputs
# results/merged.gtf             - merged annotation
# results/merged_tracking.tsv    - full provenance
# results/merged_sample_stats.tsv
# results/merged_summary.txt
```

---

## Toy Example

A complete working example is provided in the `example/` directory.

### Example data

```
example/
├── sample1_brain.gtf        # brain tissue, full-length transcripts
├── sample2_liver.gtf        # liver tissue, 5-prime truncations + fuzzy junctions
├── sample3_ESC.gtf          # ESC, novel isoforms + TE mono-exon transcripts
├── reference_annotation.gtf # mini GENCODE-like reference (3 known genes)
└── filelist.tsv             # input file list
```

The three samples cover the main merge scenarios:

| Scenario | Gene | Description |
|----------|------|-------------|
| **Exact match** | *Actb*, *Gapdh*, *Tnf* | Same junction chain across all samples |
| **5-prime truncation** | *Actb* liver | liver transcript missing first junction |
| **Fuzzy junction** | *Gapdh* liver | 1 bp coordinate drift in exon boundary |
| **3-prime truncation** | *Tnf* liver | liver transcript missing last junction |
| **Novel isoform** | *Actb* ESC | alternative 3-prime end, partial junction match |
| **Mono-exon overlap** | TE transcripts | 3 samples, slightly different 5-prime starts |
| **Singleton** | Novel gene | only in liver, unique structure kept |

### Run the example

```bash
cd example/

gtf-merge \
  -i filelist.tsv \
  -o output/merged \
  --ref_gtf reference_annotation.gtf \
  --platform pacbio \
  --threads 1 \
  --verbose
```

### Expected summary

```
============================================================
GTF Merge Summary
============================================================
Input samples        : 3
Input transcripts    : 16

--- Multi-exon ---
  Exact merged       : 8
  Fuzzy merged       : 5

--- Mono-exon ---
  Merged             : 2

Total output transcripts : 15
Total genes              : 6
============================================================
```

### Inspect results

```bash
# count transcripts per merge type
grep -P "\ttranscript\t" output/merged.gtf | \
  grep -oP 'merge_type "\K[^"]+' | sort | uniq -c

# which transcripts matched GENCODE?
grep -P "\ttranscript\t" output/merged.gtf | \
  grep -oP 'ref_match_type "\K[^"]+' | sort | uniq -c

# tracking: transcripts supported by 2+ samples
awk -F'\t' 'NR>1 && $5>=2 {print $1, $3, $5}' output/merged_tracking.tsv
```

```r
# R: post-merge filtering using tracking file
library(data.table)
tracking <- fread("output/merged_tracking.tsv")

# transcripts present in all 3 samples
tracking[num_samples == 3]

# novel transcripts not in reference
tracking[ref_match_type == "novel"]

# TE-derived mono-exon transcripts
tracking[merge_type %in% c("monoexon_exact", "monoexon_overlap")]

# high-confidence for reference annotation building
tracking[num_samples >= 2 & merge_type %in% c("exact_match", "truncation_5", "truncation_3")]
```

---

## Testing

GTF Merge includes a test suite covering core functionality.

### Run all tests

```bash
cd gtf-merge
pip install pytest
python -m pytest tests/ -v
```

### Expected output

```
tests/test_core.py::TestSimilarity::test_identical_junction_chains PASSED
tests/test_core.py::TestSimilarity::test_truncation_5prime         PASSED
tests/test_core.py::TestSimilarity::test_fuzzy_junction            PASSED
tests/test_core.py::TestSimilarity::test_no_match                  PASSED
tests/test_core.py::TestSimilarity::test_monoexon_exact            PASSED
tests/test_core.py::TestSimilarity::test_monoexon_overlap          PASSED
tests/test_core.py::TestSimilarity::test_monoexon_no_overlap       PASSED
tests/test_core.py::TestExactMerge::test_basic_exact_merge         PASSED
tests/test_core.py::TestParser::test_parse_attributes              PASSED

9 passed in 0.08s
```

### What is tested

| Test class | What it covers |
|-----------|----------------|
| `TestSimilarity` | Junction matching, truncation detection, coordinate error scoring, mono-exon overlap logic |
| `TestExactMerge` | Exact merge grouping, sample deduplication, merge_type assignment |
| `TestParser` | GTF attribute string parsing |

---

## Input File List Format

Tab-separated, first line is header:

| Column | Required | Description |
|--------|----------|-------------|
| `gtf_path` | ✅ | Full path to sample GTF |
| `sample_id` | ✅ | Unique sample identifier |
| `cap_status` | ✅ | `capped` or `no_cap` |
| `priority` | optional | Integer, lower = higher priority (default: 1) |
| `metadata` | optional | `key:value` pairs, comma-separated |

```tsv
gtf_path	sample_id	cap_status	priority	metadata
/data/s1.gtf	ENCBS004AHE	no_cap	1	tissue:brain
/data/s2.gtf	ENCBS007ZPS	no_cap	1	tissue:liver
/data/s3.gtf	ENCBS062SUO	capped	1	tissue:brain
```

---

## Merge Strategy

### Step 1: Exact Merge
Transcripts with **identical junction chains** are merged immediately. Only TSS/TES differences are resolved.

### Step 2: Fuzzy Merge (Graph Clustering)
Remaining transcripts are compared pairwise. Two transcripts are connected if:
- `similarity_score ≥ --min_similarity`
- Truncation is 5'-only or 3'-only (not internal)

Connected components form clusters → merged into representative transcripts.

**Similarity score**:
```
score = w1 × junction_match_score + w2 × coord_accuracy_score

junction_match_score = matched_junctions / shorter_transcript_junctions
coord_accuracy_score = mean(1 - coord_error / threshold) per matched junction
```

### Step 3: Mono-exon Merge
Single-exon transcripts merged by overlap ratio + 3'/5' end proximity.

### Step 4: GENCODE Comparison
Merged transcripts compared to reference:
- `full_match` → inherit GENCODE transcript_id
- `junction_match` → inherit gene_id, new transcript_id
- `partial_match` / `gene_match` → inherit gene_id
- `novel` → assigned novel ID

---

## Parameters

### Platform
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--platform` | `pacbio` | `pacbio` (3bp threshold) / `ont` (5bp) / `custom` |
| `--junction_thresh` | `3` | Junction coordinate error tolerance (bp) |

### Multi-exon Merging
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_5trunc` | `-1` | Max 5' missing junctions (-1 = unlimited) |
| `--max_3trunc` | `1` | Max 3' missing junctions |
| `--min_similarity` | `0.9` | Minimum similarity for fuzzy merge |
| `--similarity_w1` | `0.7` | Weight for junction match component |
| `--similarity_w2` | `0.3` | Weight for coordinate accuracy component |
| `--require_canonical` | `False` | Only keep GT-AG/GC-AG junctions |

### TSS/TES Strategy
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--tss_thresh` | `200` | TSS tolerance window (bp) |
| `--tes_thresh` | `50` | TES tolerance window (bp) |
| `--tss_strategy` | `longest` | `longest` / `common` / `freq` |
| `--tes_strategy` | `freq` | `longest` / `common` / `freq` |

> `longest`: pick the most upstream/downstream coordinate (best for PacBio with 5' truncation)  
> `freq`: pick the most frequent coordinate across samples (best for reliable 3' ends)

### Mono-exon
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--mono_overlap` | `0.5` | Min overlap ratio for merging |
| `--mono_3end` | `50` | Max 3' end difference (bp) |
| `--mono_5end` | `500` | Max 5' end difference (bp) |

### Filtering
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_samples` | `1` | Min supporting samples (applied post-merge) |
| `--min_transcript_len` | `100` | Min transcript length (bp) |

### Performance
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--threads` | `4` | Parallel threads (per chromosome) |
| `--tmp_dir` | `/tmp` | Temporary file directory |

### Resume
| Parameter | Description |
|-----------|-------------|
| `--resume` | Resume from checkpoint |
| `--force` | Ignore checkpoint, rerun all |
| `--checkpoint_dir` | Custom checkpoint directory |

---

## Output Files

### `merged.gtf`
Standard GTF with extended attributes on transcript lines:

```
merge_type       "exact_match"          # merge category
similarity       "1.0000"               # similarity score
num_samples      "5"                    # supporting samples
sample_ids       "S1,S2,S3,S4,S5"
freq_score       "0.8500"               # fraction of all samples
ref_match_type   "junction_match"       # vs GENCODE
ref_gene_id      "ENSMUSG00000051951"
ref_transcript_id "ENSMUST00000159265"
tss_diff_range   "0-150bp"
tes_diff_range   "0-20bp"
truncation_5_count "2"
truncation_5_ratio "0.4000"
junction_coords  "100-200,500-700"
anchor_junctions "anchor,anchor"
per_junction_support "j0:S1:S2:S3,j1:S1:S2"
```

### `merged_tracking.tsv`
35-column TSV with complete per-sample, per-junction provenance. Key columns:

| Column | Description |
|--------|-------------|
| `merged_id` | New transcript ID |
| `merge_type` | Merge category |
| `similarity_score` | Weighted similarity |
| `num_samples` | Supporting samples |
| `per_sample_merge_type` | How each sample was merged |
| `junction_coords` | All junction coordinates |
| `per_junction_sample_support` | Which samples support each junction |
| `anchor_junctions` | anchor/variable per junction |
| `per_junction_max_coord_diff` | Max coordinate error per junction |
| `truncation_5_count` / `truncation_3_count` | Truncation counts |
| `missing_junctions_5` / `_3` | Which samples lost which junctions |
| `ref_match_type` | GENCODE comparison result |

### Merge Types

| Type | Meaning | Confidence |
|------|---------|-----------|
| `exact_match` | Identical junction chains | ★★★★★ |
| `truncation_5` | 5' junctions missing | ★★★★ |
| `truncation_3` | 3' junctions missing | ★★★★ |
| `truncation_both` | Both ends truncated | ★★★ |
| `fuzzy_junction` | Junction coordinate errors | ★★★ |
| `fuzzy_truncation_5/3` | Truncation + errors | ★★★ |
| `monoexon_exact` | Mono-exon, same coordinates | ★★★★ |
| `monoexon_overlap` | Mono-exon, overlapping | ★★ |
| `singleton` | Unique structure, no merge | N/A |

---

## Recommended Parameters

### PacBio Iso-seq (polyA+, multiple tissues)
```bash
gtf-merge \
  -i filelist.tsv \
  -o merged \
  --ref_gtf gencode.vM35.annotation.gtf \
  --platform pacbio \
  --junction_thresh 3 \
  --max_5trunc -1 \
  --max_3trunc 1 \
  --tss_strategy longest \
  --tes_strategy freq \
  --tss_thresh 200 \
  --tes_thresh 50 \
  --mono_3end 50 \
  --min_similarity 0.9 \
  --threads 8
```

### Oxford Nanopore
```bash
gtf-merge \
  -i filelist.tsv \
  -o merged \
  --platform ont \
  --junction_thresh 5 \
  --tss_thresh 300 \
  --threads 8
```

---

## Downstream Analysis

Post-merge filtering with R:
```r
library(data.table)
tracking <- fread("merged_tracking.tsv")

# Keep transcripts present in ≥1% of samples (e.g., 5 out of 500)
high_conf <- tracking[num_samples >= 5]

# Keep only exact/truncation merges
reliable <- tracking[merge_type %in% c("exact_match","truncation_5","truncation_3")]

# TE-related: lower threshold
te_transcripts <- tracking[num_samples >= 2 & ref_match_type == "novel"]
```

---

## Validate Inputs (dry run)
```bash
gtf-merge -i filelist.tsv -o out --validate_only
```

---

## Running on HPC (SLURM)

```bash
#!/usr/bin/env bash
#SBATCH --job-name=gtf_merge
#SBATCH --cpus-per-task=16
#SBATCH --mem=120g
#SBATCH --time=2-00:00:00

gtf-merge \
  -i filelist.tsv \
  -o /scratch/merged/mm39_merged \
  --ref_gtf gencode.vM35.annotation.gtf \
  --platform pacbio \
  --threads 16 \
  --checkpoint_dir /scratch/merged/checkpoint \
  --resume \
  --log_file gtf_merge.log
```

---

## Citation

If you use GTF Merge in your research, please cite:
> [manuscript in preparation]

---

## License

MIT License
