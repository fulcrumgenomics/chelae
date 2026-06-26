# chelae

<p align="center">
  <a href="https://github.com/fulcrumgenomics/chelae/actions?query=workflow%3ACheck"><img src="https://github.com/fulcrumgenomics/chelae/actions/workflows/build_and_test.yml/badge.svg" alt="Build Status"></a>
  <a href="https://github.com/fulcrumgenomics/chelae/blob/main/LICENSE"><img src="https://img.shields.io/github/license/fulcrumgenomics/chelae.svg" alt="license"></a>
  <a href="https://crates.io/crates/chelae"><img src="https://img.shields.io/crates/v/chelae.svg?colorB=319e8c" alt="Version info"></a>
  <a href="https://bioconda.github.io/recipes/chelae/README.html"><img src="https://img.shields.io/conda/vn/bioconda/chelae.svg?label=bioconda" alt="Bioconda"></a>
  <br>
</p>

A fast, accurate, multi-threaded toolkit for trimming and filtering short-read FASTQ data, written in Rust.

The name *chelae* is the plural of [*chela*](https://en.wikipedia.org/wiki/Chela_(organ)) — the pincer-like claws of crustaceans — a nod to what this tool does to FASTQ reads.

<p>
<a href="https://fulcrumgenomics.com">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset=".github/logos/fulcrumgenomics-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset=".github/logos/fulcrumgenomics-light.svg">
  <img alt="Fulcrum Genomics" src=".github/logos/fulcrumgenomics-light.svg" height="100">
</picture>
</a>
</p>

[Visit us at Fulcrum Genomics](https://www.fulcrumgenomics.com) to learn more about how we can power your Bioinformatics with chelae and beyond.

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-%2338b44a.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-%2326a8e0.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

*This README is user-facing documentation. Contributors working on `chelae` itself should see [CONTRIBUTING.md](CONTRIBUTING.md) for build conventions, the pre-push checks, and the release process.*

## Contents

- [Overview](#overview)
- [`chelae trim` — examples](#chelae-trim--examples)
- [`chelae trim` — options](#chelae-trim--options)
- [`chelae detect` — examples](#chelae-detect--examples)
- [`chelae detect` — options](#chelae-detect--options)
- [Performance](#performance)
- [Installing](#installing)
- [Build Targeting and Portability](#build-targeting-and-portability)

## Overview

`chelae` ships two subcommands:

- **`chelae trim`** — short-read FASTQ trim and filter pipeline (see below).
- **`chelae detect`** — identify the adapter sequence(s) present in a FASTQ. PE input discovers adapters via R1/R2 overlap; SE input scores reads against every built-in kit plus optional user candidates. Optional FASTA output is designed to feed straight back into `chelae trim --adapter-fasta`.

`chelae trim` performs common short-read preprocessing tasks in one pass, in the following order:

1. Poly-G 3' trim (on by default)
2. Adapter trimming — PE-overlap evidence mode, plus `--kit`, `--adapter-sequence`, and `--adapter-fasta` for SE or deep trimming
3. [Read-structure](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures) based hard-trim and UMI extraction (runs after adapter trim so tail-skip segments operate on the cleaned template)
4. Optional poly-X 3' trim (`--trim-polyx`)
5. Optional 5'→3' and/or 3'→5' sliding-window quality trim
6. Length filter post-trimming (`--filter-length MIN[:MAX]`)
7. Optional N-base filter, mean-quality filter, and low-quality-fraction filter

Outputs are BGZF-compressed FASTQ plus a fastp-compatible JSON report suitable for MultiQC.

`chelae` uses paired-read overlap detection combined with adapter-sequence confirmation to rapidly and confidently identify adapter sequence in paired-end reads.  This repository includes a benchmark suite in `benchmark-pipeline/`; `chelae` is the **fastest** tool tested across all experimental setups, while also providing the **highest accuracy** trimming. See the [Performance](#performance) section below.

## `chelae trim` — examples

PE-overlap adapter detection is on by default; for best accuracy it is also recommended to supply your adapter sequences via `--kit` or `--adapter-sequence`.

#### Trim paired-end reads with Illumina's truseq adapters
```bash
chelae trim \
    -i sample.r1.fq.gz sample.r2.fq.gz \
    -o trimmed.r1.fq.gz trimmed.r2.fq.gz \
    --kit truseq
```

#### Add in a 3' quality trim with an 8bp sliding window at Q20
```bash
chelae trim \
    -i sample.r1.fq.gz sample.r2.fq.gz \
    -o trimmed.r1.fq.gz trimmed.r2.fq.gz \
    --kit truseq \
    --quality-trim-3p 8:20
```

#### Trim paired-end reads while extracting an 8bp UMI and dropping 4 fixed bases in R1
```bash
chelae trim \
    -i sample.r1.fq.gz sample.r2.fq.gz \
    -o trimmed.r1.fq.gz trimmed.r2.fq.gz \
    --kit truseq \
    --read-structures 8M4S+T +T
```

## `chelae trim` — options

The tables below summarize every option accepted by `chelae trim`. For longer
explanations (rationale, units, edge cases) run `chelae trim --help`.

### Inputs, outputs, and runtime

| Option                          | Description                                                                                                  | Default |
|---------------------------------|--------------------------------------------------------------------------------------------------------------|---------|
| `-i, --inputs <PATHS>...`       | One (SE) or two (PE) FASTQ files; plain, gzip, or bgzf (auto-detected)                                       | —       |
| `-o, --outputs <PATHS>...`      | Output FASTQ path(s); count must match `--inputs`; always BGZF-compressed                                    | —       |
| `-t, --threads <N>`             | Number of threads to use                                                                                     | `4`     |
| `-c, --compression-level <1-12>`| BGZF compression level for output files                                                                      | `5`     |
| `-m, --metrics <PATH>`          | Optional path for the trimming metrics TSV; stdout summary is always emitted                                 | —       |
| `-j, --json <PATH>`             | Optional fastp-shape JSON report; consumed by MultiQC's `fastp` module unchanged                             | —       |

### Read-structure (hard-trim + UMI extraction)

| Option                                | Description                                                                                                | Default |
|---------------------------------------|------------------------------------------------------------------------------------------------------------|---------|
| `-r, --read-structures <RS>...`       | Optional [read-structures](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures) per input; supports `T` (template), `M` (UMI → read name), `S` (skip); applied after adapter trim | —       |
| `--discard-unsupported-segments`      | Treat `B` (sample barcode) and `C` (cellular barcode) segments as `S` (skip) instead of erroring          | off     |

### Adapter trimming

| Option                              | Description                                                                                                  | Default |
|-------------------------------------|--------------------------------------------------------------------------------------------------------------|---------|
| `-k, --kit <NAME>...`               | Built-in kit preset; repeatable. Known: `truseq`, `nextera`, `small-rna`, `aviti`, `mgi` (alias `dnbseq`), `all` | —       |
| `-a, --adapter-sequence <SEQ>...`   | 3' adapter sequence(s); 1 for SE, 1 or 2 for PE (R1, R2); ACGT or IUPAC                                      | —       |
| `-f, --adapter-fasta <PATH>`        | FASTA of adapter sequences; best match is trimmed                                                            | —       |
| `--adapter-min-length <N>`          | Minimum match length when searching the 3' end for an adapter sequence (SE mode; PE mode only for inserts < `overlap-min-length`) | `6`     |
| `--adapter-mismatch-rate <0..1>`    | Max fraction of mismatches when matching adapter against the 3' end (default ≈ 1 mismatch / 8 bases)         | `0.125` |

### Paired-end overlap detection

| Option                                | Description                                                                                                | Default |
|---------------------------------------|------------------------------------------------------------------------------------------------------------|---------|
| `--no-overlap-detection`              | Disable PE-overlap trim-point detection (ignored for SE); rely on sequence matching alone                  | on (PE) |
| `--overlap-min-length <N>`            | Minimum overlap (bp) required to declare R1/R2 overlap                                                     | `30`    |
| `--overlap-max-mismatch-rate <0..1>`  | Max fraction of mismatches in the overlap probe window                                                     | `0.10`  |
| `--overlap-diagnostic-length <N>`     | When evaluating PE overlap, only examine this many overlapping bases. Multiples of 16 ideal.               | `64`    |
| `--expected-insert-size <BP>`         | Hint for typical insert size; seeds the overlap candidate-walk order so the right overlap is found sooner  | —       |
| `--insert-size-stats`                 | Emit a fastp-shape per-pair insert-size histogram under `insert_size` in the JSON (extends overlap probing to I > R configurations) | off     |

### Poly-G / poly-X trimming

| Option                  | Description                                                                                                          | Default |
|-------------------------|----------------------------------------------------------------------------------------------------------------------|---------|
| `--trim-polyg <N>`      | 3' poly-G trim minimum run length; pass `0` to disable                                                               | `10`    |
| `--trim-polyx [<N>]`    | Enable 3' poly-X trim (A/C/T homopolymer tails, e.g. poly-A from RNA-seq) with the given minimum run length         | off     |

### Quality trimming

Both quality-trim modes shorten the read at the 3' end — the `-3p` / `-5p` suffix
indicates the scan direction, not the trim location. `-3p` is conservative
(keeps everything up to the last good window from the 3' end); `-5p` is
aggressive (cuts at the first bad window encountered from the 5' end).

| Option                          | Description                                                                                                  | Default      |
|---------------------------------|--------------------------------------------------------------------------------------------------------------|--------------|
| `--quality-trim-3p [<W:Q>]`     | Scan 3'→5'; trim trailing bases until a window of size `W` has mean quality ≥ `Q` (fastp `--cut_tail`)        | off (`8:20`) |
| `--quality-trim-5p [<W:Q>]`     | Scan 5'→3'; truncate at the first window of size `W` with mean quality < `Q` (fastp `--cut_right`)            | off (`8:20`) |

### Filters (applied after trimming; pair dropped if either mate fails)

| Option                            | Description                                                                                              | Default |
|-----------------------------------|----------------------------------------------------------------------------------------------------------|---------|
| `-l, --filter-length <MIN[:MAX]>` | Drop reads/pairs with post-trim length below `MIN` (or above `MAX`)                                      | `15`    |
| `--filter-max-ns <N>`             | Drop reads/pairs whose per-mate count of ambiguous (N) bases exceeds `N`                                 | off     |
| `--filter-mean-qual <Q>`          | Drop reads/pairs whose post-trim mean Phred quality is below `Q` (runs last, after every trim stage)     | off     |
| `--filter-low-qual <Q:F>`         | Drop reads/pairs where the fraction of bases below quality `Q` exceeds `F` (e.g. `15:0.4`)               | off     |

## `chelae detect` — examples

`chelae detect` samples reads from one or two FASTQ files and reports the
adapter sequence(s) present. PE input discovers adapters via R1/R2 overlap (no
kit knowledge required); SE input scores reads against every built-in kit plus
any user-supplied candidate. Discovered/winning adapter sequences can be
written as FASTA for direct re-use by `chelae trim --adapter-fasta`.

#### Discover the adapters in a paired-end library

```bash
chelae detect \
    -i sample.r1.fq.gz sample.r2.fq.gz \
    -o adapters.fa
```

#### Identify which built-in kit a single-end library is using

```bash
chelae detect -i sample.fq.gz
```

#### Score a single-end library against extra user-supplied candidates

```bash
chelae detect \
    -i sample.fq.gz \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -f custom-adapters.fa
```

#### Detect adapters, then trim using the discovered FASTA

```bash
chelae detect -i r1.fq.gz r2.fq.gz -o adapters.fa
chelae trim -i r1.fq.gz r2.fq.gz -o t.r1.fq.gz t.r2.fq.gz --adapter-fasta adapters.fa
```

## `chelae detect` — options

Most flags are PE- or SE-only; the help-text annotation in each row says which
mode the flag applies to. Run `chelae detect --help` for the full rationale.

| Option                                | Description                                                                                                | Default      |
|---------------------------------------|------------------------------------------------------------------------------------------------------------|--------------|
| `-i, --inputs <PATHS>...`             | One (SE) or two (PE) FASTQ files; plain, gzip, or bgzf (auto-detected)                                     | —            |
| `-o, --output-fasta <PATH>`           | Optional FASTA output of discovered/winning adapter(s); ready to feed back into `chelae trim --adapter-fasta` | —            |
| `-a, --adapter-sequence <SEQ>...`     | (SE only) Extra adapter candidate(s) to score against, in addition to every built-in kit                    | —            |
| `-f, --adapter-fasta <PATH>`          | (SE only) FASTA of extra adapter candidates; record names are preserved in the report                       | —            |
| `-n, --num-detections <N>`            | Target number of usable detections before stopping. Higher = more confident composition estimate            | `5000`       |
| `--max-reads <N>`                     | Hard cap on records scanned even if `--num-detections` isn't reached                                        | `1000000`    |
| `--min-detections-for-report <N>`     | Refuse to report if final detection count falls below this floor (avoids confident-looking tiny samples)    | `20`         |
| `--min-fraction <0..1>`               | Minimum share of detections an adapter must account for to be reported                                      | `0.05`       |
| `--min-tail-length <N>`               | Minimum length of adapter evidence (bp) per detection (PE post-template tail; SE matched alignment)         | `8`          |
| `--overlap-min-length <N>`            | (PE) Minimum overlap (bp) required for PE-overlap detection                                                 | `30`         |
| `--overlap-max-mismatch-rate <0..1>`  | (PE) Max fraction of mismatches in the overlap probe                                                        | `0.10`       |
| `--overlap-diagnostic-length <N>`     | (PE) Upper bound on the probe length per overlap-length candidate (bp)                                      | `64`         |
| `--adapter-min-length <N>`            | (SE) Minimum match length (bp) when scoring a candidate against a read's 3' end                             | `10`         |
| `--adapter-mismatch-rate <0..1>`      | (SE) Max fraction of mismatches when matching a candidate against a read's 3' end                           | `0.125`      |

## Performance

A Snakemake pipeline in [`benchmark-pipeline/`](benchmark-pipeline/) runs
`chelae` against eight other FASTQ trimmers on simulated short-read libraries
spanning cfDNA, WGS at several insert sizes, exome, miRNA, and a high-error
condition — scoring both runtime and adapter-trim accuracy against the
simulator's ground-truth boundaries.

**`chelae` is the fastest tool tested on every dataset, and posts the most
accurate adapter trimming on the majority of them.** The two tables below
show, per dataset, the top three tools by wall-clock time (`wgs` config:
adapter + quality + length filter) and the top three tools by adapter-trim
exact-match rate (`adapter_only` config). Every dataset is ~2× human WGS
worth of simulated reads (read counts vary with read length — roughly 21–84 M
pairs/reads per dataset). All numbers are 8-thread median over 3 replicates
on an EC2 `c8id.2xlarge` (Intel Xeon 6975P-C, Granite Rapids).

Tool name abbreviations used in the tables: **ar** = adapterremoval, **tg** =
trim-galore, **tg-rs** = trim-galore-rs, **tmatic** = trimmomatic.

### Runtime — top 3 fastest tools per dataset, wall seconds @ 8 threads

| ID | Layout | Insert | Err rate | #1              | #2              | #3              |
|---:|--------|------------:|---------:|-----------------|-----------------|-----------------|
| 1  | 2×150  |    150 ± 30 |   0.1–1% | **chelae** 37.4 | cutadapt 78.2   | tg-rs 78.2      |
| 2  | 2×150  |    250 ± 40 |   0.1–1% | **chelae** 41.9 | tg-rs 48.4      | cutadapt 68.0   |
| 3  | 2×150  |    350 ± 60 |   0.1–1% | **chelae** 42.2 | tg-rs 48.6      | cutadapt 65.8   |
| 4  | 2×150  |    450 ± 80 |   0.1–1% | **chelae** 42.0 | tg-rs 49.2      | cutadapt 65.7   |
| 5  | 2×250  |    450 ± 80 |   0.1–1% | **chelae** 37.9 | tg-rs 44.3      | cutadapt 57.5   |
| 6  | 2×150  |    250 ± 60 |     1–5% | **chelae** 20.3 | tg-rs 32.9      | fastp-nfcore 39.4 |
| 7  | 2×150  |    170 ± 30 |   0.1–1% | **chelae** 39.9 | tg-rs 60.3      | cutadapt 71.2   |
| 8  | 2×76   |    140 ± 25 |   0.1–1% | **chelae** 50.5 | tg-rs 61.5      | cutadapt 84.8   |
| 9  | 1×150  |    300 ± 80 |   0.1–1% | **chelae** 42.8 | cutadapt 68.8   | tg-rs 72.3      |
| 10 | 1×150  |    120 ± 30 |   0.1–1% | **chelae** 39.5 | bbduk 67.9      | fastp 77.0      |
| 11 | 1×76   |     30 ± 2  |   0.1–1% | **chelae** 38.1 | bbduk 43.3      | fastp 61.2      |

`chelae` is **#1 on every dataset**, with the runner-up trailing by 15–95 %.
The lead is widest on short-insert PE libraries where more trimming occurs and on the high error rate dataset.

### Accuracy — top 3 by adapter-trim RMSE per dataset (lower is better)

We score adapter-trim accuracy as **RMSE of the difference in trim point (in bases) vs. the simulator's ground truth**, aggregated across reads. RMSE is used so as to penalize large under- and over-trim errors more than small ones under the premise that each additional base lost (or adapter base retained) in a read is more consequential than the last.

| ID | Layout | Insert | Err rate | #1               | #2              | #3              |
|---:|--------|------------:|---------:|------------------|-----------------|-----------------|
| 1  | 2×150  |    150 ± 30 |   0.1–1% | **chelae** 0.045 | ar 0.079        | fastp 0.169     |
| 2  | 2×150  |    250 ± 40 |   0.1–1% | **chelae** 0.013 | cutadapt 0.189  | tg 0.189        |
| 3  | 2×150  |    350 ± 60 |   0.1–1% | **chelae** 0.009 | cutadapt 0.098  | tg 0.098        |
| 4  | 2×150  |    450 ± 80 |   0.1–1% | **chelae** 0.010 | cutadapt 0.090  | tg 0.090        |
| 5  | 2×250  |    450 ± 80 |   0.1–1% | **chelae** 0.011 | cutadapt 0.166  | tg 0.166        |
| 6  | 2×150  |    250 ± 60 |     1–5% | **chelae** 0.393 | ar 0.482        | cutadapt 1.310  |
| 7  | 2×150  |    170 ± 30 |   0.1–1% | **chelae** 0.040 | ar 0.139        | fastp 0.207     |
| 8  | 2×76   |    140 ± 25 |   0.1–1% | **chelae** 0.016 | cutadapt 0.284  | ar 0.312        |
| 9  | 1×150  |    300 ± 80 |   0.1–1% | **chelae** 0.267 | cutadapt 0.272  | fastp 0.284     |
| 10 | 1×150  |    120 ± 30 |   0.1–1% | fastp 0.652      | **chelae** 0.803 | cutadapt 0.964 |
| 11 | 1×76   |     30 ± 2  |   0.1–1% | tmatic 0.024        | **chelae** 0.062 | fastp 0.062    |

`chelae` wins on **every PE dataset and one SE dataset**, and comes second on
the other two SE cases.  On SE datasets, most tools perform similarly where it is possible to configure key parameters the same (i.e. min match length and maximum allowed error rate).

### Tool versions

All tools were at the latest version available on bioconda as of the
benchmarking date (2026-05-06), with one deliberate exception:
`fastp-nfcore` runs a second `fastp` environment pinned to the version
[nf-core/modules](https://github.com/nf-core/modules) currently ships, so
MultiQC-via-nf-core users can see what they would actually get downstream.

| Tool             | Version tested            |
|------------------|---------------------------|
| chelae           | 0.1.0                     |
| fastp            | 1.3.2                     |
| fastp-nfcore     | 1.1.0 (nf-core/modules pin) |
| cutadapt         | 5.2                       |
| trim-galore      | 0.6.11                    |
| trim-galore-rs   | 2.1.0 (Oxidized Edition)  |
| trimmomatic      | 0.40                      |
| bbduk (bbmap)    | 39.81                     |
| adapterremoval   | 2.3.4                     |

### More detail, raw data, and reproduction

For per-dataset tables with **every** tool (not just the top 3) and additional
metrics — RSS, user CPU, parallel efficiency, MAE, false-positive and
false-negative counts, etc. — see
[**benchmark-pipeline/RESULTS.md**](benchmark-pipeline/RESULTS.md).

The full per-row source data (one row per `(sample, trim_config, tool, threads,
replicate)`) lives in
`benchmark-pipeline/.benchmark-outputs/perf-20260508/results/{bench_summary,accuracy_summary}.tsv`.
See [`benchmark-pipeline/README.md`](benchmark-pipeline/README.md) for the
reproduction recipe.

## Installing

### Install from bioconda

Using [pixi](https://pixi.sh), after adding the `bioconda` channel:

```console
pixi add chelae
```

Or using your favorite conda client (`conda`, `mamba`, `micromamba`, …):

```console
conda install -c bioconda chelae
```

### Installing with `cargo`

To install with cargo you must first [install rust](https://doc.rust-lang.org/cargo/getting-started/installation.html). Which (on macOS and Linux) can be done with:

```console
curl https://sh.rustup.rs -sSf | sh
```

Then, to install `chelae` run:

```console
cargo install chelae
```

### Building From Source

First, clone the git repo:

```console
git clone https://github.com/fulcrumgenomics/chelae.git
```

If you do not already have rust development tools installed, install via [rustup](https://rustup.rs/):

```console
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Then build in release mode:

```console
cd chelae
cargo build --release
./target/release/chelae --help
```

## Build Targeting and Portability

x86_64 release binaries ship as a single `cargo multivers` launcher that embeds three CPU-specific builds and picks the best match at startup:

- `x86-64` — SSE2 baseline, runs on any 64-bit x86 CPU (2003+)
- `x86-64-v2` — SSE4.2 + POPCNT (2008+); captures nearly all of the historical "v3 wins 6%" codegen benefit
- `x86-64-v4` — AVX-512F/BW/CD/DQ/VL for Ice Lake / Sapphire Rapids / Granite Rapids / Zen 4+

The launcher is ~3.7 MB total and adds ~0.2 s of startup for decompression + `memfd_create + exec`. v3 is intentionally skipped — on chelae's workload v2 and v3 are within measurement noise, and v4 picks up what little additional win AVX-512 gives (~1% on our benchmarks).

aarch64 release binaries (Apple Silicon, AWS Graviton, GCP Axion, Azure Cobalt) are a single build with generic ARMv8-A / NEON baseline. Benchmarks showed Neoverse-specific tuning yields only ~1-2% over generic and cross-tuning penalty is near zero, so multivers isn't worth the complexity on aarch64.

For local development, `cargo build --release` uses `target-cpu=native` (see `.cargo/config.toml`) for fastest local runs.
