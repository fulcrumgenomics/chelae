<h1 align="center">
  <img src="docs/chelae-banner.png" alt="chelae: Ferris the crab clipping the adapter off the end of a paired-end read" width="800">
</h1>

<p align="center">
  <a href="https://github.com/fulcrumgenomics/chelae/actions?query=workflow%3ACheck"><img src="https://github.com/fulcrumgenomics/chelae/actions/workflows/build_and_test.yml/badge.svg" alt="Build Status"></a>
  <a href="https://github.com/fulcrumgenomics/chelae/blob/main/LICENSE"><img src="https://img.shields.io/github/license/fulcrumgenomics/chelae.svg" alt="license"></a>
  <a href="https://crates.io/crates/chelae"><img src="https://img.shields.io/crates/v/chelae.svg?colorB=319e8c" alt="Version info"></a>
  <a href="https://bioconda.github.io/recipes/chelae/README.html"><img src="https://img.shields.io/conda/vn/bioconda/chelae.svg?label=bioconda" alt="Bioconda"></a>
  <a href="https://doi.org/10.5281/zenodo.21445782"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.21445782.svg" alt="DOI"></a>
  <a href="https://www.fulcrumgenomics.com"><img src=".github/logos/fulcrumgenomics-badge.svg" alt="Fulcrum Genomics"></a>
  <br>
</p>

A fast, accurate, multi-threaded toolkit for trimming and filtering short-read FASTQ data, written in Rust. Its name is the plural of [*chela*](https://en.wikipedia.org/wiki/Chela_(organ)), the pincer-like claws of crustaceans: a nod to what it does to FASTQ reads.

## Contents

- [Overview](#overview)
- [Examples](#examples)
- [Performance](#performance)
- [Installing](#installing)
- [About Fulcrum Genomics](#about-fulcrum-genomics)

Every option, the input and output rules, and more examples are in [docs/usage.md](docs/usage.md).

## Overview

`chelae` trims and filters short-read FASTQ in one multi-threaded pass, and identifies the adapters in a library you know nothing about.

- **Accurate paired-end adapter trimming.** `chelae` finds each pair's insert from where R1 and R2 overlap, then checks the implied adapter against known adapter sequences. That catches adapters too short for sequence matching alone to find, without mistaking adapter-like sequence inside the insert for adapter.
- **Fast.** SIMD kernels and a pipeline built for many cores trimmed about 1.7 M read pairs per second on 8 cores in our benchmark.
- **Benchmarked.** Against six other trimmers, `chelae` was the fastest, 1.25× faster than the runner-up and 2.6–5.7× faster than cutadapt, trim-galore-rs and fastp, and the most accurate on 8 of 11 simulated datasets. See [Performance](#performance).
- **One pass does it all:** poly-G, adapter, [read-structure](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures) hard-trimming with UMI extraction, poly-X and quality trimming, then length, N-base and quality filters.
- **Adapter detection.** `chelae detect` reports the adapters in a library and writes them as FASTA for `chelae trim --adapter-fasta`.
- **Fits into pipelines.** Split or interleaved paired-end files, stdin and stdout, gzip or BGZF input detected automatically, and a fastp-compatible JSON report for MultiQC.

Every option, the input and output rules, and more examples are in [docs/usage.md](docs/usage.md).

## Examples

Paired-end reads can be trimmed without naming any adapters. `chelae` finds adapters from where R1 and R2 overlap, and checks each one against every adapter in its built-in database (TruSeq, Nextera, small RNA, AVITI and MGI/DNBSEQ):

```bash
chelae trim \
    -i sample.r1.fq.gz sample.r2.fq.gz \
    -o trimmed.r1.fq.gz trimmed.r2.fq.gz
```

When you know the kit, name it - specify the right kit will make chelae a little faster _and_ a little more accurate (since it can't match to the wrong kits)  The following example also moves an 8 bp UMI from the start of R1 into the read name, skips the next 4 bases, and quality-trims 3' ends:

```bash
chelae trim \
    -i sample.r1.fq.gz sample.r2.fq.gz \
    -o trimmed.r1.fq.gz trimmed.r2.fq.gz \
    --kit truseq \
    --read-structures 8M4S+T +T \
    --quality-trim-3p 8:20
```

`-i` and `-o` default to stdin and stdout, and interleaved paired-end input is recognized from its read names, so `chelae` can sit in a pipeline with no intermediate files:

```bash
samtools fastq sample.unmapped.bam \
    | chelae trim --kit truseq \
    | bwa mem -p ref.fa - \
    | samtools sort -o sample.bam
```

For a library of unknown provenance, `chelae detect` finds the adapters, and its FASTA feeds straight back into `chelae trim`:

```bash
chelae detect -i sample.r1.fq.gz sample.r2.fq.gz -o adapters.fa
chelae trim \
    -i sample.r1.fq.gz sample.r2.fq.gz \
    -o trimmed.r1.fq.gz trimmed.r2.fq.gz \
    --adapter-fasta adapters.fa
```

For single-end reads, `chelae detect` reports which built-in kit matches best:

```bash
chelae detect -i sample.fq.gz
```

## Performance

`chelae` was benchmarked against six other FASTQ trimmers on simulated short-read libraries, for runtime and for adapter-trim accuracy against the simulator's ground truth. [RESULTS.md](benchmark-pipeline/RESULTS.md) has the method, the [dataset definitions](benchmark-pipeline/RESULTS.md#datasets) and every tool's numbers.

### Runtime

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://github.com/fulcrumgenomics/chelae/raw/HEAD/docs/throughput-dark.svg?sanitize=true">
  <source media="(prefers-color-scheme: light)" srcset="https://github.com/fulcrumgenomics/chelae/raw/HEAD/docs/throughput-light.svg?sanitize=true">
  <img alt="Throughput at 8 threads, million read pairs per second: chelae 1.71, adapterremoval 1.36, cutadapt 0.63, trim-galore-rs 0.58, fastp 0.34" src="https://github.com/fulcrumgenomics/chelae/raw/HEAD/docs/throughput-light.svg?sanitize=true" width="720">
</picture>

<details>
<summary>**Wall times, setup, CPU time and memory**</summary>

Wall seconds for 50 M read pairs at 8 threads, median of 3 runs ([WGS](benchmark-pipeline/RESULTS.md#dataset-t1---pe-wgs-350bp-2x150-50m), [cfDNA](benchmark-pipeline/RESULTS.md#dataset-t2---pe-cfdna-170bp-2x150-50m)):

| Tool           | WGS, insert 350 ± 60 | cfDNA, insert 170 ± 30 |
|----------------|---------------------:|-----------------------:|
| **chelae**     |             **29.7** |               **29.0** |
| adapterremoval |         37.4 (1.26×) |           36.1 (1.25×) |
| cutadapt       |         77.6 (2.61×) |           81.9 (2.83×) |
| trim-galore-rs |         80.0 (2.69×) |           93.9 (3.24×) |
| fastp          |        169.6 (5.71×) |          132.5 (4.58×) |

- Both libraries are simulated 2×150 paired-end reads with TruSeq adapters, trimmed with a realistic configuration: adapter, poly-G, sliding-window quality and length trimming, and an N filter.
- The host was an EC2 `r8a.2xlarge` (AMD EPYC 9R45, Zen 5, 8 cores). Inputs and outputs were on a RAM disk so the timings measure the tools rather than storage, and every tool with a compression option wrote gzip at level 1.
- The three runs of each tool agreed to within 6%.
- bbduk and trimmomatic weren't timed at this scale: on the paired-end accuracy datasets below, they took 8–16× and 4–27× as long as `chelae` (see [broad toolset timings](benchmark-pipeline/RESULTS.md#tier-1-accuracy-and-screen)).

| Tool           | User CPU s, WGS | User CPU s, cfDNA |  Max RSS |
|----------------|----------------:|------------------:|---------:|
| **chelae**     |             207 |               202 |   117 MB |
| adapterremoval |             257 |               237 |    55 MB |
| cutadapt       |             358 |               418 |    76 MB |
| trim-galore-rs |             618 |               731 |   140 MB |
| fastp          |           1,251 |               974 | 1,274 MB |

</details>

### Accuracy

Eight paired-end and three single-end datasets of 1–4 M reads each simulate common library types: WGS at several insert sizes, cfDNA, exome with Nextera adapters, miRNA, and a high error rate ([definitions](benchmark-pipeline/RESULTS.md#datasets)). Accuracy is the [RMSE of each read's trim point](benchmark-pipeline/RESULTS.md#scoring-accuracy) against the truth, in bases (lower is better).

| Tool           | Most accurate on | Median RMSE, paired-end | Median RMSE, single-end |
|----------------|-----------------:|------------------------:|------------------------:|
| **chelae**     |      **8 of 11** |               **0.013** |               **0.265** |
| adapterremoval |          2 of 11 |                   0.100 |                   0.412 |
| cutadapt       |          0 of 11 |                   0.237 |                   0.270 |
| trim-galore-rs |          0 of 11 |                   0.237 |                   0.270 |
| fastp          |          1 of 11 |                   1.687 |                   0.280 |
| trimmomatic    |          0 of 11 |                   1.063 |                   1.172 |
| bbduk          |          0 of 11 |                   1.599 |                   0.902 |

<details>
<summary>Top 3 per dataset</summary>

**ar** = adapterremoval, **tg-rs** = trim-galore-rs, **tmatic** = trimmomatic. Each ID links to that dataset's full results.

| ID | Layout | Insert   | Err rate | #1               | #2               | #3               |
|---:|--------|---------:|---------:|------------------|------------------|------------------|
| [1](benchmark-pipeline/RESULTS.md#dataset-1---pe-wgs-150bp-2x150) | 2×150  | 150 ± 30 |   0.1–1% | **chelae** 0.039 | ar 0.050         | fastp 0.169      |
| [2](benchmark-pipeline/RESULTS.md#dataset-2---pe-wgs-250bp-2x150) | 2×150  | 250 ± 40 |   0.1–1% | **chelae** 0.013 | ar 0.113         | cutadapt 0.193   |
| [3](benchmark-pipeline/RESULTS.md#dataset-3---pe-wgs-350bp-2x150) | 2×150  | 350 ± 60 |   0.1–1% | **chelae** 0.008 | cutadapt 0.099   | tg-rs 0.099      |
| [4](benchmark-pipeline/RESULTS.md#dataset-4---pe-wgs-450bp-2x150) | 2×150  | 450 ± 80 |   0.1–1% | **chelae** 0.011 | cutadapt 0.089   | tg-rs 0.089      |
| [5](benchmark-pipeline/RESULTS.md#dataset-5---pe-wgs-450bp-2x250) | 2×250  | 450 ± 80 |   0.1–1% | **chelae** 0.010 | ar 0.130         | cutadapt 0.165   |
| [6](benchmark-pipeline/RESULTS.md#dataset-6---pe-wgs-higherr-250bp-2x150) | 2×150  | 250 ± 60 |     1–5% | ar 0.095         | **chelae** 0.417 | cutadapt 1.314   |
| [7](benchmark-pipeline/RESULTS.md#dataset-7---pe-cfdna-170bp-2x150) | 2×150  | 170 ± 30 |   0.1–1% | **chelae** 0.041 | ar 0.077         | fastp 0.220      |
| [8](benchmark-pipeline/RESULTS.md#dataset-8---pe-exome-140bp-2x76-nextera) | 2×76   | 140 ± 25 |   0.1–1% | **chelae** 0.013 | ar 0.074         | cutadapt 0.280   |
| [9](benchmark-pipeline/RESULTS.md#dataset-9---se-wgs-300bp-1x150) | 1×150  | 300 ± 80 |   0.1–1% | **chelae** 0.265 | cutadapt 0.270   | tg-rs 0.270      |
| [10](benchmark-pipeline/RESULTS.md#dataset-10---se-short-120bp-1x150) | 1×150  | 120 ± 30 |   0.1–1% | fastp 0.650      | ar 0.742         | **chelae** 0.801 |
| [11](benchmark-pipeline/RESULTS.md#dataset-11---se-mirna-30bp-1x76) | 1×76   |   30 ± 2 |   0.1–1% | ar 0.000         | tmatic 0.040     | **chelae** 0.065 |

Dataset 8 is an exome library with Nextera adapters and dataset 11 is miRNA; the rest are WGS or cfDNA with TruSeq adapters. On dataset 11 the top six tools differ by at most 5 of 4.2 M reads.

</details>

[Versions](benchmark-pipeline/RESULTS.md#tools): `chelae` 0.2.0, adapterremoval 3.0.2, bbduk 40.02, cutadapt 5.2, fastp 1.3.7, trim-galore-rs 2.3.0 and trimmomatic 0.41, the latest on bioconda as of 2026-09-25.

## Installing

### From bioconda

Using [pixi](https://pixi.sh), after adding the `bioconda` channel:

```console
pixi add chelae
```

Or using your favorite conda client (`conda`, `mamba`, `micromamba`, …):

```console
conda install -c bioconda chelae
```

### With `cargo`

With [Rust](https://doc.rust-lang.org/cargo/getting-started/installation.html) 1.89 or newer installed:

```console
cargo install chelae
```

To build from source, see [CONTRIBUTING.md](CONTRIBUTING.md).

## About Fulcrum Genomics

[Visit us at Fulcrum Genomics](https://www.fulcrumgenomics.com) to learn more about how we can power your Bioinformatics with chelae and beyond.

<p>
<a href="https://www.fulcrumgenomics.com"><picture><source media="(prefers-color-scheme: dark)" srcset="https://github.com/fulcrumgenomics/chelae/raw/HEAD/.github/logos/fulcrumgenomics-dark.svg?sanitize=true"><source media="(prefers-color-scheme: light)" srcset="https://github.com/fulcrumgenomics/chelae/raw/HEAD/.github/logos/fulcrumgenomics-light.svg?sanitize=true"><img alt="Fulcrum Genomics" src="https://github.com/fulcrumgenomics/chelae/raw/HEAD/.github/logos/fulcrumgenomics-light.svg?sanitize=true" height="36" align="middle"></picture></a>&nbsp;&nbsp;&nbsp;
<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-%2338b44a.svg?&style=for-the-badge&logo=gmail&logoColor=white" alt="Email us" align="middle"></a>&nbsp;&nbsp;
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-%2326a8e0.svg?&style=for-the-badge&logo=wordpress&logoColor=white" alt="Visit us" align="middle"></a>
</p>
