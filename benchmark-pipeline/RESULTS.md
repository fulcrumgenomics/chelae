# Benchmark results - September 25th 2026

Companion to the [Performance](../README.md#performance) section in the project README, with every tool's numbers from the same runs.

## Table of contents

- [Run metadata](#run-metadata)
- [Tools](#tools)
- [Scoring accuracy](#scoring-accuracy)
- [Datasets](#datasets)
- [Tier 1: accuracy and screen](#tier-1-accuracy-and-screen)
  - [Dataset 1 - pe-wgs-150bp-2x150](#dataset-1---pe-wgs-150bp-2x150)
  - [Dataset 2 - pe-wgs-250bp-2x150](#dataset-2---pe-wgs-250bp-2x150)
  - [Dataset 3 - pe-wgs-350bp-2x150](#dataset-3---pe-wgs-350bp-2x150)
  - [Dataset 4 - pe-wgs-450bp-2x150](#dataset-4---pe-wgs-450bp-2x150)
  - [Dataset 5 - pe-wgs-450bp-2x250](#dataset-5---pe-wgs-450bp-2x250)
  - [Dataset 6 - pe-wgs-higherr-250bp-2x150](#dataset-6---pe-wgs-higherr-250bp-2x150)
  - [Dataset 7 - pe-cfdna-170bp-2x150](#dataset-7---pe-cfdna-170bp-2x150)
  - [Dataset 8 - pe-exome-140bp-2x76-nextera](#dataset-8---pe-exome-140bp-2x76-nextera)
  - [Dataset 9 - se-wgs-300bp-1x150](#dataset-9---se-wgs-300bp-1x150)
  - [Dataset 10 - se-short-120bp-1x150](#dataset-10---se-short-120bp-1x150)
  - [Dataset 11 - se-mirna-30bp-1x76](#dataset-11---se-mirna-30bp-1x76)
- [Tier 2: runtime](#tier-2-runtime)
  - [Dataset T1 - pe-wgs-350bp-2x150-50m](#dataset-t1---pe-wgs-350bp-2x150-50m)
  - [Dataset T2 - pe-cfdna-170bp-2x150-50m](#dataset-t2---pe-cfdna-170bp-2x150-50m)

## Run metadata

- **Date**: 2026-09-25
- **Host**: EC2 r8a.2xlarge - AMD EPYC 9R45 (Zen 5), 8 cores (no SMT) sharing one 32 MiB L3, 64 GiB RAM
- **OS**: Ubuntu 24.04 (kernel 7.0)
- **chelae**: 0.2.0, the x86-64 cargo-multivers launcher (which runs its x86-64-v4 build on this CPU)
- **Threads**: 8 for every tool
- **Output compression**: gzip level 1 for every tool with the option; trimmomatic has none and uses its default
- **Storage**: each dataset's inputs were copied to a tmpfs before its runs, and every tool wrote its output there, so no timed run touched the disk
- **Replicates**: 1 per cell in tier 1; 3 per cell in tier 2
- **Data**: simulated with holodeck 0.2.1 from GRCh38 (1000 Genomes analysis set with decoys and HLA)

The per-run data behind every table is in [`published/2026-09-25/`](published/2026-09-25/).

## Tools

Every tool was at the latest version on bioconda as of 2026-09-25. `tools/<name>/render.py` holds the exact command line each tool was run with.

| Tool             | Version                   |
|------------------|---------------------------|
| chelae           | 0.2.0                     |
| adapterremoval   | 3.0.2                     |
| bbduk (bbmap)    | 40.02                     |
| cutadapt         | 5.2                       |
| fastp            | 1.3.7                     |
| trim-galore-rs   | 2.3.0 (Oxidized Edition)  |
| trimmomatic      | 0.41                      |

The Perl Trim Galore (0.6.x) is not included: it runs cutadapt, so its accuracy is cutadapt's, and it is far slower than its Rust rewrite. fastp is no longer benchmarked a second time at the version nf-core/modules pins, since that is now 1.3.6.

## Scoring accuracy

Accuracy is scored as the **RMSE of the difference in trim point (in bases) vs. the simulator's ground truth**, aggregated across reads. RMSE penalizes large under- and over-trim errors more than small ones, on the premise that each additional base lost (or adapter base retained) in a read is more consequential than the last.

> **Naming note**: dataset names are written with hyphens throughout this document (e.g. `pe-wgs-150bp-2x150`). In the underlying TSVs and pipeline configs the same datasets are named with underscores (e.g. `pe_wgs_150bp_2x150`). The hyphenated form is purely a display choice; Markdown italicizes stretches of `_text_text_`, which is noisy in tables.

## Datasets

Simulated with holodeck from GRCh38. Tier 1 scores accuracy on every tool; tier 2 times the fastest tools on larger data. Each ID links to its results.

| ID | Name | Tier | Layout | Insert (bp) | Error rate | Adapters | Size |
|---:|------|-----:|--------|------------:|-----------:|----------|-----:|
| [1](#dataset-1---pe-wgs-150bp-2x150) | pe-wgs-150bp-2x150 | 1 | PE 2×150 | 150 ± 30 | 0.1–1% | TruSeq | 1.07 M pairs |
| [2](#dataset-2---pe-wgs-250bp-2x150) | pe-wgs-250bp-2x150 | 1 | PE 2×150 | 250 ± 40 | 0.1–1% | TruSeq | 1.07 M pairs |
| [3](#dataset-3---pe-wgs-350bp-2x150) | pe-wgs-350bp-2x150 | 1 | PE 2×150 | 350 ± 60 | 0.1–1% | TruSeq | 1.07 M pairs |
| [4](#dataset-4---pe-wgs-450bp-2x150) | pe-wgs-450bp-2x150 | 1 | PE 2×150 | 450 ± 80 | 0.1–1% | TruSeq | 1.07 M pairs |
| [5](#dataset-5---pe-wgs-450bp-2x250) | pe-wgs-450bp-2x250 | 1 | PE 2×250 | 450 ± 80 | 0.1–1% | TruSeq | 0.64 M pairs |
| [6](#dataset-6---pe-wgs-higherr-250bp-2x150) | pe-wgs-higherr-250bp-2x150 | 1 | PE 2×150 | 250 ± 60 | 1–5% | TruSeq | 1.07 M pairs |
| [7](#dataset-7---pe-cfdna-170bp-2x150) | pe-cfdna-170bp-2x150 | 1 | PE 2×150 | 170 ± 30 | 0.1–1% | TruSeq | 1.07 M pairs |
| [8](#dataset-8---pe-exome-140bp-2x76-nextera) | pe-exome-140bp-2x76-nextera | 1 | PE 2×76 | 140 ± 25 | 0.1–1% | Nextera | 2.12 M pairs |
| [9](#dataset-9---se-wgs-300bp-1x150) | se-wgs-300bp-1x150 | 1 | SE 1×150 | 300 ± 80 | 0.1–1% | TruSeq | 2.14 M reads |
| [10](#dataset-10---se-short-120bp-1x150) | se-short-120bp-1x150 | 1 | SE 1×150 | 120 ± 30 | 0.1–1% | TruSeq | 2.14 M reads |
| [11](#dataset-11---se-mirna-30bp-1x76) | se-mirna-30bp-1x76 | 1 | SE 1×76 | 30 ± 2 | 0.1–1% | small-RNA | 4.23 M reads |
| [T1](#dataset-t1---pe-wgs-350bp-2x150-50m) | pe-wgs-350bp-2x150-50m | 2 | PE 2×150 | 350 ± 60 | 0.1–1% | TruSeq | 50.13 M pairs |
| [T2](#dataset-t2---pe-cfdna-170bp-2x150-50m) | pe-cfdna-170bp-2x150-50m | 2 | PE 2×150 | 170 ± 30 | 0.1–1% | TruSeq | 50.13 M pairs |

## Tier 1: accuracy and screen

Every tool on 11 datasets at 0.1× human WGS each, one run per cell. The `adapter_only` config (adapter trimming and a 30 bp length filter only) feeds the accuracy tables; accuracy counts sum both mates of paired-end data, and RMSE, MAE and exact-match rate pool them. The `wgs` config (adapter, poly-G, sliding-window quality and length trimming, N filter) feeds the screen timings. These runs last about a second for the fastest tools, so the timings only decide which tools go forward to tier 2: those no slower than fastp on the paired-end datasets.

### Dataset 1 - pe-wgs-150bp-2x150

PE 2×150, insert 150 ± 30 bp, error 0.1–1%, TruSeq adapters, 1.07 M pairs.

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.0391 | 0.0005 | 0.99968 | 36 | 580 | 0 | 0 |
| adapterremoval | 0.0505 | 0.0003 | 0.99992 | 104 | 0 | 0 | 0 |
| fastp | 0.1692 | 0.0047 | 0.99905 | 1,377 | 273 | 0 | 304 |
| bbduk | 1.1257 | 0.1940 | 0.96961 | 62,686 | 336 | 2,072 | 0 |
| cutadapt | 1.1713 | 0.2525 | 0.92804 | 156 | 154,043 | 35 | 13 |
| trim-galore-rs | 1.1713 | 0.2525 | 0.92804 | 156 | 154,043 | 35 | 13 |
| trimmomatic | 15.7109 | 4.1732 | 0.92604 | 862 | 157,646 | 106 | 0 |

#### Screen timings (wgs config, 8 threads, 1 run)

| Tool | Wall (s) | ×fastest | Pairs/s (M) | User CPU (s) | CPU % | Max RSS (MB) |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.64 | 1.00× | 1.68 | 4.3 | 695 | 115 |
| adapterremoval | 0.84 | 1.31× | 1.28 | 5.0 | 634 | 47 |
| cutadapt | 2.11 | 3.30× | 0.51 | 10.6 | 567 | 56 |
| trim-galore-rs | 2.51 | 3.92× | 0.43 | 18.3 | 736 | 133 |
| fastp | 2.91 | 4.55× | 0.37 | 18.9 | 679 | 1,198 |
| bbduk | 8.32 | 13.00× | 0.13 | 62.4 | 772 | 2,105 |
| trimmomatic | 15.90 | 24.84× | 0.07 | 100.4 | 637 | 3,984 |

### Dataset 2 - pe-wgs-250bp-2x150

PE 2×150, insert 250 ± 40 bp, error 0.1–1%, TruSeq adapters, 1.07 M pairs.

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.0133 | 0.0001 | 0.99996 | 70 | 18 | 0 | 0 |
| adapterremoval | 0.1125 | 0.0014 | 0.99977 | 502 | 0 | 0 | 0 |
| cutadapt | 0.1927 | 0.0072 | 0.99795 | 289 | 4,110 | 1 | 0 |
| trim-galore-rs | 0.1927 | 0.0072 | 0.99795 | 289 | 4,110 | 1 | 0 |
| trimmomatic | 0.7491 | 0.0162 | 0.99911 | 1,576 | 322 | 2 | 0 |
| fastp | 1.4259 | 0.0323 | 0.99823 | 3,787 | 6 | 0 | 6 |
| bbduk | 1.5557 | 0.3655 | 0.94347 | 121,164 | 8 | 68 | 0 |

#### Screen timings (wgs config, 8 threads, 1 run)

| Tool | Wall (s) | ×fastest | Pairs/s (M) | User CPU (s) | CPU % | Max RSS (MB) |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.67 | 1.00× | 1.60 | 4.5 | 694 | 127 |
| adapterremoval | 0.97 | 1.45× | 1.11 | 5.7 | 606 | 51 |
| cutadapt | 1.87 | 2.79× | 0.57 | 7.7 | 481 | 50 |
| trim-galore-rs | 1.99 | 2.97× | 0.54 | 14.2 | 725 | 112 |
| fastp | 3.75 | 5.60× | 0.29 | 25.8 | 706 | 1,195 |
| bbduk | 8.73 | 13.03× | 0.12 | 65.6 | 775 | 2,460 |
| trimmomatic | 14.27 | 21.30× | 0.08 | 103.5 | 733 | 4,585 |

### Dataset 3 - pe-wgs-350bp-2x150

PE 2×150, insert 350 ± 60 bp, error 0.1–1%, TruSeq adapters, 1.07 M pairs.

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.0084 | 0.0000 | 0.99998 | 42 | 2 | 0 | 0 |
| cutadapt | 0.0988 | 0.0014 | 0.99975 | 283 | 247 | 0 | 0 |
| trim-galore-rs | 0.0988 | 0.0014 | 0.99975 | 283 | 247 | 0 | 0 |
| adapterremoval | 0.1059 | 0.0013 | 0.99978 | 476 | 0 | 0 | 0 |
| trimmomatic | 0.4771 | 0.0092 | 0.99926 | 1,510 | 72 | 0 | 0 |
| bbduk | 1.7161 | 0.3739 | 0.94296 | 122,320 | 0 | 0 | 0 |
| fastp | 4.2507 | 0.1975 | 0.99635 | 7,817 | 0 | 0 | 0 |

#### Screen timings (wgs config, 8 threads, 1 run)

| Tool | Wall (s) | ×fastest | Pairs/s (M) | User CPU (s) | CPU % | Max RSS (MB) |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.67 | 1.00× | 1.60 | 4.7 | 701 | 121 |
| adapterremoval | 1.03 | 1.54× | 1.04 | 6.0 | 616 | 51 |
| cutadapt | 1.80 | 2.69× | 0.60 | 7.5 | 486 | 50 |
| trim-galore-rs | 1.86 | 2.78× | 0.58 | 13.4 | 731 | 114 |
| fastp | 3.96 | 5.91× | 0.27 | 26.9 | 700 | 1,196 |
| bbduk | 6.34 | 9.46× | 0.17 | 46.4 | 765 | 2,492 |
| trimmomatic | 18.10 | 27.01× | 0.06 | 98.7 | 551 | 4,215 |

### Dataset 4 - pe-wgs-450bp-2x150

PE 2×150, insert 450 ± 80 bp, error 0.1–1%, TruSeq adapters, 1.07 M pairs.

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.0109 | 0.0001 | 0.99997 | 58 | 0 | 0 | 0 |
| cutadapt | 0.0891 | 0.0010 | 0.99985 | 284 | 39 | 0 | 0 |
| trim-galore-rs | 0.0891 | 0.0010 | 0.99985 | 284 | 39 | 0 | 0 |
| adapterremoval | 0.1130 | 0.0015 | 0.99972 | 592 | 0 | 0 | 0 |
| trimmomatic | 0.3978 | 0.0094 | 0.99907 | 1,974 | 18 | 0 | 0 |
| bbduk | 1.7176 | 0.3772 | 0.94249 | 123,324 | 0 | 0 | 0 |
| fastp | 4.0589 | 0.1982 | 0.99577 | 9,063 | 0 | 0 | 0 |

#### Screen timings (wgs config, 8 threads, 1 run)

| Tool | Wall (s) | ×fastest | Pairs/s (M) | User CPU (s) | CPU % | Max RSS (MB) |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.66 | 1.00× | 1.62 | 4.5 | 701 | 119 |
| adapterremoval | 1.03 | 1.56× | 1.04 | 6.0 | 616 | 51 |
| trim-galore-rs | 1.80 | 2.73× | 0.60 | 12.9 | 729 | 122 |
| cutadapt | 1.87 | 2.83× | 0.57 | 7.7 | 482 | 48 |
| fastp | 3.91 | 5.92× | 0.27 | 26.5 | 703 | 1,195 |
| bbduk | 5.74 | 8.70× | 0.19 | 41.6 | 757 | 2,318 |
| trimmomatic | 15.48 | 23.45× | 0.07 | 98.8 | 645 | 4,712 |

### Dataset 5 - pe-wgs-450bp-2x250

PE 2×250, insert 450 ± 80 bp, error 0.1–1%, TruSeq adapters, 0.64 M pairs.

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.0096 | 0.0000 | 0.99997 | 32 | 4 | 0 | 0 |
| adapterremoval | 0.1301 | 0.0014 | 0.99981 | 238 | 0 | 0 | 0 |
| cutadapt | 0.1650 | 0.0045 | 0.99879 | 167 | 1,384 | 0 | 0 |
| trim-galore-rs | 0.1650 | 0.0045 | 0.99879 | 167 | 1,384 | 0 | 0 |
| trimmomatic | 1.3774 | 0.0305 | 0.99908 | 1,136 | 50 | 0 | 0 |
| bbduk | 1.7271 | 0.3687 | 0.94340 | 72,778 | 4 | 32 | 0 |
| fastp | 4.4790 | 0.1684 | 0.99699 | 3,866 | 2 | 0 | 4 |

#### Screen timings (wgs config, 8 threads, 1 run)

| Tool | Wall (s) | ×fastest | Pairs/s (M) | User CPU (s) | CPU % | Max RSS (MB) |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.59 | 1.00× | 1.09 | 4.1 | 697 | 127 |
| adapterremoval | 1.03 | 1.75× | 0.62 | 5.8 | 586 | 67 |
| trim-galore-rs | 1.75 | 2.97× | 0.37 | 12.4 | 715 | 155 |
| cutadapt | 1.76 | 2.98× | 0.37 | 6.2 | 420 | 48 |
| fastp | 3.64 | 6.17× | 0.18 | 24.5 | 693 | 1,224 |
| bbduk | 9.45 | 16.02× | 0.07 | 71.3 | 776 | 2,398 |
| trimmomatic | 15.10 | 25.59× | 0.04 | 97.9 | 655 | 4,238 |

### Dataset 6 - pe-wgs-higherr-250bp-2x150

PE 2×150, insert 250 ± 60 bp, error 1–5%, TruSeq adapters, 1.07 M pairs.

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| adapterremoval | 0.0946 | 0.0010 | 0.99972 | 356 | 0 | 0 | 0 |
| **chelae** | 0.4167 | 0.0079 | 0.99896 | 36 | 1,950 | 0 | 0 |
| cutadapt | 1.3142 | 0.1030 | 0.98696 | 344 | 27,342 | 8 | 30 |
| trim-galore-rs | 1.3142 | 0.1030 | 0.98696 | 344 | 27,342 | 8 | 30 |
| bbduk | 1.5822 | 0.3689 | 0.94211 | 119,356 | 4,092 | 360 | 96 |
| fastp | 1.9491 | 0.0514 | 0.99767 | 4,180 | 565 | 2 | 16 |
| trimmomatic | 4.0157 | 0.2009 | 0.99618 | 4,418 | 3,680 | 80 | 4 |

#### Screen timings (wgs config, 8 threads, 1 run)

| Tool | Wall (s) | ×fastest | Pairs/s (M) | User CPU (s) | CPU % | Max RSS (MB) |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.60 | 1.00× | 1.79 | 1.7 | 291 | 56 |
| adapterremoval | 0.70 | 1.17× | 1.53 | 3.7 | 531 | 35 |
| trim-galore-rs | 1.32 | 2.20× | 0.81 | 3.5 | 275 | 51 |
| cutadapt | 1.76 | 2.93× | 0.61 | 4.0 | 260 | 34 |
| fastp | 1.80 | 3.00× | 0.60 | 10.1 | 598 | 1,176 |
| trimmomatic | 2.41 | 4.02× | 0.44 | 14.7 | 634 | 2,363 |
| bbduk | 7.41 | 12.35× | 0.14 | 56.3 | 776 | 2,406 |

### Dataset 7 - pe-cfdna-170bp-2x150

PE 2×150, insert 170 ± 30 bp, error 0.1–1%, TruSeq adapters, 1.07 M pairs.

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.0408 | 0.0005 | 0.99975 | 60 | 464 | 0 | 0 |
| adapterremoval | 0.0771 | 0.0007 | 0.99989 | 236 | 0 | 0 | 0 |
| fastp | 0.2201 | 0.0063 | 0.99887 | 2,048 | 226 | 0 | 150 |
| cutadapt | 0.9506 | 0.1792 | 0.94693 | 204 | 113,587 | 14 | 11 |
| trim-galore-rs | 0.9506 | 0.1792 | 0.94693 | 204 | 113,587 | 14 | 11 |
| bbduk | 1.3636 | 0.2846 | 0.95565 | 93,184 | 310 | 1,612 | 0 |
| trimmomatic | 7.1962 | 0.9321 | 0.98225 | 1,284 | 36,726 | 62 | 0 |

#### Screen timings (wgs config, 8 threads, 1 run)

| Tool | Wall (s) | ×fastest | Pairs/s (M) | User CPU (s) | CPU % | Max RSS (MB) |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.63 | 1.00× | 1.70 | 4.2 | 697 | 117 |
| adapterremoval | 0.90 | 1.43× | 1.19 | 5.2 | 616 | 51 |
| cutadapt | 2.05 | 3.25× | 0.52 | 9.3 | 526 | 49 |
| trim-galore-rs | 2.14 | 3.40× | 0.50 | 15.7 | 741 | 133 |
| fastp | 3.19 | 5.06× | 0.34 | 21.3 | 691 | 1,197 |
| bbduk | 8.69 | 13.79× | 0.12 | 65.3 | 772 | 2,017 |
| trimmomatic | 15.94 | 25.30× | 0.07 | 101.6 | 642 | 3,550 |

### Dataset 8 - pe-exome-140bp-2x76-nextera

PE 2×76, insert 140 ± 25 bp, error 0.1–1%, Nextera adapters, 2.12 M pairs.

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.0128 | 0.0001 | 0.99996 | 130 | 34 | 0 | 0 |
| adapterremoval | 0.0744 | 0.0010 | 0.99974 | 1,070 | 0 | 0 | 0 |
| cutadapt | 0.2805 | 0.0137 | 0.99671 | 4,032 | 9,858 | 7 | 0 |
| trim-galore-rs | 0.2805 | 0.0137 | 0.99671 | 4,032 | 9,858 | 7 | 0 |
| trimmomatic | 0.6469 | 0.0190 | 0.99867 | 3,026 | 2,606 | 12 | 0 |
| fastp | 1.0139 | 0.0489 | 0.99559 | 18,623 | 10 | 0 | 6 |
| bbduk | 1.6165 | 0.3937 | 0.93942 | 255,558 | 36 | 842 | 0 |

#### Screen timings (wgs config, 8 threads, 1 run)

| Tool | Wall (s) | ×fastest | Pairs/s (M) | User CPU (s) | CPU % | Max RSS (MB) |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.82 | 1.00× | 2.58 | 5.7 | 715 | 91 |
| adapterremoval | 0.98 | 1.20× | 2.16 | 6.3 | 676 | 38 |
| trim-galore-rs | 2.27 | 2.77× | 0.93 | 16.7 | 738 | 92 |
| cutadapt | 2.50 | 3.05× | 0.85 | 11.7 | 540 | 53 |
| fastp | 3.48 | 4.24× | 0.61 | 22.7 | 685 | 1,175 |
| bbduk | 6.26 | 7.63× | 0.34 | 44.8 | 748 | 2,488 |
| trimmomatic | 19.51 | 23.79× | 0.11 | 106.2 | 550 | 4,430 |

### Dataset 9 - se-wgs-300bp-1x150

SE 1×150, insert 300 ± 80 bp, error 0.1–1%, TruSeq adapters, 2.14 M reads.

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.2652 | 0.0149 | 0.99530 | 641 | 8,667 | 2 | 0 |
| cutadapt | 0.2701 | 0.0147 | 0.99536 | 259 | 8,913 | 1 | 0 |
| trim-galore-rs | 0.2701 | 0.0147 | 0.99536 | 259 | 8,913 | 1 | 0 |
| fastp | 0.2805 | 0.0161 | 0.99514 | 2,587 | 7,061 | 2 | 0 |
| adapterremoval | 0.4122 | 0.0308 | 0.99275 | 6,291 | 8,465 | 16 | 0 |
| bbduk | 0.9019 | 0.1236 | 0.97849 | 36,634 | 8,597 | 73 | 59 |
| trimmomatic | 1.1722 | 0.0603 | 0.99270 | 4,513 | 10,242 | 27 | 0 |

#### Screen timings (wgs config, 8 threads, 1 run)

| Tool | Wall (s) | ×fastest | Reads/s (M) | User CPU (s) | CPU % | Max RSS (MB) |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 1.13 | 1.00× | 1.90 | 4.4 | 400 | 92 |
| adapterremoval | 1.23 | 1.09× | 1.74 | 4.8 | 428 | 29 |
| cutadapt | 2.10 | 1.86× | 1.02 | 7.5 | 414 | 42 |
| bbduk | 2.14 | 1.89× | 1.00 | 11.2 | 598 | 3,009 |
| trim-galore-rs | 2.38 | 2.11× | 0.90 | 12.5 | 530 | 51 |
| fastp | 3.12 | 2.76× | 0.69 | 20.2 | 681 | 1,109 |
| trimmomatic | 20.42 | 18.07× | 0.11 | 94.6 | 466 | 2,312 |

### Dataset 10 - se-short-120bp-1x150

SE 1×150, insert 120 ± 30 bp, error 0.1–1%, TruSeq adapters, 2.14 M reads.

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| fastp | 0.6499 | 0.1053 | 0.96176 | 445 | 78,757 | 52 | 7 |
| adapterremoval | 0.7421 | 0.1413 | 0.95383 | 1,089 | 95,008 | 180 | 0 |
| **chelae** | 0.8007 | 0.1491 | 0.95301 | 128 | 97,868 | 45 | 0 |
| cutadapt | 0.9635 | 0.1702 | 0.95112 | 54 | 102,000 | 19 | 13 |
| trim-galore-rs | 0.9635 | 0.1702 | 0.95112 | 54 | 102,000 | 19 | 13 |
| trimmomatic | 1.2042 | 0.2240 | 0.94228 | 763 | 119,845 | 357 | 0 |
| bbduk | 1.6140 | 0.2273 | 0.94817 | 6,058 | 99,849 | 658 | 1,858 |

#### Screen timings (wgs config, 8 threads, 1 run)

| Tool | Wall (s) | ×fastest | Reads/s (M) | User CPU (s) | CPU % | Max RSS (MB) |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 1.13 | 1.00× | 1.90 | 3.8 | 350 | 89 |
| adapterremoval | 1.18 | 1.04× | 1.82 | 4.4 | 410 | 27 |
| bbduk | 2.01 | 1.78× | 1.07 | 10.4 | 591 | 2,819 |
| cutadapt | 2.38 | 2.11× | 0.90 | 13.9 | 640 | 39 |
| fastp | 2.48 | 2.19× | 0.86 | 14.3 | 610 | 1,111 |
| trim-galore-rs | 3.29 | 2.91× | 0.65 | 25.2 | 768 | 65 |
| trimmomatic | 20.65 | 18.27× | 0.10 | 85.9 | 420 | 3,828 |

### Dataset 11 - se-mirna-30bp-1x76

SE 1×76, insert 30 ± 2 bp, error 0.1–1%, small-RNA adapters, 4.23 M reads.

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| adapterremoval | 0.0000 | 0.0000 | 0.59912 | 0 | 0 | 0 | 0 |
| trimmomatic | 0.0404 | 0.0000 | 0.59912 | 0 | 2 | 0 | 0 |
| **chelae** | 0.0655 | 0.0001 | 0.59912 | 0 | 5 | 0 | 0 |
| cutadapt | 0.0655 | 0.0001 | 0.59912 | 0 | 5 | 0 | 0 |
| fastp | 0.0655 | 0.0001 | 0.59912 | 0 | 5 | 0 | 0 |
| trim-galore-rs | 0.0655 | 0.0001 | 0.59912 | 0 | 5 | 0 | 0 |
| bbduk | 0.8428 | 0.0155 | 0.59900 | 0 | 858 | 0 | 0 |

#### Screen timings (wgs config, 8 threads, 1 run)

| Tool | Wall (s) | ×fastest | Reads/s (M) | User CPU (s) | CPU % | Max RSS (MB) |
|---|---:|---:|---:|---:|---:|---:|
| adapterremoval | 1.17 | 1.00× | 3.62 | 3.4 | 314 | 19 |
| **chelae** | 1.31 | 1.12× | 3.23 | 3.0 | 244 | 47 |
| bbduk | 1.69 | 1.44× | 2.50 | 8.2 | 562 | 2,851 |
| fastp | 2.33 | 1.99× | 1.82 | 8.8 | 422 | 1,100 |
| cutadapt | 2.48 | 2.12× | 1.71 | 16.1 | 694 | 31 |
| trim-galore-rs | 2.95 | 2.52× | 1.44 | 21.1 | 716 | 51 |
| trimmomatic | 6.95 | 5.94× | 0.61 | 41.3 | 605 | 2,509 |

## Tier 2: runtime

The tools that passed tier 1's screen, on two paired-end datasets of ~50 M pairs each, `wgs` config, three runs per cell. Wall time is the median, with the range in brackets; the other columns are medians too.

### Dataset T1 - pe-wgs-350bp-2x150-50m

PE 2×150, insert 350 ± 60 bp, error 0.1–1%, TruSeq adapters, 50.13 M pairs.

#### Runtime (wgs config, 8 threads, median of 3)

| Tool | Wall (s) | ×fastest | Pairs/s (M) | User CPU (s) | CPU % | Max RSS (MB) | Cycles (G) |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 29.68 (29.3–29.9) | 1.00× | 1.69 | 206.5 | 723 | 117 | 831 |
| adapterremoval | 37.37 (37.1–38.4) | 1.26× | 1.34 | 256.7 | 717 | 55 | 1,056 |
| cutadapt | 77.56 (77.3–81.7) | 2.61× | 0.65 | 358.3 | 526 | 52 | 1,567 |
| trim-galore-rs | 79.98 (79.4–80.4) | 2.69× | 0.63 | 617.8 | 783 | 131 | 2,479 |
| fastp | 169.57 (168.6–169.9) | 5.71× | 0.30 | 1251.0 | 758 | 1,274 | 5,001 |

### Dataset T2 - pe-cfdna-170bp-2x150-50m

PE 2×150, insert 170 ± 30 bp, error 0.1–1%, TruSeq adapters, 50.13 M pairs.

#### Runtime (wgs config, 8 threads, median of 3)

| Tool | Wall (s) | ×fastest | Pairs/s (M) | User CPU (s) | CPU % | Max RSS (MB) | Cycles (G) |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 28.97 (28.9–29.1) | 1.00× | 1.73 | 201.8 | 719 | 113 | 814 |
| adapterremoval | 36.14 (36.0–36.1) | 1.25× | 1.39 | 237.2 | 702 | 54 | 974 |
| cutadapt | 81.94 (81.4–84.1) | 2.83× | 0.61 | 418.0 | 572 | 76 | 1,782 |
| trim-galore-rs | 93.94 (93.2–94.1) | 3.24× | 0.53 | 730.9 | 787 | 140 | 2,922 |
| fastp | 132.54 (132.4–134.4) | 4.58× | 0.38 | 973.8 | 757 | 1,223 | 3,901 |
