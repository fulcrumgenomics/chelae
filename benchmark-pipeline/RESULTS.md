# Benchmark results - May 8th 2026
Companion to the [Performance](../README.md#performance) section in the project README.  This document provides more in depth infromation on all tools from the same set of benchmarking runs referenced in the project README.

## Run metadata

- **Date**: 2026-05-08
- **Host**: EC2 c8id.2xlarge - Intel Xeon 6975P-C (Granite Rapids), 8 logical CPUs, 16 GB RAM
- **OS**: Amazon Linux 2023 (kernel 6.1)
- **chelae version**: 0.1.0
- **Replicates**: 1 per cell in the cross-tool suite, up to 3 per cell in the throughput suite

For tool versions see [the project README section "Tool versions"](../README.md#tool-versions).

> **Naming note**: dataset names are written with hyphens throughout this
> document (e.g. `pe-wgs-150bp-2x150`). In the underlying TSVs and pipeline
> configs the same datasets are named with underscores (e.g. `pe_wgs_150bp_2x150`).
> The hyphenated form is purely a display choice — Markdown italicizes
> stretches of `_text_text_`, which is noisy in tables.

## Table of contents

- [Cross-tool comparison suite](#cross-tool-comparison-suite)
  - [Dataset 1 - pe-wgs-150bp-2x150](#dataset-1---pe-wgs-150bp-2x150) - PE 2x150, insert 150 +/- 30, TruSeq adapters
  - [Dataset 2 - pe-wgs-250bp-2x150](#dataset-2---pe-wgs-250bp-2x150) - PE 2x150, insert 250 +/- 40, TruSeq adapters
  - [Dataset 3 - pe-wgs-350bp-2x150](#dataset-3---pe-wgs-350bp-2x150) - PE 2x150, insert 350 +/- 60, TruSeq adapters
  - [Dataset 4 - pe-wgs-450bp-2x150](#dataset-4---pe-wgs-450bp-2x150) - PE 2x150, insert 450 +/- 80, TruSeq adapters
  - [Dataset 5 - pe-wgs-450bp-2x250](#dataset-5---pe-wgs-450bp-2x250) - PE 2x250, insert 450 +/- 80, TruSeq adapters
  - [Dataset 6 - pe-wgs-higherr-250bp-2x150](#dataset-6---pe-wgs-higherr-250bp-2x150) - PE 2x150, insert 250 +/- 60, TruSeq adapters
  - [Dataset 7 - pe-cfdna-170bp-2x150](#dataset-7---pe-cfdna-170bp-2x150) - PE 2x150, insert 170 +/- 30, TruSeq adapters
  - [Dataset 8 - pe-exome-140bp-2x76-nextera](#dataset-8---pe-exome-140bp-2x76-nextera) - PE 2x76, insert 140 +/- 25, Nextera adapters
  - [Dataset 9 - se-wgs-300bp-1x150](#dataset-9---se-wgs-300bp-1x150) - SE 1x150, insert 300 +/- 80, TruSeq adapters
  - [Dataset 10 - se-short-120bp-1x150](#dataset-10---se-short-120bp-1x150) - SE 1x150, insert 120 +/- 30, TruSeq adapters
  - [Dataset 11 - se-mirna-30bp-1x76](#dataset-11---se-mirna-30bp-1x76) - SE 1x76, insert 30 +/- 2, small-RNA adapters
- [Throughput / scalability suite](#throughput--scalability-suite)
  - [Dataset T1 - pe-wgs-150](#dataset-t1---pe-wgs-150) - PE 2x150, insert 350 +/- 80, TruSeq adapters
  - [Dataset T2 - pe-cfdna-150](#dataset-t2---pe-cfdna-150) - PE 2x150, insert 170 +/- 40, TruSeq adapters

## Cross-tool comparison suite

The cross-tool comparison suite is comprised of 11 datasets designed to exercise adapter trimming across a broad range of scenarios.  Approximately 2X human WGS data was simulated for each dataset - large enough for prelimenary performance comparisons, but small enough that slower tools could still complete in a reasoanble time.  The primary goal of this suite is to assess **accuracy**, with the secondary goal of classifying tools into high vs. lower performance tiers.

Only one (1) replicate per cell was run. Every tool runs on every dataset. A `wgs` config (adapter + quality + length filter) feeds the runtime tables; an `adapter_only` config feeds the accuracy tables.

### Dataset 1 - pe-wgs-150bp-2x150

PE 2x150, insert 150 +/- 30 bp, error 0.1-1%, TruSeq adapters, 21.4 M pairs.

#### Runtime (wgs config, 8 threads)

| Tool | Wall (s) | xfastest | Reads/s (M) | Max RSS (MB) | User CPU (s) | CPU % |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 37.42 | 1.00x | 1.14 | 70 | 289.5 | 781 |
| cutadapt | 78.15 | 2.09x | 0.55 | 76 | 529.1 | 711 |
| trim-galore-rs | 78.23 | 2.09x | 0.55 | 150 | 605.2 | 778 |
| fastp | 95.71 | 2.56x | 0.45 | 1,233 | 704.2 | 742 |
| fastp-nfcore | 106.67 | 2.85x | 0.40 | 1,380 | 712.1 | 674 |
| adapterremoval | 168.80 | 4.51x | 0.25 | 93 | 565.6 | 337 |
| bbduk | 221.77 | 5.93x | 0.19 | 1,548 | 1752.7 | 795 |
| trimmomatic | 465.52 | 12.44x | 0.09 | 4,954 | 2478.8 | 535 |
| trim-galore | 794.93 | 21.24x | 0.05 | 38 | 4606.5 | 587 |

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.0448 | 0.0006 | 0.99969 | 878 | 11,202 | 4 | 0 |
| adapterremoval | 0.0791 | 0.0005 | 0.99990 | 2,974 | 0 | 0 | 0 |
| fastp | 0.1685 | 0.0047 | 0.99905 | 27,916 | 5,442 | 3 | 6,082 |
| fastp-nfcore | 0.1685 | 0.0047 | 0.99905 | 27,916 | 5,442 | 3 | 6,082 |
| bbduk | 1.1247 | 0.1939 | 0.96962 | 1,247,952 | 7,040 | 43,688 | 0 |
| cutadapt | 1.1746 | 0.2533 | 0.92810 | 2,896 | 3,071,517 | 585 | 429 |
| trim-galore | 1.1746 | 0.2533 | 0.92810 | 2,896 | 3,071,517 | 585 | 429 |
| trim-galore-rs | 1.1746 | 0.2533 | 0.92810 | 2,896 | 3,071,517 | 585 | 429 |
| trimmomatic | 15.7372 | 4.1834 | 0.92591 | 17,462 | 3,150,896 | 1,826 | 0 |

### Dataset 2 - pe-wgs-250bp-2x150

PE 2x150, insert 250 +/- 40 bp, error 0.1-1%, TruSeq adapters, 21.4 M pairs.

#### Runtime (wgs config, 8 threads)

| Tool | Wall (s) | xfastest | Reads/s (M) | Max RSS (MB) | User CPU (s) | CPU % |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 41.86 | 1.00x | 1.02 | 75 | 321.8 | 776 |
| trim-galore-rs | 48.39 | 1.16x | 0.88 | 91 | 358.2 | 750 |
| cutadapt | 67.96 | 1.62x | 0.63 | 77 | 435.1 | 677 |
| fastp-nfcore | 125.43 | 3.00x | 0.34 | 1,347 | 913.8 | 735 |
| fastp | 125.58 | 3.00x | 0.34 | 1,210 | 921.1 | 739 |
| adapterremoval | 184.83 | 4.42x | 0.23 | 84 | 647.8 | 352 |
| bbduk | 234.01 | 5.59x | 0.18 | 1,600 | 1848.3 | 795 |
| trimmomatic | 472.31 | 11.28x | 0.09 | 4,634 | 2498.8 | 532 |
| trim-galore | 808.38 | 19.31x | 0.05 | 40 | 4578.5 | 575 |

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.0127 | 0.0001 | 0.99996 | 1,468 | 302 | 0 | 0 |
| cutadapt | 0.1891 | 0.0071 | 0.99795 | 5,739 | 81,780 | 23 | 2 |
| trim-galore | 0.1891 | 0.0071 | 0.99795 | 5,739 | 81,780 | 23 | 2 |
| trim-galore-rs | 0.1891 | 0.0071 | 0.99795 | 5,739 | 81,780 | 23 | 2 |
| adapterremoval | 0.5729 | 0.0143 | 0.99879 | 51,936 | 0 | 0 | 0 |
| trimmomatic | 0.7903 | 0.0176 | 0.99905 | 33,346 | 7,186 | 64 | 0 |
| fastp | 1.4151 | 0.0321 | 0.99826 | 74,269 | 164 | 0 | 72 |
| fastp-nfcore | 1.4151 | 0.0321 | 0.99826 | 74,269 | 164 | 0 | 72 |
| bbduk | 1.5653 | 0.3688 | 0.94293 | 2,440,488 | 206 | 1,296 | 0 |

### Dataset 3 - pe-wgs-350bp-2x150

PE 2x150, insert 350 +/- 60 bp, error 0.1-1%, TruSeq adapters, 21.4 M pairs.

#### Runtime (wgs config, 8 threads)

| Tool | Wall (s) | xfastest | Reads/s (M) | Max RSS (MB) | User CPU (s) | CPU % |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 42.24 | 1.00x | 1.01 | 73 | 319.4 | 764 |
| trim-galore-rs | 48.62 | 1.15x | 0.88 | 85 | 354.9 | 741 |
| cutadapt | 65.75 | 1.56x | 0.65 | 83 | 445.0 | 715 |
| fastp | 129.72 | 3.07x | 0.33 | 1,239 | 968.3 | 751 |
| fastp-nfcore | 129.95 | 3.08x | 0.33 | 1,392 | 960.4 | 745 |
| bbduk | 168.85 | 4.00x | 0.25 | 972 | 1328.3 | 792 |
| adapterremoval | 185.28 | 4.39x | 0.23 | 88 | 672.8 | 365 |
| trimmomatic | 456.52 | 10.81x | 0.09 | 5,615 | 2504.4 | 552 |
| trim-galore | 807.21 | 19.11x | 0.05 | 39 | 4566.7 | 574 |

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.0093 | 0.0000 | 0.99998 | 1,030 | 12 | 0 | 0 |
| cutadapt | 0.0976 | 0.0014 | 0.99975 | 5,784 | 4,909 | 0 | 0 |
| trim-galore | 0.0976 | 0.0014 | 0.99975 | 5,784 | 4,909 | 0 | 0 |
| trim-galore-rs | 0.0976 | 0.0014 | 0.99975 | 5,784 | 4,909 | 0 | 0 |
| trimmomatic | 0.4583 | 0.0093 | 0.99925 | 31,006 | 1,052 | 2 | 0 |
| adapterremoval | 0.7625 | 0.0232 | 0.99813 | 79,862 | 0 | 0 | 0 |
| bbduk | 1.7995 | 0.3802 | 0.94243 | 2,463,164 | 12 | 82 | 0 |
| fastp | 4.5922 | 0.2293 | 0.99595 | 173,475 | 6 | 0 | 6 |
| fastp-nfcore | 4.5922 | 0.2293 | 0.99595 | 173,475 | 6 | 0 | 6 |

### Dataset 4 - pe-wgs-450bp-2x150

PE 2x150, insert 450 +/- 80 bp, error 0.1-1%, TruSeq adapters, 21.4 M pairs.

#### Runtime (wgs config, 8 threads)

| Tool | Wall (s) | xfastest | Reads/s (M) | Max RSS (MB) | User CPU (s) | CPU % |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 42.02 | 1.00x | 1.02 | 63 | 321.4 | 773 |
| trim-galore-rs | 49.25 | 1.17x | 0.87 | 92 | 352.4 | 726 |
| cutadapt | 65.72 | 1.56x | 0.65 | 82 | 440.8 | 711 |
| fastp | 129.65 | 3.09x | 0.33 | 1,240 | 970.9 | 754 |
| fastp-nfcore | 131.46 | 3.13x | 0.33 | 1,393 | 964.0 | 739 |
| bbduk | 159.02 | 3.78x | 0.27 | 979 | 1246.7 | 790 |
| adapterremoval | 184.81 | 4.40x | 0.23 | 90 | 674.1 | 366 |
| trimmomatic | 441.73 | 10.51x | 0.10 | 5,678 | 2503.1 | 570 |
| trim-galore | 803.59 | 19.12x | 0.05 | 40 | 4565.2 | 577 |

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.0103 | 0.0001 | 0.99997 | 1,244 | 6 | 0 | 0 |
| cutadapt | 0.0903 | 0.0011 | 0.99984 | 5,795 | 895 | 0 | 0 |
| trim-galore | 0.0903 | 0.0011 | 0.99984 | 5,795 | 895 | 0 | 0 |
| trim-galore-rs | 0.0903 | 0.0011 | 0.99984 | 5,795 | 895 | 0 | 0 |
| trimmomatic | 0.4285 | 0.0110 | 0.99893 | 45,228 | 392 | 2 | 0 |
| adapterremoval | 1.0189 | 0.0405 | 0.99722 | 119,144 | 0 | 0 | 0 |
| bbduk | 1.7864 | 0.3815 | 0.94225 | 2,470,992 | 4 | 22 | 0 |
| fastp | 4.3102 | 0.2245 | 0.99537 | 198,118 | 4 | 0 | 2 |
| fastp-nfcore | 4.3102 | 0.2245 | 0.99537 | 198,118 | 4 | 0 | 2 |

### Dataset 5 - pe-wgs-450bp-2x250

PE 2x250, insert 450 +/- 80 bp, error 0.1-1%, TruSeq adapters, 12.8 M pairs.

#### Runtime (wgs config, 8 threads)

| Tool | Wall (s) | xfastest | Reads/s (M) | Max RSS (MB) | User CPU (s) | CPU % |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 37.92 | 1.00x | 0.68 | 110 | 297.0 | 790 |
| trim-galore-rs | 44.30 | 1.17x | 0.58 | 121 | 325.2 | 743 |
| cutadapt | 57.52 | 1.52x | 0.45 | 77 | 373.9 | 694 |
| fastp-nfcore | 119.04 | 3.14x | 0.22 | 1,453 | 890.6 | 754 |
| fastp | 119.11 | 3.14x | 0.22 | 1,256 | 896.2 | 757 |
| adapterremoval | 173.18 | 4.57x | 0.15 | 127 | 654.3 | 379 |
| bbduk | 250.58 | 6.61x | 0.10 | 1,142 | 1975.3 | 792 |
| trimmomatic | 516.14 | 13.61x | 0.05 | 6,494 | 2362.2 | 460 |
| trim-galore | 798.46 | 21.06x | 0.03 | 40 | 4484.2 | 569 |

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.0109 | 0.0001 | 0.99997 | 678 | 100 | 0 | 0 |
| cutadapt | 0.1663 | 0.0047 | 0.99879 | 3,550 | 27,517 | 6 | 8 |
| trim-galore | 0.1663 | 0.0047 | 0.99879 | 3,550 | 27,517 | 6 | 8 |
| trim-galore-rs | 0.1663 | 0.0047 | 0.99879 | 3,550 | 27,517 | 6 | 8 |
| trimmomatic | 1.4231 | 0.0335 | 0.99896 | 25,454 | 1,254 | 42 | 0 |
| bbduk | 1.7161 | 0.3724 | 0.94289 | 1,465,860 | 42 | 432 | 0 |
| adapterremoval | 4.1232 | 0.1963 | 0.99645 | 91,166 | 0 | 0 | 0 |
| fastp | 4.7069 | 0.1897 | 0.99672 | 84,265 | 40 | 0 | 24 |
| fastp-nfcore | 4.7069 | 0.1897 | 0.99672 | 84,265 | 40 | 0 | 24 |

### Dataset 6 - pe-wgs-higherr-250bp-2x150

PE 2x150, insert 250 +/- 60 bp, error 1-5%, TruSeq adapters, 21.4 M pairs.

#### Runtime (wgs config, 8 threads)

| Tool | Wall (s) | xfastest | Reads/s (M) | Max RSS (MB) | User CPU (s) | CPU % |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 20.27 | 1.00x | 2.11 | 31 | 88.7 | 443 |
| trim-galore-rs | 32.88 | 1.62x | 1.30 | 21 | 87.1 | 270 |
| fastp-nfcore | 39.44 | 1.95x | 1.08 | 1,223 | 297.3 | 761 |
| fastp | 39.71 | 1.96x | 1.08 | 1,182 | 289.5 | 738 |
| cutadapt | 43.57 | 2.15x | 0.98 | 36 | 97.2 | 257 |
| trimmomatic | 45.86 | 2.26x | 0.93 | 2,641 | 188.1 | 421 |
| adapterremoval | 135.28 | 6.67x | 0.32 | 64 | 436.1 | 323 |
| trim-galore | 165.98 | 8.19x | 0.26 | 32 | 724.4 | 453 |
| bbduk | 176.54 | 8.71x | 0.24 | 1,142 | 1393.7 | 791 |

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.3932 | 0.0074 | 0.99897 | 820 | 38,268 | 0 | 0 |
| adapterremoval | 0.4819 | 0.0112 | 0.99882 | 45,208 | 0 | 0 | 0 |
| cutadapt | 1.3104 | 0.1019 | 0.98704 | 6,552 | 542,091 | 146 | 642 |
| trim-galore | 1.3104 | 0.1019 | 0.98704 | 6,552 | 542,091 | 146 | 642 |
| trim-galore-rs | 1.3104 | 0.1019 | 0.98704 | 6,552 | 542,091 | 146 | 642 |
| bbduk | 1.5933 | 0.3687 | 0.94221 | 2,377,962 | 81,568 | 6,648 | 1,434 |
| fastp | 1.9895 | 0.0524 | 0.99766 | 83,838 | 10,754 | 50 | 328 |
| fastp-nfcore | 1.9895 | 0.0524 | 0.99766 | 83,838 | 10,754 | 50 | 328 |
| trimmomatic | 4.0191 | 0.2007 | 0.99616 | 89,752 | 72,828 | 1,804 | 82 |

### Dataset 7 - pe-cfdna-170bp-2x150

PE 2x150, insert 170 +/- 30 bp, error 0.1-1%, TruSeq adapters, 21.4 M pairs.

#### Runtime (wgs config, 8 threads)

| Tool | Wall (s) | xfastest | Reads/s (M) | Max RSS (MB) | User CPU (s) | CPU % |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 39.88 | 1.00x | 1.07 | 79 | 305.9 | 774 |
| trim-galore-rs | 60.33 | 1.51x | 0.71 | 135 | 460.6 | 770 |
| cutadapt | 71.17 | 1.78x | 0.60 | 79 | 482.6 | 713 |
| fastp | 105.61 | 2.65x | 0.41 | 1,218 | 781.9 | 747 |
| fastp-nfcore | 112.86 | 2.83x | 0.38 | 1,362 | 786.3 | 703 |
| adapterremoval | 177.51 | 4.45x | 0.24 | 89 | 593.8 | 336 |
| bbduk | 237.09 | 5.95x | 0.18 | 1,251 | 1866.8 | 792 |
| trimmomatic | 478.04 | 11.99x | 0.09 | 4,567 | 2453.8 | 516 |
| trim-galore | 805.75 | 20.20x | 0.05 | 39 | 4591.6 | 578 |

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.0395 | 0.0005 | 0.99977 | 1,330 | 8,580 | 0 | 2 |
| adapterremoval | 0.1391 | 0.0013 | 0.99983 | 7,314 | 0 | 0 | 0 |
| fastp | 0.2066 | 0.0063 | 0.99886 | 41,525 | 4,226 | 4 | 2,892 |
| fastp-nfcore | 0.2066 | 0.0063 | 0.99886 | 41,525 | 4,226 | 4 | 2,892 |
| cutadapt | 0.9544 | 0.1805 | 0.94677 | 4,351 | 2,272,684 | 454 | 215 |
| trim-galore | 0.9544 | 0.1805 | 0.94677 | 4,351 | 2,272,684 | 454 | 215 |
| trim-galore-rs | 0.9544 | 0.1805 | 0.94677 | 4,351 | 2,272,684 | 454 | 215 |
| bbduk | 1.3630 | 0.2841 | 0.95576 | 1,854,554 | 5,264 | 33,344 | 0 |
| trimmomatic | 7.2354 | 0.9415 | 0.98206 | 27,362 | 739,162 | 1,294 | 0 |

### Dataset 8 - pe-exome-140bp-2x76-nextera

PE 2x76, insert 140 +/- 25 bp, error 0.1-1%, Nextera adapters, 42.2 M pairs.

#### Runtime (wgs config, 8 threads)

| Tool | Wall (s) | xfastest | Reads/s (M) | Max RSS (MB) | User CPU (s) | CPU % |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 50.46 | 1.00x | 1.67 | 76 | 375.5 | 754 |
| trim-galore-rs | 61.46 | 1.22x | 1.37 | 62 | 430.0 | 705 |
| cutadapt | 84.75 | 1.68x | 1.00 | 83 | 587.7 | 732 |
| fastp | 113.60 | 2.25x | 0.74 | 1,201 | 832.7 | 739 |
| fastp-nfcore | 127.23 | 2.52x | 0.66 | 1,321 | 835.0 | 661 |
| bbduk | 175.02 | 3.47x | 0.48 | 1,294 | 1367.9 | 788 |
| adapterremoval | 208.66 | 4.14x | 0.40 | 64 | 705.8 | 340 |
| trimmomatic | 478.36 | 9.48x | 0.18 | 3,112 | 2661.7 | 559 |
| trim-galore | 862.19 | 17.09x | 0.10 | 39 | 4952.4 | 585 |

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.0156 | 0.0001 | 0.99995 | 2,975 | 756 | 2 | 0 |
| cutadapt | 0.2836 | 0.0139 | 0.99672 | 81,842 | 194,846 | 268 | 0 |
| trim-galore | 0.2836 | 0.0139 | 0.99672 | 81,842 | 194,846 | 268 | 0 |
| trim-galore-rs | 0.2836 | 0.0139 | 0.99672 | 81,842 | 194,846 | 268 | 0 |
| adapterremoval | 0.3119 | 0.0079 | 0.99885 | 96,782 | 0 | 0 | 0 |
| trimmomatic | 0.6280 | 0.0182 | 0.99870 | 59,058 | 50,392 | 180 | 0 |
| fastp | 0.9937 | 0.0483 | 0.99553 | 376,236 | 209 | 5 | 74 |
| fastp-nfcore | 0.9937 | 0.0483 | 0.99553 | 376,236 | 209 | 5 | 74 |
| bbduk | 1.6227 | 0.3964 | 0.93901 | 5,131,460 | 752 | 18,654 | 0 |

### Dataset 9 - se-wgs-300bp-1x150

SE 1x150, insert 300 +/- 80 bp, error 0.1-1%, TruSeq adapters, 42.8 M reads.

#### Runtime (wgs config, 8 threads)

| Tool | Wall (s) | xfastest | Reads/s (M) | Max RSS (MB) | User CPU (s) | CPU % |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 42.77 | 1.00x | 1.00 | 39 | 328.9 | 774 |
| cutadapt | 68.85 | 1.61x | 0.62 | 83 | 409.2 | 628 |
| trim-galore-rs | 72.29 | 1.69x | 0.59 | 27 | 328.7 | 460 |
| bbduk | 77.40 | 1.81x | 0.55 | 1,322 | 586.2 | 771 |
| fastp | 108.71 | 2.54x | 0.39 | 1,119 | 802.7 | 744 |
| fastp-nfcore | 174.57 | 4.08x | 0.25 | 1,205 | 761.6 | 441 |
| adapterremoval | 320.93 | 7.50x | 0.13 | 59 | 512.1 | 160 |
| trim-galore | 492.23 | 11.51x | 0.09 | 40 | 2151.2 | 445 |
| trimmomatic | 505.84 | 11.83x | 0.08 | 3,014 | 2431.3 | 484 |

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| **chelae** | 0.2670 | 0.0153 | 0.99519 | 13,595 | 176,973 | 63 | 1 |
| cutadapt | 0.2717 | 0.0150 | 0.99526 | 5,605 | 182,076 | 36 | 10 |
| trim-galore | 0.2717 | 0.0150 | 0.99526 | 5,605 | 182,076 | 36 | 10 |
| trim-galore-rs | 0.2717 | 0.0150 | 0.99526 | 5,605 | 182,076 | 36 | 10 |
| fastp | 0.2843 | 0.0166 | 0.99499 | 53,203 | 145,798 | 70 | 9 |
| fastp-nfcore | 0.2843 | 0.0166 | 0.99499 | 53,203 | 145,798 | 70 | 9 |
| bbduk | 0.9010 | 0.1234 | 0.97845 | 729,046 | 175,568 | 1,239 | 969 |
| trimmomatic | 1.1426 | 0.0596 | 0.99260 | 87,579 | 211,529 | 523 | 0 |
| adapterremoval | 1.3832 | 0.1479 | 0.98330 | 525,726 | 168,008 | 4,823 | 0 |

### Dataset 10 - se-short-120bp-1x150

SE 1x150, insert 120 +/- 30 bp, error 0.1-1%, TruSeq adapters, 42.8 M reads.

#### Runtime (wgs config, 8 threads)

| Tool | Wall (s) | xfastest | Reads/s (M) | Max RSS (MB) | User CPU (s) | CPU % |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 39.55 | 1.00x | 1.08 | 29 | 261.7 | 669 |
| bbduk | 67.91 | 1.72x | 0.63 | 1,049 | 508.0 | 763 |
| fastp | 77.00 | 1.95x | 0.56 | 1,138 | 565.2 | 740 |
| cutadapt | 81.55 | 2.06x | 0.52 | 59 | 607.5 | 769 |
| trim-galore-rs | 118.96 | 3.01x | 0.36 | 77 | 924.1 | 780 |
| fastp-nfcore | 138.73 | 3.51x | 0.31 | 1,197 | 533.6 | 390 |
| adapterremoval | 262.86 | 6.65x | 0.16 | 58 | 440.2 | 168 |
| trim-galore | 435.82 | 11.02x | 0.10 | 36 | 2141.1 | 500 |
| trimmomatic | 438.56 | 11.09x | 0.10 | 2,878 | 2158.6 | 496 |

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| fastp | 0.6515 | 0.1059 | 0.96153 | 8,794 | 1,581,756 | 935 | 107 |
| fastp-nfcore | 0.6515 | 0.1059 | 0.96153 | 8,794 | 1,581,756 | 935 | 107 |
| **chelae** | 0.8030 | 0.1500 | 0.95270 | 2,222 | 1,966,396 | 796 | 4 |
| cutadapt | 0.9642 | 0.1711 | 0.95081 | 926 | 2,048,529 | 383 | 231 |
| trim-galore | 0.9642 | 0.1711 | 0.95081 | 926 | 2,048,529 | 383 | 231 |
| trim-galore-rs | 0.9642 | 0.1711 | 0.95081 | 926 | 2,048,529 | 383 | 231 |
| adapterremoval | 0.9684 | 0.1670 | 0.95208 | 88,075 | 1,853,895 | 53,573 | 0 |
| trimmomatic | 1.2216 | 0.2294 | 0.94159 | 14,733 | 2,421,546 | 6,868 | 1 |
| bbduk | 1.6225 | 0.2288 | 0.94782 | 122,172 | 2,005,582 | 12,899 | 37,458 |

### Dataset 11 - se-mirna-30bp-1x76

SE 1x76, insert 30 +/- 2 bp, error 0.1-1%, small-RNA adapters, 84.5 M reads.

#### Runtime (wgs config, 8 threads)

| Tool | Wall (s) | xfastest | Reads/s (M) | Max RSS (MB) | User CPU (s) | CPU % |
|---|---:|---:|---:|---:|---:|---:|
| **chelae** | 38.15 | 1.00x | 2.21 | 24 | 141.2 | 375 |
| bbduk | 43.33 | 1.14x | 1.95 | 1,361 | 312.3 | 739 |
| fastp | 61.17 | 1.60x | 1.38 | 1,108 | 368.3 | 611 |
| fastp-nfcore | 76.38 | 2.00x | 1.11 | 1,125 | 341.8 | 454 |
| cutadapt | 77.31 | 2.03x | 1.09 | 42 | 557.8 | 752 |
| trim-galore-rs | 100.42 | 2.63x | 0.84 | 39 | 767.1 | 765 |
| adapterremoval | 122.93 | 3.22x | 0.69 | 30 | 424.8 | 348 |
| trimmomatic | 169.70 | 4.45x | 0.50 | 3,390 | 971.1 | 575 |
| trim-galore | 180.44 | 4.73x | 0.47 | 34 | 1199.5 | 679 |

#### Accuracy (adapter_only config, 8 threads)

| Tool | RMSE | MAE | Exact match | FP | FN | Over-trim | Under-trim |
|---|---:|---:|---:|---:|---:|---:|---:|
| adapterremoval | 0.0000 | 0.0000 | 0.59882 | 0 | 0 | 0 | 0 |
| trimmomatic | 0.0238 | 0.0000 | 0.59882 | 0 | 14 | 0 | 0 |
| **chelae** | 0.0624 | 0.0001 | 0.59882 | 0 | 95 | 0 | 0 |
| fastp | 0.0624 | 0.0001 | 0.59882 | 0 | 95 | 0 | 0 |
| fastp-nfcore | 0.0624 | 0.0001 | 0.59882 | 0 | 95 | 0 | 0 |
| cutadapt | 0.0625 | 0.0001 | 0.59882 | 0 | 95 | 5 | 24 |
| trim-galore | 0.0625 | 0.0001 | 0.59882 | 0 | 95 | 5 | 24 |
| trim-galore-rs | 0.0625 | 0.0001 | 0.59882 | 0 | 95 | 5 | 24 |
| bbduk | 0.8616 | 0.0161 | 0.59869 | 0 | 17,745 | 0 | 0 |

## Throughput / scalability suite

Two PE datasets at 10× depth (~107 M pairs each), four fast-tier tools, two
replicates per cell, threads {4, 8}. 

### Dataset T1 - pe-wgs-150

PE 2×150, insert 350 ± 80 bp, 10× depth (~107.0 M PE pairs), TruSeq adapters.

#### Wall time vs thread count (s, median over 2 replicates)

| Tool             | 4 threads | 8 threads | Scaling 4→8 |
|------------------|----------:|----------:|------------:|
| **chelae**       |    301.31 |    212.34 |      1.42×  |
| trim-galore-rs   |    332.39 |    238.93 |      1.39×  |
| cutadapt         |    356.76 |    350.46 |      1.02×  |
| fastp            |   1103.74 |    621.30 |      1.78×  |

#### Resource use at 8 threads (median over 2 replicates)

| Tool             | Wall (s) | ×fastest | Reads/s (M) | Max RSS (MB) | User CPU (s) | CPU % |
|------------------|---------:|---------:|------------:|-------------:|-------------:|------:|
| **chelae**       |   212.34 |    1.00× |        1.01 |           57 |       1578.0 |   751 |
| trim-galore-rs   |   238.93 |    1.13× |        0.90 |           98 |       1780.7 |   754 |
| cutadapt         |   350.46 |    1.65× |        0.61 |           79 |       2161.2 |   653 |
| fastp            |   621.30 |    2.93× |        0.34 |        1,298 |       4773.9 |   774 |

### Dataset T2 - pe-cfdna-150

PE 2×150, insert 170 ± 40 bp, 10× depth (~107.0 M PE pairs), TruSeq adapters.

#### Wall time vs thread count (s, median over 2 replicates)

| Tool             | 4 threads | 8 threads | Scaling 4→8 |
|------------------|----------:|----------:|------------:|
| **chelae**       |    276.89 |    201.37 |      1.38×  |
| trim-galore-rs   |    484.83 |    324.98 |      1.49×  |
| cutadapt         |    383.20 |    371.18 |      1.03×  |
| fastp            |    862.38 |    497.75 |      1.73×  |

#### Resource use at 8 threads (median over 2 replicates)

| Tool             | Wall (s) | ×fastest | Reads/s (M) | Max RSS (MB) | User CPU (s) | CPU % |
|------------------|---------:|---------:|------------:|-------------:|-------------:|------:|
| **chelae**       |   201.37 |    1.00× |        1.06 |           65 |       1476.6 |   745 |
| trim-galore-rs   |   324.98 |    1.61× |        0.66 |          148 |       2529.7 |   784 |
| cutadapt         |   371.18 |    1.84× |        0.58 |           79 |       2419.5 |   684 |
| fastp            |   497.75 |    2.47× |        0.43 |        1,257 |       3782.8 |   766 |


