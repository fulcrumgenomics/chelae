# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

When `cargo release` runs it promotes the `[Unreleased]` section below to a
versioned entry stamped with the release date; new entries should go under
`[Unreleased]` between releases.

## [Unreleased]

### Added

- `chelae detect` subcommand: identifies the 3' adapter sequence(s) present
  in one or two FASTQ files by sampling a modest number of records.
  - Paired-end input discovers adapters via R1/R2 overlap detection (no kit
    knowledge required) and builds a position-by-position consensus across
    reads that landed in the same near-identical k-mer cluster.
  - Single-end input scores each read against every built-in kit adapter
    plus any user `--adapter-sequence` / `--adapter-fasta` candidate.
  - 3'-end cleanups (poly-G default-on at min run 10; poly-X A/C/T default-on
    at min run 5; cut-right quality trim default-on at `4:20`) are applied
    to every read before the overlap probe / candidate scan, so 2-color
    chemistry artifacts and quality-degraded tails don't corrupt detection.
    All three can be tuned or disabled (`0` for the homopolymer flags,
    `off`/`none`/`no` for the quality flag).
  - Sampling stops once `--num-detections` usable detections (default 5000)
    or `--max-reads` records (default 1M) are reached. A
    `--min-detections-for-report` floor (default 20) refuses to report on
    samples too small to be informative.
  - Console report is two sections: matched kit(s) with their published
    adapter sequences, and the full-length discovered consensus per mate
    with **uppercase** marking the kit-stable region and **lowercase**
    marking any per-sample extension (typically an i7/i5 barcode tail).
  - Optional `--output-fasta` writes a cross-sample-portable FASTA: when
    a discovered sequence matches a known kit (within 1 mismatch over the
    first 16 bp), the kit's full published adapter is emitted in place of
    the sample-specific consensus, so the FASTA round-trips cleanly through
    `chelae trim --adapter-fasta` on every sample in a batch.

### Fixed

- `chelae trim --expected-insert-size` is now honored. The hint is stored
  in I-space (insert size) rather than shift-space, so it takes effect on
  the first pair regardless of variable read length, and `--insert-size-stats`
  reports anchor against the same I-space estimate.

### Internal

- `BUFFER_SIZE` and the FASTQ-reader opening helper are now shared between
  `chelae trim` and `chelae detect` via `commands/utils.rs`.
- The benchmark pipeline gained a `trim-galore-rs` tool wrapper.

## [0.1.0] - 2026-05-13

### Added

- Initial public release of `chelae`.
- `chelae trim` subcommand: single-pass short-read FASTQ trimming and
  filtering for SE and PE input.
  - Poly-G 3' trim (on by default).
  - Paired-end overlap-based adapter detection (on by default for PE),
    plus 3' adapter trimming by user-supplied sequence (`--adapter-sequence`),
    FASTA (`--adapter-fasta`), and built-in kit presets (`--kit`: `truseq`,
    `nextera`, `small-rna`, `aviti`, `mgi`/`dnbseq`, and `all`).
  - Read-structure based hard-trim and UMI extraction (`--read-structures`),
    applied after adapter trim so tail-skip segments operate on the cleaned
    template.
  - Optional poly-X 3' trim (`--trim-polyx`).
  - Optional 5'→3' and 3'→5' sliding-window quality trimming.
  - Length filter, N-base filter, mean-quality filter, and low-quality-fraction
    filter.
  - BGZF-compressed FASTQ output.
  - Optional fastp-compatible JSON report (`--json`) consumed by MultiQC's
    `fastp` module unchanged.
  - Optional paired-end insert-size histogram (`--insert-size-stats`),
    fastp-shape.
- Portable SIMD kernels (via the `wide` crate) for adapter detection,
  overlap detection, poly-X trim, and sliding-window quality trim.
- x86_64 release binaries packaged via `cargo-multivers` with three CPU
  variants (`x86-64`, `x86-64-v2`, `x86-64-v4`) and runtime dispatch.
- aarch64 release binaries as a single generic ARMv8-A / NEON build.

### Origin

`chelae` was extracted from [`fqtk`](https://github.com/fulcrumgenomics/fqtk)
on 2026-04-21; the entire `chelae trim` implementation was developed as
`fqtk trim` on the `tf_trim` branch and split into this standalone crate.
The pre-split incremental history (design decisions, performance work,
benchmarks) lives in the `fqtk` repo.

[Unreleased]: https://github.com/fulcrumgenomics/chelae/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/fulcrumgenomics/chelae/releases/tag/v0.1.0
