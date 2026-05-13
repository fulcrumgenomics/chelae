# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

When `cargo release` runs it promotes the `[Unreleased]` section below to a
versioned entry stamped with the release date; new entries should go under
`[Unreleased]` between releases.

## [Unreleased]

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

[Unreleased]: https://github.com/fulcrumgenomics/chelae/compare/HEAD...HEAD
