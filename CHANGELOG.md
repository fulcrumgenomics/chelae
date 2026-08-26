# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

When `cargo release` runs it promotes the `[Unreleased]` section below to a
versioned entry stamped with the release date; new entries should go under
`[Unreleased]` between releases.

## [Unreleased]

### Added

- `chelae trim` and `chelae detect` gain interleaved paired-end FASTQ I/O,
  stdin/stdout, and (for `trim`) uncompressed output:
  - **Interleaved PE** input and/or output, inferred from input/output counts
    with no new flag: a single input is sniffed for an interleaved pair by
    peeking up to its first 4 records and selecting a `PairingRule` (Casava
    1.8+ comment markers or bare-name equality; a trailing `/1`/`/2`; or a
    trailing `.1`/`.2`/`_1`/`_2` — covering SRA `fastq-dump -I --split-spot`
    output) confirmed across two probe pairs, so a standard SE SRA file
    (`@SRR.1`, `@SRR.2`, `@SRR.3`, …) doesn't misdetect as interleaved. The
    same rule then enforces mate orientation (a reversed `/2`-then-`/1` pair
    is rejected) for the rest of the run; a single output interleaves both
    mates. Two `--inputs` files are always split R1/R2 by position and never
    sniffed for interleaving — but their read names are now checked to
    correspond pair-by-pair, using the same rule selection as the
    interleaved case. An out-of-sync, odd-length, or non-corresponding input
    fails loudly, naming the offending pair (and, for interleaved input, the
    underlying file-record indices).
  - **stdin/stdout via `-`**: `-i`/`-o` on `trim` and `-i` on `detect` are
    now optional and default to `-`. Reading FASTQ from an interactive
    terminal is refused with an actionable error; writing to one is always
    allowed (e.g. `chelae trim | head`). If a downstream reader closes
    stdout early, chelae stops promptly and exits successfully with
    whatever partial output it had produced, rather than erroring.
  - **`chelae trim --output-compression {auto,bgzf,none}`** (default `auto`):
    `auto` writes BGZF for a `.gz`/`.bgz`-suffixed path (case-insensitive)
    and plain text otherwise; `bgzf`/`none` force the encoding on every
    output regardless of extension.
  - Input gzip/BGZF detection switched from file-extension to magic-byte
    sniffing (required for stdin; also fixes misnamed files).
  - `chelae trim` rejects two outputs (or an output and `--metrics`/`--json`)
    that resolve to the same path, and rejects `--metrics -`/`--json -`
    (neither ever had a `-` default; passing `-` previously created a
    literal file named `-`).
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
  - Per-position consensus uses a discontinuity-aware cut: a running-min
    baseline of the per-column majority fraction stops the consensus at
    the first sharp drop, at an absolute floor of 50% majority, or at
    the column coverage floor — whichever comes first. This prevents
    plurality-noise bases (e.g. an imbalanced sample-index pool where one
    index dominates at 25%) from being emitted as if they were real.
  - Optional `--output-fasta` writes a cross-sample-portable FASTA: when
    a discovered sequence matches a known kit (within 1 mismatch over the
    first 16 bp), the kit's full published adapter is emitted in place of
    the sample-specific consensus, so the FASTA round-trips cleanly through
    `chelae trim --adapter-fasta` on every sample in a batch. `.gz`-
    extensioned output paths are transparently gzip-compressed. If any
    mate (PE) or candidate (SE) fails to clear `--min-fraction`, detect
    hard-fails with an actionable error rather than writing a silently
    incomplete FASTA.

### Changed

- **Breaking**: `chelae trim -o out.fq` (no `.gz` suffix) now writes plain
  text instead of silently writing BGZF-compressed bytes to a misleadingly-
  named file. Pass `--output-compression bgzf` to force BGZF on any path, or
  name the output `*.gz` for the previous default behavior.
- `chelae trim -i`/`-o` and `chelae detect -i` no longer require an argument
  (`required = true` dropped); both default to `-` (stdin/stdout) when
  omitted.

### Fixed

- `chelae trim --expected-insert-size` is now honored. The hint is stored
  in I-space (insert size) rather than shift-space, so it takes effect on
  the first pair regardless of variable read length, and `--insert-size-stats`
  reports anchor against the same I-space estimate.
- Gzip content-sniffing now reads a full 2-byte magic-number probe (looping
  until 2 bytes are buffered or EOF) instead of a single `fill_buf()` call,
  which could legally return just 1 byte on a pipe and misdetect a gzipped
  stdin stream as plain text; I/O errors during the probe now propagate
  instead of being silently treated as "not gzip".
- A completely empty lone input (e.g. `printf '' | chelae trim -o r1.fq -o
  r2.fq`) no longer forces single-end layout regardless of what the rest of
  the CLI implies; the mate count is inferred from `--outputs` /
  `--read-structures` / `--adapter-sequence` instead, so a genuinely paired
  empty run succeeds with valid (empty) outputs on every requested file.
- A reversed interleaved pair (mate 2's record immediately followed by mate
  1's) is now rejected instead of silently accepted as a valid pair.

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

### Fixed

- Usage line now always reads `Usage: chelae ...`. clap derives the displayed
  program name from `argv[0]`'s basename when `bin_name` is unset; the
  cargo-multivers fexecve/memfd launcher passes `/proc/self/fd/N` as `argv[0]`
  under binfmt emulation (e.g. an amd64 biocontainer running on Apple Silicon),
  so the usage line printed `Usage: 11 ...` in that environment. Fixed by
  setting `name`/`bin_name = "chelae"` explicitly on the CLI. Also pinned
  `cargo-multivers` to `>=0.12.0` (the release carrying the fexecve/memfd fix)
  in the benchmark pipeline installer. Mirrors fulcrumgenomics/riker#38.

### Origin

`chelae` was extracted from [`fqtk`](https://github.com/fulcrumgenomics/fqtk)
on 2026-04-21; the entire `chelae trim` implementation was developed as
`fqtk trim` on the `tf_trim` branch and split into this standalone crate.
The pre-split incremental history (design decisions, performance work,
benchmarks) lives in the `fqtk` repo.

[Unreleased]: https://github.com/fulcrumgenomics/chelae/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/fulcrumgenomics/chelae/releases/tag/v0.1.0
