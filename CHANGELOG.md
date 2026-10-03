# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

When `cargo release` runs it promotes the `[Unreleased]` section below to a
versioned entry stamped with the release date; new entries should go under
`[Unreleased]` between releases.

## [Unreleased]

### Added

- `chelae trim --umi-tag <TAG>` writes the UMI extracted by `M` read-structure
  segments to the FASTQ comment as a SAM tag (e.g. `RX:Z:AAA-GGG`) instead of
  appending it to the read name, so `bwa mem -C` carries it onto every
  alignment with no FASTQ to uBAM round trip. `--umi-qual-tag <TAG>` adds the
  UMI base qualities (e.g. `QX:Z:FFF FFF`, segments joined by a space as the
  SAM specification recommends). The comment is rewritten to hold only SAM
  tags: existing `TAG:TYPE:VALUE` fields are kept, a Casava 1.8 index becomes
  `BC:Z:` as with `samtools import -i`, and other text is dropped, since
  `bwa mem -C` would otherwise copy it in as an invalid SAM field. Because
  that drops a Casava `1:N:0`/`2:N:0` read number, interleaved output then ends
  the read-ids of mates that share one in `/1` and `/2`, so tools that pair on
  the suffix, such as `samtools import -s`, still tell R1 from R2; `bwa mem`
  strips it. A `/1` or `/2` already on the input names is not doubled.

### Fixed

- `chelae trim` now appends a read-structure UMI to the read-id before a
  trailing `/1` or `/2` mate suffix (`@frag/1` becomes `@frag:AAAA/1`, not
  `@frag/1:AAAA`), so once `bwa mem` or fgumi strips the suffix both mates
  have the same name ending in the UMI field. The read-id now also ends at a
  tab, so a tab-separated comment (`@frag<TAB>BC:Z:ACGT`) no longer gets the
  UMI appended after it or has its colons counted as read-id fields.

## [0.2.1] - 2026-10-03

### Changed

- When a pair reads through and has read-structures, `chelae trim` now also
  trims from each read's 3' end the mate's fixed segments before its first
  template, so a short insert no longer keeps the reverse complement of the
  mate's UMI and skip bases. This includes inserts only a few bases longer than
  the reads, which leave no adapter in either read: when a mate has such a
  prefix, the overlap search also covers them, with or without
  `--insert-size-stats`. Pairs whose insert is longer than the reads plus the
  mate's prefix are unchanged; the removed bases are counted under
  `bases_trimmed_read_structure`.
- `chelae trim --expected-insert-size` now fixes where the paired-end overlap search starts for the whole run, rather than seeding an estimate that each worker thread went on adjusting. Without it, the insert size is estimated once, from the overlaps detected in the first 65,536 pairs, and used from then on.
- `chelae detect` likewise fixes its overlap-search start once it has detected 64 overlaps, instead of adjusting it as it goes. Its output was already deterministic; the tails it harvests from a few repeat pairs may differ from 0.2.0.

### Removed

- `chelae trim`'s hidden `--batch-size` option; batches are always 1,024 records.

### Fixed

- `chelae trim` now trims read-through pairs in which one mate is too degraded
  to read its adapter or to align, as low-complexity or poly-G 3' ends often
  are. The overlap is accepted when either mate's tail alone is unlikely to
  match adapter by chance, rather than only when both tails together are. When
  no overlap is found and a mate's last 30 bases are at least half below Q20, a
  match of at least 20 adapter bases at most 10% mismatched anywhere in either
  mate (with `--kit`, `--adapter-sequence` or `--adapter-fasta`) marks the
  insert's end and both mates are cut there. Such pairs previously kept their
  adapter.

- `chelae trim` metrics: when a read was too short for its read-structure and
  the pair was dropped, bases that the read-structure step had already cut from
  the pair were counted nowhere, so `bases_in` exceeded `bases_out` plus the
  `bases_trimmed_*` and `bases_filtered` counts. They are now counted under
  `bases_trimmed_read_structure`.
- `chelae trim` and `chelae detect` no longer hang when their two `--inputs` are named pipes (FIFOs) fed by one process, such as `k2tools filter`, that opens both before writing or writes ahead on one of them. Either input can now run at least 1 MB ahead of the other.
- `chelae trim` output is now fully deterministic: the same input and options give byte-identical output on every run and at every `--threads` value, including the `--insert-size-stats` histogram. The remaining cases where a pair's result depended on thread scheduling (post-cut tails that look like adapter at two overlap shifts, and long tandem repeats under `--insert-size-stats`) came from each worker thread starting its overlap search from its own running insert-size estimate. The starting point is now the same for every thread (see above). A few such pairs may be trimmed differently than in 0.2.0.
- `chelae trim`'s progress line (every 5M records) is logged again. Batches are 1,024 records, so the running count was almost never an exact multiple of 5,000,000 and the line appeared only every 80M records.

## [0.2.0] - 2026-09-25

### Added

- `chelae trim` and `chelae detect` gain interleaved paired-end FASTQ I/O,
  stdin/stdout, and (for `trim`) uncompressed output:
  - **Interleaved PE** input and/or output, inferred from input/output counts
    with no new flag: a single input is sniffed for an interleaved pair by
    peeking up to its first 4 records and selecting a mate-naming convention
    (Casava 1.8+ `1:`/`2:` comment markers, or ENA-style `/1`/`/2` at the
    end of the comment's original read name; identical names, as in SRA's
    default `fastq-dump`/`fasterq-dump` defline; a trailing `/1`/`/2`; or a
    trailing `.1`/`.2`/`_1`/`_2`, as in `fastq-dump -I --split-spot` output)
    confirmed across two probe pairs, so a standard SE SRA file (`@SRR.1`,
    `@SRR.2`, `@SRR.3`, …) doesn't misdetect as interleaved. The same
    convention then enforces pairing and mate orientation for the rest of
    the run: a reversed `/2`-then-`/1` pair is rejected, including in the
    first records, rather than read as single-end; an
    out-of-sync, odd-length, or non-corresponding stream fails loudly,
    naming the offending pair and its file-record indices. A single output
    interleaves both mates. Two `--inputs` files are always split R1/R2 by
    position and never sniffed for interleaving.
  - **stdin/stdout via `-`**: `-i`/`-o` on `trim` and `-i` on `detect` are
    now optional and default to `-`. Reading FASTQ from an interactive
    terminal is refused with an actionable error; writing to one is always
    allowed (e.g. `chelae trim | head`). If a downstream reader closes
    stdout early, chelae stops promptly and exits successfully with
    whatever partial output it had produced, rather than erroring (this also
    covers `chelae detect -o -`). With `--metrics`/`--json`, a warning at
    that moment notes that their counts may include reads chelae processed
    but that never made it out before the pipe closed.
  - **`chelae trim --output-compression {auto,bgzf,none}`** (default `auto`):
    `auto` writes BGZF for a `.gz`/`.bgz`-suffixed path (case-insensitive)
    and plain text otherwise; `bgzf`/`none` force the encoding on every
    output regardless of extension. `--compression-level` applies only to
    BGZF outputs; setting it when every output is plain text logs a warning.
  - Input gzip/BGZF detection switched from file-extension to magic-byte
    sniffing (required for stdin; also fixes misnamed files).
  - `chelae trim` rejects two outputs (or an output and `--metrics`/`--json`)
    that name the same file, and `trim`/`detect` reject two `--inputs` that
    name the same file; paths are compared after resolving symlinks, `.`/`..`
    and relative components. `chelae trim` also rejects `--metrics -`/`--json
    -` (neither ever had a `-` default; passing `-` previously created a
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
- `chelae trim`'s default `--compression-level` is now 1 (was 5). On
  simulated 2×100–2×250 and single-end data, level 1 used ~60% less CPU and
  wall time than level 5 for BGZF output that was 3–5% larger. Pass `-c 5`
  for the previous default.
- `chelae trim`'s paired-end overlap search and sequence-based adapter
  matching now screen 16 candidate offsets at a time with SIMD and fully
  compare only those that could match, with identical output. On simulated
  data this cut total CPU at compression level 1 by 5–11% for paired-end
  2×100–2×250 (most with long reads or inserts) and by 9–23% for
  single-end (most with several adapters, e.g. `--kit all`); with
  uncompressed output, by 6–29% and 36–40%.
- `chelae trim` spends less CPU moving records around, with identical
  output: each batch holds its records in one list rather than one list per
  pair, input is no longer copied through a redundant buffer, the paired
  read-name check and per-read base statistics are cheaper, and
  uncompressed-output buffers no longer regrow every batch. On simulated
  data this cut CPU cycles by 3–4% at compression level 1, and with
  uncompressed output by 17–21% for paired-end and 10% for single-end.
- The performance figures above were measured on Apple M4 (aarch64). On
  x86 (an AMD EPYC Zen 5 EC2 instance, `-t 8 -c 1`), `chelae trim` 0.2.0
  uses 11–26% less CPU than 0.1.0 on the same simulated paired-end and
  single-end data, saving the most with long reads and inserts, and its
  output is identical to that of aarch64 builds.
- `chelae trim`'s split paired-end input (two `--inputs` files) now has its
  read names checked pair-by-pair, using the same mate-naming conventions as
  interleaved input: once the first pair establishes a convention, a later
  pair whose names don't correspond fails the run, naming the offending
  record. A first pair in mate-2/mate-1 order (swapped `--inputs`), or
  whose names are both mate-marked but don't correspond, also fails. If the
  first pair's names aren't both mate-marked in a recognized way, a warning
  is logged and records are paired by position only, as in 0.1.0 (which only
  checked that both files had the same number of records).
- `chelae` now declares Rust 1.89 as its minimum supported version
  (`rust-version` in `Cargo.toml`), so building it with an older toolchain
  fails up front with a clear error.

### Fixed

- `chelae trim` paired-end output is now identical across runs and thread
  counts, except in the rare pair whose post-cut tails look convincingly
  like adapter at two different overlap shifts. With `--insert-size-stats`,
  the insert-size histogram can also still vary for tandem repeats longer
  than the reads; the trimmed reads don't. In tandem repeats (e.g.
  satellites, telomeres) R1 and R2 can overlap acceptably at several
  shifts, and each worker thread's search started from its own running
  insert-size estimate, so which shift won depended on thread scheduling.
  A first-found overlap is now kept only if its post-cut tails are unlikely
  to match adapter by chance; otherwise every shift is evaluated and the
  best kept, preferring one whose tails pass that test, then the
  best-aligned. On simulated 2×150 data with ~145 bp inserts this also cuts
  over-trimmed reads by ~17%, for ~2% more CPU.
- `chelae trim` rejects more than two `--inputs` or `--outputs` given across
  repeated flags (e.g. `-i a.fq b.fq -i c.fq`); clap's per-flag limit
  didn't catch the extra paths.
- A truncated or failing input (e.g. a cut-off gzip file, or an upstream
  process dying mid-stream) no longer prints a spurious parser panic ahead
  of the real `FASTQ read error`.
- `chelae trim --expected-insert-size` is now honored. The hint is stored
  in I-space (insert size) rather than shift-space, so it takes effect on
  the first pair regardless of variable read length, and `--insert-size-stats`
  reports anchor against the same I-space estimate.

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

[Unreleased]: https://github.com/fulcrumgenomics/chelae/compare/v0.2.1...HEAD
[0.2.1]: https://github.com/fulcrumgenomics/chelae/releases/tag/v0.2.1
[0.2.0]: https://github.com/fulcrumgenomics/chelae/releases/tag/v0.2.0
[0.1.0]: https://github.com/fulcrumgenomics/chelae/releases/tag/v0.1.0
