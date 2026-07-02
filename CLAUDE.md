# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`chelae` is a FASTQ trimming and filtering toolkit written in Rust. The name is the plural of *chela*, the pincer-like claws of crustaceans — apt for a tool that clips reads.

The crate ships two subcommands:

- **`chelae trim`** runs the common short-read preprocessing pipeline (poly-G trim → adapter trim → read-structure hard-trim + UMI extraction → poly-X trim → sliding-window quality trim → length/N-base/mean-quality/low-qual-fraction filters) on SE or PE FASTQ data. Outputs are BGZF-compressed FASTQ plus a fastp-compatible JSON report for MultiQC consumption.
- **`chelae detect`** identifies the adapter sequence(s) present in a FASTQ. PE input discovers adapters by harvesting the post-template tails from overlap-detected pairs; SE input scores each read against every built-in kit (plus optional user `--adapter-sequence` / `--adapter-fasta` entries). Reports the discovered sequences and the kit (if any) they match; optional FASTA output is designed to feed straight back into `chelae trim --adapter-fasta`.

Read-structure (in `trim`) runs after adapter trim so tail-skip segments (like the `10S` in `+T10S`) drop bases from the post-adapter template, not from the raw read where the adapter step would usually have removed them already. If a read is shorter than the read-structure's fixed-length segments require (either originally or after adapter trim), the pair is dropped and counted under `reads_filtered_length`.

The crate structure (multi-subcommand CLI with `enum_dispatch`) was deliberately built open-ended so additional trimming/filtering/QC utilities can land as siblings of `trim` and `detect` without a CLI shape change.

## Origin — very important

**`chelae` was extracted from [`fqtk`](https://github.com/fulcrumgenomics/fqtk) on 2026-04-21.** The entire `chelae trim` implementation — code, design decisions, performance work, tests — was developed as `fqtk trim` on the `tf_trim` branch in that repo and then split into this standalone crate. If you are working on `chelae` and need context that isn't in this file:

1. The full design history is in the [`fqtk`](https://github.com/fulcrumgenomics/fqtk) repo's git log and GitHub PR history, branch `tf_trim`.
2. Tim keeps a local checkout at `/Users/tfenne/work/open-source/fqtk` (on branch `tf_trim`) — all the incremental commits that became the initial chelae import are visible there (pre-squash, on the `tf_trim_backup` branch locally), along with their commit messages. If you need to understand *why* a design choice was made, that's the source of truth.
3. Auto-memory from fqtk development was pre-seeded into this project's memory directory at the time of extraction. See `~/.claude/projects/-Users-tfenne-work-open-source-chelae/memory/` for feedback rules, project notes, and references.

The `chelae` initial import is a three-commit squash of 30+ fqtk commits. The squashed history omits detail intentionally — when you need "why did we do it this way", look in the fqtk backup branch, not here.

## Build & Test Commands

```bash
# Full verification (format + clippy + tests) — run before pushing.
bash ci/check.sh

# Individual steps
cargo fmt --all
cargo clippy --all-features --all-targets -- -D warnings
cargo test

# Run a single test
cargo test <test_name>

# Build release binary
cargo build --release
```

CI (`.github/workflows/build_and_test.yml`) additionally runs `src/scripts/precommit.sh` (the same checks with `--locked`).

The README has a hand-curated **Options** table that summarizes every `chelae trim` flag. When you add, remove, rename, or materially change a CLI option, update that table in `README.md` to match. The `chelae trim --help` output remains the authoritative reference; the README table is the short-form pointer.

## Toolchain

Pinned to Rust 1.95 via `rust-toolchain.toml`. Format settings: `max_width = 100`, `use_small_heuristics = "max"` (see `rustfmt.toml`).

## Build Targeting

### x86_64: cargo-multivers

x86_64 release binaries are packaged via [`cargo-multivers`](https://github.com/ronnychevalier/cargo-multivers) into a single launcher that embeds three CPU-specific builds and dispatches to the best match at startup. See `[package.metadata.multivers.x86_64]` in `Cargo.toml`:

```toml
cpus = ["x86-64", "x86-64-v2", "x86-64-v4"]
```

- `x86-64` — SSE2 baseline, any 64-bit x86 (2003+)
- `x86-64-v2` — SSE4.2 + POPCNT (2008+). Captures nearly all of the scalar codegen win
- `x86-64-v4` — AVX-512F/BW/CD/DQ/VL (2017+ server / 2022+ consumer). ~1% additional win

We intentionally skip `x86-64-v3`: our Granite Rapids benchmark showed v2 and v3 within measurement noise on chelae's workload. The historical "x86-64-v3 wins 6% over baseline" finding is actually a v1→v2 win; v2→v3 contributes ~0. Including v3 would bloat the binary without buying anything.

Variants are delta-compressed (`gdelta`) + lz4. Total binary is ~3.7 MB (vs ~2.9 MB for a single-variant build). Startup adds ~0.2 s for decompression + `memfd_create + exec` — negligible for chelae's batch workload.

The cargo-multivers runner sorts variants by feature count descending and picks the first match — so v4 runs on capable hardware, falling back to v2 on pre-AVX-512 systems and v1 on pre-SSE4.2 systems.

### Build commands

- x86_64: `cargo multivers --profile dist`
- aarch64: `cargo build --profile dist`

`.cargo/config.toml` carries no build-wide rustflags. Plain `cargo build --release` therefore produces a portable baseline (target-cpu = the target triple's default — `x86-64` v1 on x86_64). `cargo multivers` reads `[package.metadata.multivers.x86_64].cpus` from `Cargo.toml` and supplies the per-variant `target-cpu` itself, so no file-swapping is needed.

For local profiling on x86_64 where you want `target-cpu=native`, set `RUSTFLAGS="-C target-cpu=native"` on the command line. The repo doesn't make that the default because the primary dev hosts are aarch64 (Apple Silicon), where the target-triple default is already current-generation.

### Release profile

`[profile.dist]` in `Cargo.toml` inherits `release` with `incremental = false` for deterministic delta compression across multivers variants. Use `cargo multivers --profile dist` for x86_64 releases and `cargo build --profile dist` for aarch64 releases.

### aarch64

Single binary, no multivers. Benchmarks showed Neoverse-specific tuning yields only ~1-2% over generic ARMv8-A with near-zero cross-tuning penalty, so multivers infrastructure isn't justified. The generational upgrade (Graviton3 → Graviton4 = +24%) dwarfs any tuning delta anyway. Release build target-cpu is whatever we land on post-benchmarking; the dev default (`target-cpu=native`) is fine locally.

## Architecture

The crate produces both a library (`chelae_lib`) and a binary (`chelae`).

### Library (`src/lib/`)

- `mod.rs` — Exposes `IUPAC_MASKS`, a 256-entry ASCII-indexed LUT mapping base characters (ACGT + IUPAC ambiguity codes, case-insensitive) to 4-bit masks. Two bases are IUPAC-compatible iff `IUPAC_MASKS[a] & IUPAC_MASKS[b] != 0`. The adapter matcher's IUPAC slow path uses this.
- `adapter_db.rs` — `KitAdapter` preset table (`truseq`, `nextera`, `small-rna`, `aviti`, `mgi` / `dnbseq` alias) + the `all` union. Exposed to `trim` via `--kit`.

### Binary (`src/bin/`)

- `main.rs` — CLI entry point using clap derive. `Subcommand` enum is `enum_dispatch`-wired to the `Command` trait. `mimalloc` is the global allocator. `ensure_avx2_or_die` runs pre-argv on x86_64.
- `commands/command.rs` — `Command` trait that each subcommand implements.
- `commands/utils.rs` — Small shared utilities. Currently only `fmt_count` (comma-separated u64 formatter). Keep this as the home for cross-subcommand helpers.
- `commands/trim.rs` — The `Trim` subcommand. Large file (~5700 lines) organized into sections per the conventions below. Key types (roughly top-down in the hot path):
  - `Trim` — the clap-derived CLI options struct; `impl Command for Trim` holds `execute()`.
  - `Pipeline<'a>` — per-worker mutable state (scratch buffers, per-worker counters). Workers hold a `Pipeline` and call its methods for each batch.
  - `Batch` / `Packet` — the unit passed reader → worker → writer. Each `Packet` carries a `oneshot::Sender` per output file; workers `.send()` serialized+compressed bytes on those senders in their input order, and the writer thread consumes them in submit order to preserve FASTQ record order.
  - `OverlapStats` / `detect_pe_overlap` — PE-overlap evidence acceptance via a single signed-shift outward walk centered on the running mean detected insert; `try_shift_neg` / `try_shift_pos` split the probe by sign so the stats-off hot path skips a runtime branch. With `--insert-size-stats`, the walk extends to positive shifts (the I > R inner-overlap geometry) and emits a fastp-shape histogram in the JSON. 5% hysteresis on center updates.
  - `AdapterSet` / `Adapter` — the active adapter set (user sequences + kit union + FASTA). `Adapter::pure_acgt` gates the SIMD fast path in `find_adapter_3prime`.
  - `FastpJsonReport` + supporting structs — fastp-shape JSON emitted when `--json` is set; key names mirror fastp so MultiQC's fastp module parses chelae's output unchanged.
  - Module-level functions hold SIMD kernels (`find_polyx_tail_len`, `trim_quality_sliding_3prime`, `count_mismatches_ci_bounded`, `observe_stats`, `count_bases_below_q`, …) and byte-slice helpers that don't belong to a type we own.

### Adding a New Subcommand

1. Create `src/bin/commands/<name>.rs` with a struct deriving `Parser` that implements `Command`.
2. Add the module to `src/bin/commands/mod.rs`.
3. Add the variant to the `Subcommand` enum in `src/bin/main.rs`.

## Code Organization

### Module Ordering

Command modules (`src/bin/commands/*.rs`) must follow this ordering convention:

1. **`use` statements**
2. **Constants and type aliases**
3. **Structs/enums, each immediately followed by all its impl blocks:**
   - CLI options / command struct (the `clap::Parser`-derived struct) + `impl` + `impl Command`
   - Per-run config structs (e.g. `PipelineConfig`) + `impl`
   - Metric / output-row structs + `impl`
   - Helper structs/enums — ordered higher-to-lower level (if A uses B, A comes first), then by importance to the overall implementation
4. **Module-level functions** (functions operating on primitives or external types — SIMD kernels, byte-slice helpers, parsers that don't belong to a type we own)
5. **`#[cfg(test)] mod tests`**

### Impl Block Rules

- Consolidate all inherent methods into one `impl` block per struct; keep each trait impl as a separate block.
- Trait impls go immediately after the struct's own `impl` block.
- Within an impl block, order methods callers-before-callees (higher in the call stack first); constructors (`new`, `from_str`) always come first.
- Functions that naturally belong to a type we own should be methods, not standalone functions.

### Keep trim-related code in trim.rs

While developing `trim` in fqtk, we kept all trim-specific code in the single `trim.rs` file rather than preemptively splitting into submodules. That rule continues here: do not split `trim.rs` into per-section files without a concrete maintenance win. Cross-subcommand helpers go in `src/bin/commands/utils.rs` or `src/lib/`.

## Performance Notes (established in fqtk, carried over intact)

The following baselines were established during fqtk development. They're summary reference points — the work and measurements happened in fqtk and are recorded in that repo's history.

- **PE benchmark** (27M pairs, 6 worker threads, compression-level 4, `--detect-adapter-for-pe`): ~63 s wall / ~412 s user CPU on Apple Silicon, ~66 s wall / ~485 s user CPU on EC2 Granite Rapids. fastp at matched config: ~2.4× slower.
- **SE benchmark with single adapter** (27M reads): fqtk ~47 s vs fastp ~101 s (2.14× faster).
- **SE benchmark with all 5 kit adapters** (same data): chelae ~50 s (+6% vs single-adapter) vs fastp ~309 s (+206% vs single-adapter). fastp's adapter matching scales linearly in adapter count; chelae's SIMD ACGT kernel amortizes it.
- **AVX2 via `x86-64-v3`**: ~6% wall-time recovery vs SSE2 baseline on Granite Rapids, no measurable cost on AVX2 hardware.
- **Multiversion runtime dispatch**: tried and rejected in fqtk (see memory `project_multiversion_rejected_2026_04.md`). The v3 gain is scattered auto-vectorized scalar code, not our SIMD kernels, so runtime dispatch wouldn't help.
- **`flate2` backend**: switched from `zlib-ng` (C dep) to `zlib-rs` (pure Rust) in the final fqtk commit. Ties or edges out zlib-ng in wall time on both Apple Silicon and EC2 x86; drops the C toolchain requirement.

## Open Work Carried Over From fqtk

These were the pending fqtk tasks at the time of split and remain open for chelae:

- **Reader-thread CPU cost vs thread count** — when `--threads` is small, the reader thread runs at near-100% CPU and may be the bottleneck; worth revisiting whether to budget a reader slot against `--threads` at low thread counts.
- **Profiling-guided optimizations beyond SIMD** — samply/perf-led; no concrete candidate yet.
- **Explore mgzip/mzip blocking as an alternative to BGZF** for output — potentially faster than BGZF at matched compression ratios.
- **SIMD-aware FASTQ parser** (replace or augment `seq_io`) — the parser is the last major scalar loop in the hot path. Speculative; my (Claude's) prior analysis was skeptical since `seq_io` already uses `memchr` for newlines.
- **SIMD FASTQ serializer** — ditto; byte-moving code is likely memory-bound, not compute-bound.

## Key Dependencies

- `read-structure` — parses read structure strings (e.g. `8B92T`) into typed segments.
- `fgoxide` — Fulcrum Genomics utilities (chunked read-ahead iterator, delimited I/O).
- `seq_io` — FASTQ parsing.
- `bgzf` — BGZF compression (via libdeflate under the hood).
- `flate2` with `zlib-rs` feature — pure-Rust gzip decode for input files. `zlib-ng` is a viable alternative but we switched to `zlib-rs` on 2026-04-21 to drop the C-compiler dependency; see the `build: switch flate2 backend` commit in the fqtk `tf_trim` branch for benchmark data.
- `wide` — portable SIMD primitives (u8x16, u8x32, i8x16, u16x8). Compile-time dispatch only on stable Rust — the `.cargo/config.toml` target-cpu controls which vector width LLVM selects.
- `crossbeam-channel` + `oneshot` — worker-pool plumbing. `oneshot` per-packet senders provide send-once ordering without a Mutex.
- `mimalloc` — global allocator.

## Things To Not Do

- **Don't split `trim.rs` into multiple files** without a specific win. Its size is a feature of the decision log.
- **Don't rename the fastp-compatible JSON keys.** MultiQC's fastp module consumes them by name.
- **Don't add runtime SIMD dispatch (`multiversion` crate or hand-rolled).** Tried and reverted in fqtk; the perf delta on AVX2 hosts comes from scalar auto-vectorization, not the SIMD kernels.
- **Don't reintroduce `zlib-ng`** unless you have a specific measurement showing it wins on a platform you care about — adding the C toolchain dependency back is a real cost.
- **Don't add comments explaining what code does.** Add comments explaining *why* only when the why is non-obvious. Most trim.rs comments meet this bar; the ones that don't are drift.
