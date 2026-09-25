# Contributing to chelae

This document is for developers working on `chelae` itself. End-user
documentation lives in [README.md](README.md).

## Pre-push checks

Before pushing, run the full verification suite:

```console
bash ci/check.sh
```

This runs the cargo aliases defined in `.cargo/config.toml`, which CI runs as
well: `cargo ci-fmt` (a formatting check; `cargo fmt --all` fixes it),
`cargo ci-lint` (clippy with `-D warnings`) and `cargo ci-test`, all `--locked`
where it applies.

The pinned toolchain (`rust-toolchain.toml`) is what CI uses; the minimum
supported version is the `rust-version` in `Cargo.toml`. Clippy's
`incompatible_msrv` lint flags standard-library APIs newer than that minimum.

## Building from source

Clone the repository and build in release mode. If you don't have Rust yet, install it with [rustup](https://rustup.rs/); `rust-toolchain.toml` pins the version CI uses, and rustup fetches it on the first build.

```console
git clone https://github.com/fulcrumgenomics/chelae.git
cd chelae
cargo build --release
./target/release/chelae --help
```

## Build targeting and portability

`cargo build --release` produces a portable build for the target's baseline CPU, since `.cargo/config.toml` sets no `target-cpu`. For local profiling, tune it to your machine with `RUSTFLAGS="-C target-cpu=native" cargo build --release`.

x86_64 release binaries are built with `cargo multivers --profile dist` as a single launcher that embeds three CPU-specific builds and picks the best match at startup:

- `x86-64`: SSE2 baseline, runs on any 64-bit x86 CPU (2003+)
- `x86-64-v2`: SSE4.2 + POPCNT (2008+); captures nearly all of the speed-up over the baseline
- `x86-64-v4`: AVX-512F/BW/CD/DQ/VL for Ice Lake / Sapphire Rapids / Granite Rapids / Zen 4+

The launcher is ~4 MB and adds ~0.2 s of startup for decompression + `memfd_create` + `exec`. v3 is skipped on purpose: on chelae's workload v2 and v3 are within measurement noise, and v4 adds what little AVX-512 gives (~1% on our benchmarks).

aarch64 release binaries (Apple Silicon, AWS Graviton, GCP Axion, Azure Cobalt) are a single build (`cargo build --profile dist`) for the generic ARMv8-A / NEON baseline. Neoverse-specific tuning gained only ~1-2% over generic in our benchmarks, with a near-zero penalty on other cores, so multivers isn't worth the complexity there.

## Code organization

The project follows conventions documented in [`CLAUDE.md`](CLAUDE.md) —
section ordering within command modules, impl-block collocation, and
callers-before-callees inside an impl. Please read `CLAUDE.md` before making
structural changes.

## Updating the usage docs

[`docs/usage.md`](docs/usage.md) has hand-curated **Options** tables that summarize every visible `chelae trim` and `chelae detect` option. When adding, removing, or renaming an option (or materially changing its meaning or default), update those tables so the docs stay in step with the CLI. The `--help` output remains the authoritative reference; the tables are short-form pointers to it.

## Changelog

`CHANGELOG.md` follows the [Keep a Changelog](https://keepachangelog.com/)
format and is promoted on release by `cargo-release`. As you make user-visible
changes, add a bullet under the `[Unreleased]` section heading using the
appropriate subsection (`Added`, `Changed`, `Deprecated`, `Removed`, `Fixed`,
or `Security`).

## Releasing a new version

### Pre-requisites

Install [`cargo-release`][cargo-release-link]:

```console
cargo install cargo-release
```

### Dry run

Verify the release plan without publishing:

```console
cargo release [major|minor|patch|release|rc] --no-publish
```

Dry-run is the default for `cargo-release`; add `--execute` once the output
looks right.

`release.toml` in the repo root configures the tag format (`vX.Y.Z`), the
CHANGELOG promotion (`[Unreleased]` → a versioned section), and disables
`crates.io` publishing from `cargo release` itself — publishing is done from
CI on tag.

See the [`cargo-release` reference documentation][cargo-release-docs-link]
for more.

### Semantic versioning

`chelae` follows [Semantic Versioning](https://semver.org/):

- **MAJOR** when you make incompatible API changes,
- **MINOR** when you add functionality in a backwards-compatible manner,
- **PATCH** when you make backwards-compatible bug fixes.

## Benchmark pipeline

The benchmark pipeline that produces the numbers in `README.md` and
`benchmark-pipeline/RESULTS.md` is documented separately in
[`benchmark-pipeline/README.md`](benchmark-pipeline/README.md). That's where
to look if you're modifying tool render scripts, adding a new tool to the
comparison, changing the sample matrix, or re-running the benchmark suite.

[cargo-release-link]:      https://github.com/crate-ci/cargo-release
[cargo-release-docs-link]: https://github.com/crate-ci/cargo-release/blob/master/docs/reference.md
