# Contributing to chelae

This document is for developers working on `chelae` itself. End-user
documentation lives in [README.md](README.md).

## Pre-push checks

Before pushing, run the full verification suite:

```console
bash ci/check.sh
```

This runs `cargo fmt`, `cargo clippy --all-features --all-targets -- -D warnings`,
and `cargo test`. CI also runs a `--locked` version via
`src/scripts/precommit.sh`.

For local development, `cargo build --release` uses `target-cpu=native` (see
`.cargo/config.toml`) for fastest local runs.

## Code organization

The project follows conventions documented in [`CLAUDE.md`](CLAUDE.md) —
section ordering within command modules, impl-block collocation, and
callers-before-callees inside an impl. Please read `CLAUDE.md` before making
structural changes.

## Updating the README

The README contains a hand-curated **Options** table that summarizes every
`chelae trim` flag. When adding, removing, or renaming an option (or
materially changing its meaning or default), update that table in `README.md`
so the public docs stay in step with the CLI. The full `chelae trim --help`
output remains the authoritative reference — the README table is a short-form
pointer to it.

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
