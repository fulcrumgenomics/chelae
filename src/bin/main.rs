//! `chelae` — a FASTQ trimming and filtering toolkit. This file is the CLI entry
//! point; it dispatches to the set of subcommands via the [`Command`] trait and
//! `enum_dispatch`. Installs `mimalloc` as the global allocator.

extern crate core;

pub mod commands;

use anyhow::Result;
use clap::Parser;
use commands::command::Command;
use commands::trim::Trim;
use enum_dispatch::enum_dispatch;
use env_logger::Env;

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Top-level `chelae` CLI parser. Holds a single [`Subcommand`] variant chosen by the
/// user's argv; every subcommand enum variant in turn wraps that subcommand's own
/// clap-derived options struct.
#[derive(Parser, Debug)]
struct Args {
    #[clap(subcommand)]
    subcommand: Subcommand,
}

/// Exhaustive list of `chelae` subcommands, wired to the [`Command`] trait via
/// `enum_dispatch` so `execute()` forwards to the chosen variant without a match arm
/// per-subcommand. Kept as a multi-subcommand enum even while `Trim` is the only
/// variant to leave room for future tools without a CLI shape change.
#[enum_dispatch(Command)]
#[derive(Parser, Debug)]
// `name`/`bin_name` are pinned explicitly: clap otherwise derives the displayed
// program name from argv[0]'s basename, and the cargo-multivers fexecve/memfd
// launcher passes `/proc/self/fd/N` as argv[0] under binfmt emulation (e.g. an
// amd64 biocontainer on Apple Silicon), which would print `Usage: 11 ...`.
// See fulcrumgenomics/riker#38.
#[command(name = "chelae", bin_name = "chelae", version)]
// Single enum instance per program; a large variant like Trim doesn't warrant a
// Box layer on the hot path.
#[allow(clippy::large_enum_variant)]
enum Subcommand {
    Trim(Trim),
}

/// Process entry point. Initializes env_logger with a default of `info` level,
/// parses argv, and dispatches to the selected subcommand's `execute()`. Errors
/// propagate out as `anyhow::Error` for the runtime's default `eprintln!` +
/// nonzero-exit handling.
fn main() -> Result<()> {
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
    let args: Args = Args::parse();
    args.subcommand.execute()
}
