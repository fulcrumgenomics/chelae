//! Integration tests for paired FASTQ inputs that are named pipes (FIFOs), spawning the
//! real `chelae` binary against a producer thread in this process.
#![cfg(unix)]

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};
use tempfile::TempDir;

/// How long chelae gets before a test calls it hung; a working run takes about a second.
const TIMEOUT: Duration = Duration::from_secs(60);

/// Pairs the producer writes. At ~250 bytes per record, R1's buffer fills and is
/// flushed whole after ~2,000 pairs, before anything has been written to R2.
const NUM_PAIRS: usize = 10_000;

/// Per-output buffer in the producer, matching `k2tools filter`.
const PRODUCER_BUFFER_SIZE: usize = 512 * 1024;

/// Every pair's insert: 67 bp, so each 100 bp read runs through into 33 bp of adapter.
const INSERT: &str = "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGG";
const R1_ADAPTER: &str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
const R2_ADAPTER: &str = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

/// Path to the compiled `chelae` binary under test, provided by Cargo for
/// integration tests (`tests/*.rs`).
fn chelae_bin() -> &'static str {
    env!("CARGO_BIN_EXE_chelae")
}

/// Creates a named pipe at `path`.
fn make_fifo(path: &Path) {
    let status = Command::new("mkfifo").arg(path).status().expect("failed to run mkfifo");
    assert!(status.success(), "mkfifo {path:?} failed");
}

/// Reverse complement of an ACGT sequence.
fn revcomp(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            _ => panic!("unexpected base {base}"),
        })
        .collect()
}

/// Starts a thread that writes `NUM_PAIRS` copies of the `INSERT` read pair to the named
/// pipes `r1` and `r2` the way `k2tools filter` does: it opens both pipes before writing to either, then
/// writes each pair's mates alternately through `PRODUCER_BUFFER_SIZE` buffers, so R1's
/// buffer is flushed whole before anything reaches R2. The thread is never joined: if
/// chelae doesn't open a pipe, it stays blocked in `open`.
fn spawn_producer(r1: PathBuf, r2: PathBuf, open_r2_first: bool) {
    std::thread::spawn(move || -> std::io::Result<()> {
        let open = |path: &Path| {
            File::create(path).map(|f| BufWriter::with_capacity(PRODUCER_BUFFER_SIZE, f))
        };
        let (mut w1, mut w2) = if open_r2_first {
            let w2 = open(&r2)?;
            (open(&r1)?, w2)
        } else {
            let w1 = open(&r1)?;
            (w1, open(&r2)?)
        };
        let r1_seq = format!("{INSERT}{R1_ADAPTER}");
        let r2_seq = format!("{}{R2_ADAPTER}", revcomp(INSERT));
        let quals = "I".repeat(r1_seq.len());
        for i in 0..NUM_PAIRS {
            write!(w1, "@pair{i}/1\n{r1_seq}\n+\n{quals}\n")?;
            write!(w2, "@pair{i}/2\n{r2_seq}\n+\n{quals}\n")?;
        }
        w1.flush()?;
        w2.flush()
    });
}

/// Runs `chelae` with `args`, killing it and failing the test if it hasn't exited
/// within `TIMEOUT`. Asserts a zero exit status.
fn run_chelae_with_timeout(args: &[&str], tmp: &Path) {
    let stderr_path = tmp.join("stderr.txt");
    let mut child = Command::new(chelae_bin())
        .args(args)
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(File::create(&stderr_path).unwrap())
        .spawn()
        .expect("failed to spawn chelae");
    let deadline = Instant::now() + TIMEOUT;
    let status = loop {
        if let Some(status) = child.try_wait().unwrap() {
            break status;
        }
        if Instant::now() > deadline {
            child.kill().unwrap();
            child.wait().unwrap();
            panic!(
                "chelae {args:?} still running after {TIMEOUT:?}:\n{}",
                std::fs::read_to_string(&stderr_path).unwrap()
            );
        }
        std::thread::sleep(Duration::from_millis(50));
    };
    assert!(
        status.success(),
        "chelae {args:?} failed ({status}):\n{}",
        std::fs::read_to_string(&stderr_path).unwrap()
    );
}

/// Number of FASTQ records in the plain-text file at `path`.
fn count_records(path: &Path) -> usize {
    std::fs::read_to_string(path).unwrap().lines().count() / 4
}

/// Runs `chelae trim` over two named pipes fed by [`spawn_producer`] and asserts that
/// every pair comes out.
fn assert_trim_reads_every_pair_from_named_pipes(open_r2_first: bool) {
    let tmp = TempDir::new().unwrap();
    let (r1, r2) = (tmp.path().join("r1.fq"), tmp.path().join("r2.fq"));
    let (out1, out2) = (tmp.path().join("out1.fq"), tmp.path().join("out2.fq"));
    make_fifo(&r1);
    make_fifo(&r2);
    spawn_producer(r1.clone(), r2.clone(), open_r2_first);

    run_chelae_with_timeout(
        &[
            "trim",
            "-i",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "-o",
            out1.to_str().unwrap(),
            out2.to_str().unwrap(),
        ],
        tmp.path(),
    );

    assert_eq!(count_records(&out1), NUM_PAIRS);
    assert_eq!(count_records(&out2), NUM_PAIRS);
}

#[test]
fn trim_reads_named_pipes_whose_producer_opens_r1_then_r2_before_writing() {
    assert_trim_reads_every_pair_from_named_pipes(false);
}

#[test]
fn trim_reads_named_pipes_whose_producer_opens_r2_then_r1_before_writing() {
    assert_trim_reads_every_pair_from_named_pipes(true);
}

#[test]
fn detect_reads_named_pipes_whose_producer_opens_both_before_writing() {
    let tmp = TempDir::new().unwrap();
    let (r1, r2) = (tmp.path().join("r1.fq"), tmp.path().join("r2.fq"));
    make_fifo(&r1);
    make_fifo(&r2);
    spawn_producer(r1.clone(), r2.clone(), false);

    run_chelae_with_timeout(
        &["detect", "-i", r1.to_str().unwrap(), r2.to_str().unwrap()],
        tmp.path(),
    );
}
