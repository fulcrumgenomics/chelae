//! Integration tests for stdin/stdout FASTQ I/O, spawning the real `chelae` binary.
//!
//! These paths can't be exercised through the in-process `execute()` tests in
//! `commands/trim.rs` / `commands/detect.rs` because stdin/stdout are process-level
//! resources — `std::io::stdin()` inside the library always refers to the *test
//! harness's* stdin, not a value we can substitute per-test. Spawning the compiled
//! binary with piped stdin/stdout is the only way to drive that code path.

use flate2::Compression;
use flate2::write::GzEncoder;
use std::io::{Read, Write};
use std::process::{Command, Output, Stdio};
use tempfile::TempDir;

/// Path to the compiled `chelae` binary under test, provided by Cargo for
/// integration tests (`tests/*.rs`).
fn chelae_bin() -> &'static str {
    env!("CARGO_BIN_EXE_chelae")
}

/// Runs `chelae` with `args`, piping `stdin_bytes` to its stdin and capturing
/// stdout/stderr. Asserts a zero exit status (with stderr in the panic message on
/// failure) and returns the raw stdout bytes.
fn run_chelae_ok(args: &[&str], stdin_bytes: &[u8]) -> Vec<u8> {
    let output = spawn_chelae(args, stdin_bytes);
    assert!(
        output.status.success(),
        "chelae {args:?} failed (status {:?}): {}",
        output.status.code(),
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}

/// Runs `chelae` with `args` and `stdin_bytes`, returning the full captured
/// `Output` (including a non-zero exit status) without asserting success.
///
/// Writes `stdin_bytes` on a separate scoped thread, running concurrently with
/// `wait_with_output()`'s stdout-draining read loop, rather than writing sequentially
/// before reading: a payload larger than the OS pipe buffer would otherwise deadlock
/// (child blocks writing stdout because nothing is reading it yet; we're blocked
/// writing stdin because the child isn't reading it yet, busy blocking on stdout).
/// `thread::scope` joins the writer thread only after `wait_with_output()` returns.
fn spawn_chelae(args: &[&str], stdin_bytes: &[u8]) -> Output {
    let mut child = Command::new(chelae_bin())
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to spawn chelae");
    // Dropping the `ChildStdin` when the writer thread finishes closes the pipe,
    // signaling EOF to the child — required for it to stop reading stdin.
    let mut child_stdin = child.stdin.take().unwrap();
    std::thread::scope(|s| {
        s.spawn(move || {
            // `move` so the write thread owns (and, on return, drops/closes) the
            // pipe — closing it is what signals EOF to the child. Without `move`
            // the outer `child_stdin` binding would only be dropped when this
            // function returns, which can't happen until the child sees EOF and
            // exits: a self-inflicted deadlock.
            //
            // A chelae process that exits early (e.g. its downstream stdout reader
            // went away — see `stdout_closed_early_exits_promptly`) closes its
            // stdin too, which can fail this write with BrokenPipe; that's an
            // expected race for this helper, not a bug, so it's swallowed rather
            // than unwrapped.
            let _ = child_stdin.write_all(stdin_bytes);
        });
        child.wait_with_output().expect("failed to wait on chelae")
    })
}

/// One FASTQ record's 4 lines as a string (trailing newline included). Qualities are
/// all `I` (Q40 at Phred+33).
fn fq_record(name: &str, seq: &str) -> String {
    format!("@{name}\n{seq}\n+\n{}\n", "I".repeat(seq.len()))
}

/// Plain-text FASTQ for `n` single-end reads named `read{i}`, all with sequence `seq`.
fn se_fastq_text(n: usize, seq: &str) -> String {
    (0..n).map(|i| fq_record(&format!("read{i}"), seq)).collect()
}

/// Interleaved PE FASTQ text (R1, R2, R1, R2, ...) for `n` pairs, mate names
/// distinguished by a `/1` `/2` suffix.
fn interleaved_fastq_text(n: usize, r1_seq: &str, r2_seq: &str) -> String {
    let mut s = String::new();
    for i in 0..n {
        s += &fq_record(&format!("pair{i}/1"), r1_seq);
        s += &fq_record(&format!("pair{i}/2"), r2_seq);
    }
    s
}

/// Gzip-compresses `data` in memory (default compression level).
fn gzip(data: &[u8]) -> Vec<u8> {
    let mut enc = GzEncoder::new(Vec::new(), Compression::default());
    enc.write_all(data).unwrap();
    enc.finish().unwrap()
}

#[test]
fn stdin_plain_to_file_out() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("out.fq");
    let input = se_fastq_text(5, "ACGTACGTACGTACGTACGT");

    run_chelae_ok(
        &["trim", "-i", "-", "-o", out.to_str().unwrap(), "--output-compression", "none"],
        input.as_bytes(),
    );

    let written = std::fs::read_to_string(&out).unwrap();
    assert_eq!(written.matches("@read").count(), 5);
}

#[test]
fn stdin_gzip_to_file_out() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("out.fq");
    let input = gzip(se_fastq_text(5, "ACGTACGTACGTACGTACGT").as_bytes());

    run_chelae_ok(
        &["trim", "-i", "-", "-o", out.to_str().unwrap(), "--output-compression", "none"],
        &input,
    );

    let written = std::fs::read_to_string(&out).unwrap();
    assert_eq!(written.matches("@read").count(), 5);
}

#[test]
fn no_input_flag_defaults_to_stdin() {
    // Omit `-i` entirely (rather than passing `-i -`) to exercise the default-to-`-`
    // path specifically.
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("out.fq");
    let input = se_fastq_text(3, "ACGTACGTACGTACGTACGT");

    run_chelae_ok(
        &["trim", "-o", out.to_str().unwrap(), "--output-compression", "none"],
        input.as_bytes(),
    );

    let written = std::fs::read_to_string(&out).unwrap();
    assert_eq!(written.matches("@read").count(), 3);
}

#[test]
fn file_in_to_stdout_is_plain_fastq() {
    let tmp = TempDir::new().unwrap();
    let in_path = tmp.path().join("in.fq");
    std::fs::write(&in_path, se_fastq_text(4, "ACGTACGTACGTACGTACGT")).unwrap();

    // No `-o` at all: defaults to stdout.
    let stdout = run_chelae_ok(&["trim", "-i", in_path.to_str().unwrap()], &[]);

    assert!(!stdout.starts_with(&[0x1f, 0x8b]), "expected plain text on stdout by default");
    let text = String::from_utf8(stdout).unwrap();
    assert_eq!(text.matches("@read").count(), 4);
}

#[test]
fn file_in_to_stdout_is_bgzf_with_override() {
    let tmp = TempDir::new().unwrap();
    let in_path = tmp.path().join("in.fq");
    std::fs::write(&in_path, se_fastq_text(4, "ACGTACGTACGTACGTACGT")).unwrap();

    let stdout = run_chelae_ok(
        &["trim", "-i", in_path.to_str().unwrap(), "-o", "-", "--output-compression", "bgzf"],
        &[],
    );

    assert!(stdout.starts_with(&[0x1f, 0x8b]), "expected gzip/BGZF magic bytes on stdout");
}

#[test]
fn interleaved_gz_stdin_to_interleaved_plain_stdout() {
    let input =
        gzip(interleaved_fastq_text(3, "AAAACCCCTTAAAACCCCTT", "GGGGTTTTAAGGGGTTTTAA").as_bytes());

    let stdout =
        run_chelae_ok(&["trim", "-i", "-", "-o", "-", "--output-compression", "none"], &input);

    assert!(!stdout.starts_with(&[0x1f, 0x8b]), "expected plain text on stdout");
    let text = String::from_utf8(stdout).unwrap();
    let heads: Vec<&str> = text.lines().filter(|l| l.starts_with('@')).collect();
    assert_eq!(heads, vec!["@pair0/1", "@pair0/2", "@pair1/1", "@pair1/2", "@pair2/1", "@pair2/2"]);
}

#[test]
fn detect_output_fasta_dash_writes_stdout() {
    // TruSeq R1 adapter readthrough tail on every read, enough reads to clear the
    // default `--min-detections-for-report` floor.
    let template = "ACGTGACCTGATTGCAACGATCGTAGCTAGCATCGATCGATTAGCGATCGA";
    let adapter_tail = "AGATCGGAAGAGCACACGTCTGA";
    let seq = format!("{template}{adapter_tail}");
    let input = se_fastq_text(200, &seq);

    let stdout = run_chelae_ok(&["detect", "-i", "-", "-o", "-"], input.as_bytes());

    let text = String::from_utf8(stdout).unwrap();
    assert!(text.trim_start().starts_with('>'), "expected FASTA on stdout, got:\n{text}");
    assert!(text.contains("truseq"), "expected the truseq kit name in the FASTA, got:\n{text}");
}

#[test]
fn detect_stdin_input() {
    let tmp = TempDir::new().unwrap();
    let out = tmp.path().join("adapters.fa");
    let template = "ACGTGACCTGATTGCAACGATCGTAGCTAGCATCGATCGATTAGCGATCGA";
    let adapter_tail = "AGATCGGAAGAGCACACGTCTGA";
    let seq = format!("{template}{adapter_tail}");
    let input = se_fastq_text(200, &seq);

    run_chelae_ok(&["detect", "-i", "-", "-o", out.to_str().unwrap()], input.as_bytes());

    let fasta = std::fs::read_to_string(&out).unwrap();
    assert!(fasta.contains("truseq"), "expected the truseq kit name in the FASTA, got:\n{fasta}");
}

/// Spawns `chelae` with `args`, feeds `stdin_bytes` from a writer thread, reads (and
/// then closes) the first 4 KB of the child's stdout — mimicking `| head -c 4096` —
/// and polls for exit with a 10-second bound so a regression back to the old
/// (error-out or hang) behavior fails the test instead of hanging the suite.
/// Returns the child's exit status. The stdin writer swallows a BrokenPipe: once the
/// child exits (having stopped reading stdin), the remaining `write_all` legitimately
/// fails (see `spawn_chelae`'s doc comment for the same reasoning).
fn run_until_stdout_closed(args: &[&str], stdin_bytes: Vec<u8>) -> std::process::ExitStatus {
    let mut child = Command::new(chelae_bin())
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .expect("failed to spawn chelae");

    let mut child_stdin = child.stdin.take().unwrap();
    let writer = std::thread::spawn(move || {
        let _ = child_stdin.write_all(&stdin_bytes);
    });

    let mut child_stdout = child.stdout.take().unwrap();
    let mut buf = [0u8; 4096];
    std::io::Read::read(&mut child_stdout, &mut buf).expect("failed to read initial stdout bytes");
    drop(child_stdout);

    let start = std::time::Instant::now();
    let status = loop {
        if let Some(status) = child.try_wait().expect("failed to poll chelae") {
            break status;
        }
        assert!(
            start.elapsed() < std::time::Duration::from_secs(10),
            "chelae did not exit promptly after its stdout was closed"
        );
        std::thread::sleep(std::time::Duration::from_millis(20));
    };
    writer.join().unwrap();
    status
}

#[test]
fn stdout_closed_early_exits_promptly() {
    // 250k single-end records serializes to several MB of plain-text output —
    // comfortably larger than any OS pipe buffer (tens of KB) — so the writer_loop
    // thread is still blocked writing (nobody's reading) when the helper closes its
    // end of the stdout pipe, guaranteeing the BrokenPipe path actually runs rather
    // than the child having already finished and exited on its own.
    let input = se_fastq_text(250_000, "ACGTACGTACGTACGTACGT");

    let status = run_until_stdout_closed(
        &["trim", "-i", "-", "-o", "-", "--output-compression", "none"],
        input.into_bytes(),
    );

    assert!(status.success(), "expected exit 0 after stdout closed early, got {status:?}");
}

#[test]
fn stdout_closed_early_leaves_split_file_output_valid() {
    // Split PE output: mate 1 to stdout, mate 2 to a file. Both outputs share one
    // batch stream from the reader loop, so once stdout closes and the reader loop
    // stops submitting further batches, the file writer only sees a partial record
    // set — but it must still finish cleanly (BGZF EOF block, normal flush) rather
    // than surfacing a spurious "channel closed" error just because its sibling
    // writer hit BrokenPipe.
    let tmp = TempDir::new().unwrap();
    let out2 = tmp.path().join("out2.fq.gz");
    let input = interleaved_fastq_text(250_000, "ACGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCA");

    let status = run_until_stdout_closed(
        &["trim", "-i", "-", "-o", "-", "-o", out2.to_str().unwrap()],
        input.into_bytes(),
    );

    assert!(status.success(), "expected exit 0, got {status:?}");

    let bytes = std::fs::read(&out2).unwrap();
    assert!(!bytes.is_empty(), "expected non-empty partial output on the file sink");
    let mut decoded = Vec::new();
    flate2::bufread::MultiGzDecoder::new(bytes.as_slice())
        .read_to_end(&mut decoded)
        .expect("out2.fq.gz should be a valid, non-truncated BGZF stream");
    assert!(decoded.starts_with(b"@"), "expected valid FASTQ content, got:\n{decoded:?}");
}
