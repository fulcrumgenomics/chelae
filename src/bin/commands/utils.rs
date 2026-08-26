//! Small utilities shared across subcommands.

use anyhow::{Result, anyhow};
use flate2::bufread::MultiGzDecoder;
use log::info;
use seq_io::fastq::OwnedRecord;
use seq_io::fastq::Reader as FastqReader;
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, Read};
use std::path::{Path, PathBuf};

/// BufReader / BufWriter capacity used by the FASTQ I/O paths in every
/// subcommand. 512 KiB chosen via a 32k→4MiB sweep on Graviton4 (Neoverse-V2) and
/// x86 Granite Rapids (c8i): wall time and cycle count are flat across 256k–2MiB
/// on both architectures, with instruction count showing a shallow U-shape that
/// bottoms out at 512k–1024k. 512k is at the floor everywhere tested and halves
/// resident memory per reader/writer vs 1 MiB.
pub(crate) const BUFFER_SIZE: usize = 512 * 1024;

/// The two leading bytes of every gzip (and BGZF, since BGZF is gzip-framed)
/// stream. Used to sniff compression by content rather than file extension, which
/// is required for stdin (no extension) and also fixes misnamed files.
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// How a paired FASTQ's read names encode mate identity. Selected once from the first
/// pair (see [`PairingRule::select`]) and applied to every subsequent pair via
/// [`PairingRule::check_pair`], so steady-state checking is a stem compare + marker
/// check rather than re-deriving the convention per pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PairingRule {
    /// First whitespace-delimited tokens are identical. If both headers carry a comment
    /// whose first character is `1`/`2` (Casava 1.8+), mate 1 must be `1` and mate 2 `2`;
    /// with no such comments, bare token equality suffices.
    CasavaOrBare,
    /// Token ends in `/1` (mate 1) / `/2` (mate 2); stems before the suffix are equal.
    SlashDigit,
    /// Token ends in `<sep>1` / `<sep>2` where `sep` is `.` or `_` (SRA `fastq-dump -I`
    /// and common pipeline conventions); stems before the suffix are equal.
    SepDigit(u8),
}

impl PairingRule {
    /// Picks the pairing convention that explains `head1`/`head2` as a mate-1/mate-2
    /// pair (in that order — a reversed pair never selects a rule), trying
    /// `SlashDigit`, then `SepDigit('.')`, then `SepDigit('_')`, then `CasavaOrBare`.
    /// `None` if no rule explains the pair.
    pub(crate) fn select(head1: &[u8], head2: &[u8]) -> Option<PairingRule> {
        if matches_slash_digit(head1, head2) {
            Some(PairingRule::SlashDigit)
        } else if matches_sep_digit(head1, head2, b'.') {
            Some(PairingRule::SepDigit(b'.'))
        } else if matches_sep_digit(head1, head2, b'_') {
            Some(PairingRule::SepDigit(b'_'))
        } else if matches_casava_or_bare(head1, head2) {
            Some(PairingRule::CasavaOrBare)
        } else {
            None
        }
    }

    /// Verifies that `head1`/`head2` are a valid mate-1/mate-2 pair (in that order)
    /// under this already-selected rule. A reversed pair (mate 2 first) always fails.
    pub(crate) fn check_pair(&self, head1: &[u8], head2: &[u8]) -> bool {
        match self {
            PairingRule::CasavaOrBare => matches_casava_or_bare(head1, head2),
            PairingRule::SlashDigit => matches_slash_digit(head1, head2),
            PairingRule::SepDigit(sep) => matches_sep_digit(head1, head2, *sep),
        }
    }
}

impl std::fmt::Display for PairingRule {
    /// Human-readable rule name for the sniff-decision log line (e.g. `mate suffix '.1'/'.2'`),
    /// rendering the separator as a character rather than `Debug`'s raw byte value.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PairingRule::CasavaOrBare => write!(f, "identical names (Casava comment or bare)"),
            PairingRule::SlashDigit => write!(f, "mate suffix '/1'/'/2'"),
            PairingRule::SepDigit(sep) => {
                let s = *sep as char;
                write!(f, "mate suffix '{s}1'/'{s}2'")
            }
        }
    }
}

/// Splits a FASTQ header into its first whitespace-delimited token and the remaining
/// comment (if any), matching the Casava 1.8+ `<id> <comment>` convention (and
/// degrading gracefully for headers with no comment).
fn header_token_and_comment(head: &[u8]) -> (&[u8], Option<&[u8]>) {
    match head.iter().position(|&b| b == b' ' || b == b'\t') {
        Some(i) => (&head[..i], Some(&head[i + 1..])),
        None => (head, None),
    }
}

/// If `token` ends in `<sep><digit>`, returns the stem before that suffix.
fn strip_suffix_digit(token: &[u8], sep: u8, digit: u8) -> Option<&[u8]> {
    if token.len() >= 2 && token[token.len() - 2] == sep && token[token.len() - 1] == digit {
        Some(&token[..token.len() - 2])
    } else {
        None
    }
}

/// `PairingRule::SlashDigit`: token 1 ends `/1`, token 2 ends `/2`, stems equal.
fn matches_slash_digit(head1: &[u8], head2: &[u8]) -> bool {
    let (t1, _) = header_token_and_comment(head1);
    let (t2, _) = header_token_and_comment(head2);
    matches!(
        (strip_suffix_digit(t1, b'/', b'1'), strip_suffix_digit(t2, b'/', b'2')),
        (Some(s1), Some(s2)) if s1 == s2
    )
}

/// `PairingRule::SepDigit(sep)`: token 1 ends `<sep>1`, token 2 ends `<sep>2`, stems equal.
fn matches_sep_digit(head1: &[u8], head2: &[u8], sep: u8) -> bool {
    let (t1, _) = header_token_and_comment(head1);
    let (t2, _) = header_token_and_comment(head2);
    matches!(
        (strip_suffix_digit(t1, sep, b'1'), strip_suffix_digit(t2, sep, b'2')),
        (Some(s1), Some(s2)) if s1 == s2
    )
}

/// `PairingRule::CasavaOrBare`: bare tokens equal; if both headers carry a comment
/// whose first byte is `1`/`2`, the mate-1 comment must start `1` and mate-2 `2`.
fn matches_casava_or_bare(head1: &[u8], head2: &[u8]) -> bool {
    let (t1, c1) = header_token_and_comment(head1);
    let (t2, c2) = header_token_and_comment(head2);
    if t1 != t2 {
        return false;
    }
    match (c1.and_then(|c| c.first()), c2.and_then(|c| c.first())) {
        (Some(b'1' | b'2'), Some(b'1' | b'2')) => c1.unwrap()[0] == b'1' && c2.unwrap()[0] == b'2',
        _ => true,
    }
}

/// Confirms `head1`/`head2` pair under `*rule`, selecting the rule from the first pair
/// seen if one hasn't been chosen yet. Shared by `trim`'s split-file zipper
/// ([`crate::commands::trim`]'s `pull_per_file_slot`) and `detect`'s
/// `PairSource::Split` so the lazy-rule-selection error text is identical in both.
/// `where_` is appended to the mismatch error verbatim (e.g. `" at record 5"` or `""`).
pub(crate) fn check_or_select_split_pair(
    rule: &mut Option<PairingRule>,
    head1: &[u8],
    head2: &[u8],
    where_: &str,
) -> Result<()> {
    let mismatch = || {
        anyhow!(
            "R1/R2 read names do not correspond{where_}: {:?} / {:?}",
            String::from_utf8_lossy(head1),
            String::from_utf8_lossy(head2),
        )
    };
    match rule {
        Some(r) => {
            if !r.check_pair(head1, head2) {
                return Err(mismatch());
            }
        }
        None => *rule = Some(PairingRule::select(head1, head2).ok_or_else(mismatch)?),
    }
    Ok(())
}

/// Pulls one pair from a single interleaved-input iterator, checking it against the
/// already-selected `rule`. `Ok(None)` signals a clean EOF between pairs; errors report
/// both the pair index and the underlying file-record indices. Shared by `trim`'s
/// zipper and `detect`'s `PairSource::Interleaved` so the odd-count check, pair check,
/// and error strings live in one place.
pub(crate) fn pull_pair_interleaved<I>(
    iter: &mut I,
    rule: PairingRule,
    pair_idx: u64,
) -> Result<Option<(OwnedRecord, OwnedRecord)>>
where
    I: Iterator<Item = Result<OwnedRecord>>,
{
    let r1 = match iter.next() {
        None => return Ok(None),
        Some(Ok(r)) => r,
        Some(Err(e)) => return Err(e),
    };
    let r2 = match iter.next() {
        None => anyhow::bail!(
            "interleaved input ended mid-pair (odd record count; file truncated?) at pair \
             {pair_idx} (file record {})",
            pair_idx * 2 - 1,
        ),
        Some(Ok(r)) => r,
        Some(Err(e)) => return Err(e),
    };
    anyhow::ensure!(
        rule.check_pair(&r1.head, &r2.head),
        "interleaved input out of sync at pair {pair_idx} (file records {}/{}): {:?} / {:?} are \
         not a pair",
        pair_idx * 2 - 1,
        pair_idx * 2,
        String::from_utf8_lossy(&r1.head),
        String::from_utf8_lossy(&r2.head),
    );
    Ok(Some((r1, r2)))
}

/// Wraps a `FastqReader` as an `Iterator<Item = Result<OwnedRecord>>`. Moving
/// ownership of the reader into the iterator makes it easy to hand off to a
/// background thread (e.g. via fgoxide's `read_ahead`), or to peek a few records
/// and replay them via `Iterator::chain` (see [`sniff_single_input`]).
pub(crate) struct OwnedRecordIter {
    pub(crate) reader: FastqReader<Box<dyn BufRead + Send>>,
}

impl Iterator for OwnedRecordIter {
    type Item = Result<OwnedRecord>;

    /// Advances the underlying [`FastqReader`] and materializes each `RefRecord` as an
    /// `OwnedRecord` so it can cross thread boundaries. Parse errors are wrapped in
    /// `anyhow::Error` with context.
    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.next()? {
            Ok(refrec) => Some(Ok(refrec.to_owned_record())),
            Err(e) => Some(Err(anyhow!("FASTQ read error: {e}"))),
        }
    }
}

/// Opens one FASTQ input path as a boxed, buffered reader. `-` means stdin.
/// Compression is detected by reading up to the first two bytes and checking for the
/// gzip magic number rather than by file extension. Reads in a small loop (rather than
/// a single `fill_buf`) because a single `fill_buf` call on a pipe can legally return
/// just 1 byte, which would otherwise misdetect a gzip stream as plain text.
fn open_one_fastq_input(path: &Path) -> Result<Box<dyn BufRead + Send>> {
    let mut inner: Box<dyn BufRead + Send> = if path.as_os_str() == "-" {
        Box::new(BufReader::with_capacity(BUFFER_SIZE, std::io::stdin()))
    } else {
        let file = File::open(path).map_err(|e| anyhow!("Failed to open input {path:?}: {e}"))?;
        Box::new(BufReader::with_capacity(BUFFER_SIZE, file))
    };

    let mut probe = [0u8; 2];
    let mut probed = 0usize;
    while probed < probe.len() {
        let n = inner
            .read(&mut probe[probed..])
            .map_err(|e| anyhow!("Failed to read input {path:?}: {e}"))?;
        if n == 0 {
            break; // EOF before 2 bytes accumulated; too short to be gzip.
        }
        probed += n;
    }
    let prefix = Cursor::new(probe[..probed].to_vec());
    let chained: Box<dyn BufRead + Send> =
        Box::new(BufReader::with_capacity(BUFFER_SIZE, prefix.chain(inner)));

    if probed == probe.len() && probe == GZIP_MAGIC {
        Ok(Box::new(BufReader::with_capacity(BUFFER_SIZE, MultiGzDecoder::new(chained))))
    } else {
        Ok(chained)
    }
}

/// Opens every input path as a [`FastqReader`] with `BUFFER_SIZE` capacity. `-`
/// means stdin; gzip/BGZF is detected by content, not extension (see
/// [`open_one_fastq_input`]).
pub(crate) fn open_fastq_inputs(
    paths: &[PathBuf],
) -> Result<Vec<FastqReader<Box<dyn BufRead + Send>>>> {
    paths
        .iter()
        .map(|p| open_one_fastq_input(p).map(|r| FastqReader::with_capacity(r, BUFFER_SIZE)))
        .collect()
}

/// Result of [`sniff_single_input`]: whether the lone input was detected as
/// interleaved paired-end, the pairing rule selected for ongoing enforcement (`None`
/// for single-end, or when fewer than 2 records were seen), whether the input produced
/// no records at all (drives the layout-agnostic empty-input rule — see the callers in
/// `trim`/`detect`), and an iterator over the entire stream (peeked records replayed via
/// `Iterator::chain` so nothing is lost to the peek).
pub(crate) struct Sniffed<I> {
    pub(crate) interleaved: bool,
    pub(crate) pairing_rule: Option<PairingRule>,
    pub(crate) is_empty: bool,
    pub(crate) records: I,
}

/// Sniffs whether a single FASTQ input is single-end or interleaved paired-end by
/// peeking up to its first 4 records and returns the sniff result alongside a
/// [`Sniffed`] wrapping an iterator over the *entire* stream.
///
/// Rule: 0 or 1 records sniff as single-end. Otherwise a [`PairingRule`] is selected
/// from records 1–2 ([`PairingRule::select`]); no rule matching sniffs as single-end.
/// If a rule matches and records 3–4 both exist, the input is interleaved only if
/// records 3–4 *also* pair under that same rule — this rescues a standard SE SRA file
/// (`@SRR.1`, `@SRR.2`, `@SRR.3`, …) whose first two records would otherwise
/// misdetect as an SRA-interleaved (`.1`/`.2`-suffixed) pair, since a true SE file's
/// records 3–4 (`.3`/`.4`) fail the `SepDigit` orientation check. If a rule matches
/// and the file has fewer than 4 records (2 or 3), the input is interleaved on the
/// records 1–2 match alone — an inherent, accepted ambiguity for very short SE SRA
/// files (2 records only). Loudly `info!`s the decision and the selected rule.
pub(crate) fn sniff_single_input(
    reader: FastqReader<Box<dyn BufRead + Send>>,
) -> Result<Sniffed<impl Iterator<Item = Result<OwnedRecord>> + Send + 'static>> {
    let mut iter = OwnedRecordIter { reader };
    let mut peeked: Vec<OwnedRecord> = Vec::with_capacity(4);
    for _ in 0..4 {
        match iter.next() {
            Some(Ok(r)) => peeked.push(r),
            Some(Err(e)) => return Err(e),
            None => break,
        }
    }

    let mut interleaved = false;
    let mut pairing_rule: Option<PairingRule> = None;
    if peeked.len() >= 2
        && let Some(rule) = PairingRule::select(&peeked[0].head, &peeked[1].head)
    {
        let confirmed = peeked.len() < 4 || rule.check_pair(&peeked[2].head, &peeked[3].head);
        if confirmed {
            interleaved = true;
            pairing_rule = Some(rule);
        }
    }

    if interleaved {
        info!(
            "Single input sniffed as interleaved paired-end (rule: {}; records {:?} / {:?} pair)",
            pairing_rule.expect("set alongside interleaved"),
            String::from_utf8_lossy(&peeked[0].head),
            String::from_utf8_lossy(&peeked[1].head),
        );
    } else {
        info!("Single input sniffed as single-end");
    }

    let is_empty = peeked.is_empty();
    let preface: Vec<Result<OwnedRecord>> = peeked.into_iter().map(Ok).collect();
    Ok(Sniffed { interleaved, pairing_rule, is_empty, records: preface.into_iter().chain(iter) })
}

/// Defaults an empty `--inputs`/`--outputs` list to `["-"]` (stdin/stdout), now
/// that both flags are optional on `trim` and `-i` is optional on `detect`.
pub(crate) fn default_dash(raw: &[PathBuf]) -> Vec<PathBuf> {
    if raw.is_empty() { vec![PathBuf::from("-")] } else { raw.to_vec() }
}

/// Applies [`default_dash`] to an input list, then guards against reading FASTQ
/// from an interactive terminal — whether `-` was defaulted or given explicitly.
/// There is no equivalent guard for outputs: writing FASTQ to a terminal (e.g.
/// eyeballing a few reads with `chelae trim | head`) is always allowed.
pub(crate) fn resolve_inputs(raw: &[PathBuf], stdin_is_tty: bool) -> Result<Vec<PathBuf>> {
    let resolved = default_dash(raw);
    if stdin_is_tty && resolved.iter().any(|p| p.as_os_str() == "-") {
        return Err(anyhow!("stdin is a terminal; pass --inputs or pipe data in"));
    }
    Ok(resolved)
}

/// Appends an error to `errors` if `-` (stdin/stdout) appears more than once in
/// `paths` — chelae has exactly one stdin and one stdout, so at most one input
/// (or output) may claim it.
pub(crate) fn check_dash_at_most_once(paths: &[PathBuf], label: &str, errors: &mut Vec<String>) {
    if paths.iter().filter(|p| p.as_os_str() == "-").count() > 1 {
        errors.push(format!("{label} may specify '-' (stdin/stdout) at most once."));
    }
}

/// Aggregates a list of user-facing validation-error strings into one `Result`, formatted
/// as a bulleted list. Empty `errors` returns `Ok(())`. Shared by `trim` and `detect`'s
/// `validate()` (and `trim`'s post-detection validation pass) so every subcommand reports
/// multi-error validation failures in the same shape.
pub(crate) fn aggregate_errors(errors: Vec<String>) -> Result<()> {
    if errors.is_empty() {
        return Ok(());
    }
    use std::fmt::Write;
    let detail = errors.iter().fold(String::new(), |mut s, e| {
        let _ = writeln!(s, "    - {e}");
        s
    });
    Err(anyhow!("Input validation failed:\n{detail}"))
}

/// Formats a `u64` with comma thousands-separators (e.g. `1,234,567`).
pub(crate) fn fmt_count(n: u64) -> String {
    let s = n.to_string();
    let mut result = String::with_capacity(s.len() + s.len() / 3);
    for (i, ch) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            result.push(',');
        }
        result.push(ch);
    }
    result
}

#[cfg(test)]
/// Builds `n` FASTQ records' worth of bytes (4 lines each, `I`-quality). Shared by the
/// `open_one_fastq_input` / `sniff_single_input` unit tests below.
pub(crate) fn test_fastq_bytes(records: &[(&str, &str)]) -> Vec<u8> {
    let mut out = Vec::new();
    for (name, seq) in records {
        out.extend_from_slice(format!("@{name}\n{seq}\n+\n{}\n", "I".repeat(seq.len())).as_bytes());
    }
    out
}

#[cfg(test)]
/// Gzip-compresses `data` in memory (default compression level). Shared by the
/// `open_one_fastq_input` unit tests below.
pub(crate) fn test_gzip(data: &[u8]) -> Vec<u8> {
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write as _;
    let mut enc = GzEncoder::new(Vec::new(), Compression::default());
    enc.write_all(data).unwrap();
    enc.finish().unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small() {
        assert_eq!(fmt_count(0), "0");
        assert_eq!(fmt_count(999), "999");
    }

    #[test]
    fn with_commas() {
        assert_eq!(fmt_count(1_000), "1,000");
        assert_eq!(fmt_count(1_234_567), "1,234,567");
        assert_eq!(fmt_count(1_000_000_000), "1,000,000,000");
    }

    // ---- PairingRule::select ----

    #[test]
    fn select_casava_style() {
        assert_eq!(
            PairingRule::select(b"read1 1:N:0:ATCG", b"read1 2:N:0:ATCG"),
            Some(PairingRule::CasavaOrBare)
        );
    }

    #[test]
    fn select_bare_identical_names() {
        assert_eq!(PairingRule::select(b"SRR1.1", b"SRR1.1"), Some(PairingRule::CasavaOrBare));
    }

    #[test]
    fn select_slash_suffix_style() {
        assert_eq!(PairingRule::select(b"read1/1", b"read1/2"), Some(PairingRule::SlashDigit));
    }

    #[test]
    fn select_dot_suffix_style() {
        assert_eq!(
            PairingRule::select(b"SRR000001.1.1", b"SRR000001.1.2"),
            Some(PairingRule::SepDigit(b'.'))
        );
    }

    #[test]
    fn select_underscore_suffix_style() {
        assert_eq!(PairingRule::select(b"read1_1", b"read1_2"), Some(PairingRule::SepDigit(b'_')));
    }

    #[test]
    fn select_none_on_mismatched_stems() {
        assert_eq!(PairingRule::select(b"read1/1", b"read2/2"), None);
    }

    #[test]
    fn select_none_on_reversed_slash_pair() {
        assert_eq!(PairingRule::select(b"read1/2", b"read1/1"), None);
    }

    #[test]
    fn select_none_on_reversed_casava_comment() {
        assert_eq!(PairingRule::select(b"read1 2:N:0:AT", b"read1 1:N:0:AT"), None);
    }

    #[test]
    fn select_ignores_non_mate_trailing_digit() {
        // A trailing "/3" isn't a recognized mate suffix.
        assert_eq!(PairingRule::select(b"read1/3", b"read1/4"), None);
    }

    // ---- PairingRule::check_pair ----

    #[test]
    fn check_pair_slash_digit_accepts_forward_order() {
        assert!(PairingRule::SlashDigit.check_pair(b"read1/1", b"read1/2"));
    }

    #[test]
    fn check_pair_slash_digit_rejects_reversed_order() {
        assert!(!PairingRule::SlashDigit.check_pair(b"read1/2", b"read1/1"));
    }

    #[test]
    fn check_pair_sep_digit_rejects_wrong_separator() {
        // Rule was selected as '.'-separated; a '_'-separated pair doesn't match it.
        assert!(!PairingRule::SepDigit(b'.').check_pair(b"read1_1", b"read1_2"));
    }

    #[test]
    fn check_pair_casava_or_bare_rejects_mismatched_stem() {
        assert!(!PairingRule::CasavaOrBare.check_pair(b"read1 1:N:0:AT", b"read2 2:N:0:AT"));
    }

    // ---- resolve_inputs / check_dash_at_most_once ----

    #[test]
    fn resolve_inputs_defaults_empty_to_dash() {
        assert_eq!(resolve_inputs(&[], false).unwrap(), vec![PathBuf::from("-")]);
    }

    #[test]
    fn resolve_inputs_leaves_explicit_paths_untouched() {
        let raw = vec![PathBuf::from("a.fq"), PathBuf::from("b.fq")];
        assert_eq!(resolve_inputs(&raw, false).unwrap(), raw);
    }

    #[test]
    fn resolve_inputs_errors_when_stdin_is_tty_and_defaulted() {
        assert!(resolve_inputs(&[], true).is_err());
    }

    #[test]
    fn resolve_inputs_errors_when_stdin_is_tty_and_explicit_dash() {
        assert!(resolve_inputs(&[PathBuf::from("-")], true).is_err());
    }

    #[test]
    fn resolve_inputs_ok_when_stdin_is_tty_but_dash_not_used() {
        let raw = vec![PathBuf::from("a.fq")];
        assert!(resolve_inputs(&raw, true).is_ok());
    }

    #[test]
    fn check_dash_at_most_once_allows_single_dash() {
        let mut errors = Vec::new();
        check_dash_at_most_once(&[PathBuf::from("-")], "Inputs", &mut errors);
        assert!(errors.is_empty());
    }

    #[test]
    fn check_dash_at_most_once_rejects_two_dashes() {
        let mut errors = Vec::new();
        check_dash_at_most_once(&[PathBuf::from("-"), PathBuf::from("-")], "Inputs", &mut errors);
        assert_eq!(errors.len(), 1);
    }

    fn reader_from(bytes: Vec<u8>) -> FastqReader<Box<dyn BufRead + Send>> {
        let boxed: Box<dyn BufRead + Send> = Box::new(std::io::Cursor::new(bytes));
        FastqReader::with_capacity(boxed, BUFFER_SIZE)
    }

    // ---- sniff_single_input ----

    #[test]
    fn sniff_single_input_empty_is_single_end_and_empty() {
        let s = sniff_single_input(reader_from(Vec::new())).unwrap();
        assert!(!s.interleaved);
        assert!(s.is_empty);
        assert_eq!(s.pairing_rule, None);
        assert!(s.records.count() == 0);
    }

    #[test]
    fn sniff_single_input_one_record_is_single_end() {
        let bytes = test_fastq_bytes(&[("read1", "ACGT")]);
        let s = sniff_single_input(reader_from(bytes)).unwrap();
        assert!(!s.interleaved);
        assert!(!s.is_empty);
        assert_eq!(s.records.count(), 1);
    }

    #[test]
    fn sniff_single_input_two_unrelated_reads_is_single_end() {
        let bytes = test_fastq_bytes(&[("read1", "ACGT"), ("read2", "ACGT")]);
        let s = sniff_single_input(reader_from(bytes)).unwrap();
        assert!(!s.interleaved);
        assert_eq!(s.records.count(), 2);
    }

    #[test]
    fn sniff_single_input_mate_pair_is_interleaved() {
        let bytes = test_fastq_bytes(&[("read1 1:N:0:AT", "ACGT"), ("read1 2:N:0:AT", "TGCA")]);
        let s = sniff_single_input(reader_from(bytes)).unwrap();
        assert!(s.interleaved);
        assert_eq!(s.pairing_rule, Some(PairingRule::CasavaOrBare));
        assert_eq!(s.records.count(), 2);
    }

    #[test]
    fn sniff_single_input_replays_peeked_records_in_order() {
        let bytes =
            test_fastq_bytes(&[("read1/1", "ACGT"), ("read1/2", "TGCA"), ("read2/1", "AAAA")]);
        let s = sniff_single_input(reader_from(bytes)).unwrap();
        assert!(s.interleaved);
        let heads: Vec<Vec<u8>> = s.records.map(|r| r.unwrap().head).collect();
        assert_eq!(heads, vec![b"read1/1".to_vec(), b"read1/2".to_vec(), b"read2/1".to_vec()]);
    }

    #[test]
    fn sniff_single_input_sra_interleaved_sniffs_as_pe() {
        // fastq-dump -I --split-spot style: `.1.1`/`.1.2`/`.2.1`/`.2.2`.
        let bytes = test_fastq_bytes(&[
            ("SRR1.1.1", "ACGT"),
            ("SRR1.1.2", "TGCA"),
            ("SRR1.2.1", "AAAA"),
            ("SRR1.2.2", "TTTT"),
        ]);
        let s = sniff_single_input(reader_from(bytes)).unwrap();
        assert!(s.interleaved);
        assert_eq!(s.pairing_rule, Some(PairingRule::SepDigit(b'.')));
    }

    #[test]
    fn sniff_single_input_sra_se_sniffs_as_se() {
        // Standard SE SRA naming: `.1`/`.2`/`.3`/`.4` — records 1-2 look like a
        // `.`-suffixed pair, but records 3-4 fail the same rule's orientation check.
        let bytes = test_fastq_bytes(&[
            ("SRR1.1", "ACGT"),
            ("SRR1.2", "TGCA"),
            ("SRR1.3", "AAAA"),
            ("SRR1.4", "TTTT"),
        ]);
        let s = sniff_single_input(reader_from(bytes)).unwrap();
        assert!(!s.interleaved);
        assert_eq!(s.pairing_rule, None);
        assert_eq!(s.records.count(), 4);
    }

    #[test]
    fn sniff_single_input_underscore_suffix_sniffs_as_pe() {
        let bytes = test_fastq_bytes(&[
            ("read_1", "ACGT"),
            ("read_2", "TGCA"),
            ("read2_1", "AAAA"),
            ("read2_2", "TTTT"),
        ]);
        let s = sniff_single_input(reader_from(bytes)).unwrap();
        assert!(s.interleaved);
        assert_eq!(s.pairing_rule, Some(PairingRule::SepDigit(b'_')));
    }

    #[test]
    fn sniff_single_input_reversed_pair_sniffs_as_se() {
        let bytes = test_fastq_bytes(&[("read1/2", "ACGT"), ("read1/1", "TGCA")]);
        let s = sniff_single_input(reader_from(bytes)).unwrap();
        assert!(!s.interleaved);
    }

    #[test]
    fn sniff_single_input_rule_confirmation_failure_at_3_4_sniffs_as_se() {
        // Records 1-2 pair under SlashDigit; records 3-4 don't (unrelated names).
        let bytes = test_fastq_bytes(&[
            ("pair0/1", "ACGT"),
            ("pair0/2", "TGCA"),
            ("unrelated_a", "AAAA"),
            ("unrelated_b", "TTTT"),
        ]);
        let s = sniff_single_input(reader_from(bytes)).unwrap();
        assert!(!s.interleaved);
    }

    #[test]
    fn sniff_single_input_short_pe_file_with_fewer_than_4_records_is_interleaved() {
        // Only 2 records total; a matching rule on 1-2 alone is enough (no records
        // 3-4 to confirm against).
        let bytes = test_fastq_bytes(&[("pair0/1", "ACGT"), ("pair0/2", "TGCA")]);
        let s = sniff_single_input(reader_from(bytes)).unwrap();
        assert!(s.interleaved);
    }

    // ---- open_one_fastq_input ----

    #[test]
    fn opener_detects_gzip_content_by_magic_bytes() {
        let tmp = tempfile::TempDir::new().unwrap();
        let path = tmp.path().join("misnamed.txt"); // deliberately not *.gz
        std::fs::write(&path, test_gzip(b"@r\nACGT\n+\nIIII\n")).unwrap();

        let mut reader = open_one_fastq_input(&path).unwrap();
        let mut first = [0u8; 1];
        reader.read_exact(&mut first).unwrap();
        assert_eq!(&first, b"@"); // decoded through MultiGzDecoder, not raw gzip bytes
    }

    #[test]
    fn opener_reads_plain_content() {
        let tmp = tempfile::TempDir::new().unwrap();
        let path = tmp.path().join("plain.fq");
        std::fs::write(&path, b"@r\nACGT\n+\nIIII\n").unwrap();

        let mut reader = open_one_fastq_input(&path).unwrap();
        let mut buf = String::new();
        std::io::Read::read_to_string(&mut reader, &mut buf).unwrap();
        assert_eq!(buf, "@r\nACGT\n+\nIIII\n");
    }

    #[test]
    fn opener_handles_empty_file() {
        let tmp = tempfile::TempDir::new().unwrap();
        let path = tmp.path().join("empty.fq");
        std::fs::write(&path, b"").unwrap();

        let mut reader = open_one_fastq_input(&path).unwrap();
        let mut buf = Vec::new();
        std::io::Read::read_to_end(&mut reader, &mut buf).unwrap();
        assert!(buf.is_empty());
    }

    #[test]
    fn opener_handles_one_byte_file() {
        // A single byte can never be gzip (magic is 2 bytes); must be read back as
        // plain content, not dropped or misdetected.
        let tmp = tempfile::TempDir::new().unwrap();
        let path = tmp.path().join("one_byte.fq");
        std::fs::write(&path, b"@").unwrap();

        let mut reader = open_one_fastq_input(&path).unwrap();
        let mut buf = Vec::new();
        std::io::Read::read_to_end(&mut reader, &mut buf).unwrap();
        assert_eq!(buf, b"@");
    }

    /// A `Read` that dribbles out 1 byte per call, regardless of the caller's buffer
    /// size — exercises the 2-byte gzip-sniff loop against a source that can't
    /// satisfy it in one `read()` (as a real pipe legitimately can't either).
    struct OneByteAtATime(std::io::Cursor<Vec<u8>>);

    impl Read for OneByteAtATime {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            if buf.is_empty() {
                return Ok(0);
            }
            self.0.read(&mut buf[..1])
        }
    }

    #[test]
    fn gzip_sniff_loop_accumulates_across_short_reads() {
        // Exercises the accumulation loop directly (rather than through
        // `open_one_fastq_input`, which is file/stdin-only): a source yielding 1
        // byte per `read()` call must still be correctly identified as gzip once 2
        // bytes have accumulated.
        let data = test_gzip(b"@r\nACGT\n+\nIIII\n");
        let mut src = OneByteAtATime(std::io::Cursor::new(data));
        let mut probe = [0u8; 2];
        let mut probed = 0usize;
        while probed < probe.len() {
            let n = src.read(&mut probe[probed..]).unwrap();
            if n == 0 {
                break;
            }
            probed += n;
        }
        assert_eq!(probed, 2);
        assert_eq!(probe, GZIP_MAGIC);
    }
}
