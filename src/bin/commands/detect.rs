//! The `chelae detect` subcommand: identify the adapter sequences present in one or
//! two FASTQ files by sampling a modest number of records and tallying matches.
//!
//! Two modes share the same sampler-and-reporter shell:
//!
//! * **Paired-end** runs the same overlap probe as `chelae trim` but with an empty
//!   adapter-evidence library — every probe-accepted shift becomes a candidate
//!   adapter detection regardless of whether the post-cut bases look like a known
//!   adapter. The post-template tails on each mate are then summarized as fixed-
//!   length k-mers and tallied. Near-identical k-mers are merged via Hamming-
//!   distance aggregation before reporting, so single-base sequencing errors don't
//!   fragment the count of the true adapter.
//!
//! * **Single-end** runs each candidate adapter (every built-in kit plus any user
//!   `--adapter-sequence` / `--adapter-fasta` entry) against each read's 3' end via
//!   the same scanner `chelae trim` uses. The candidate with the longest matched
//!   overhang on a given read wins that read's vote.
//!
//! In both modes, sampling stops once `--num-detections` usable detections have
//! accumulated or `--max-reads` records have been read, whichever comes first. Any
//! adapter sequence accounting for at least `--min-fraction` of usable detections
//! is reported on stdout, and (optionally) written as a FASTA file ready to feed
//! back into `chelae trim --adapter-fasta`.

use crate::commands::command::Command;
use crate::commands::trim::{
    Adapter, OverlapAdapterLibrary, OverlapStats, QualityTrim, count_mismatches_ci_bounded,
    cut_right_quality_position, detect_pe_overlap, find_adapter_3prime, find_polyx_tail_len,
    load_adapter_fasta_with_names, validate_adapter_bases,
};
use crate::commands::utils::{BUFFER_SIZE, fmt_count, open_fastq_inputs};
use anyhow::{Result, anyhow};
use chelae_lib::adapter_db::ALL_KITS;
use clap::Parser;
use fgoxide::io::Io;
use log::{info, warn};
use seq_io::fastq::{Reader as FastqReader, Record};
use std::collections::HashMap;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};
use std::str::FromStr;

/// Length of the k-mer summary taken from each PE post-template tail. 16 bp is long
/// enough to separate every adapter in [`ALL_KITS`] cleanly (their first 16 bases
/// are pairwise distinct) yet short enough that the post-collection Hamming-merge
/// pass folds single-base sequencing errors back into the true bucket.
const TAIL_KMER_LEN: usize = 16;

/// Maximum Hamming distance for two tail k-mers to be considered the same adapter
/// during aggregation. At 16 bp the per-base sequencing-error rate is ~1%, so two
/// errors in the window is realistic; the kit-to-kit minimum distance over the
/// first 16 bp is well above this, so kits never collide.
const FUZZY_MERGE_HAMMING: usize = 2;

/// Maximum Hamming distance allowed when declaring a "near match" of a discovered
/// k-mer to a known kit adapter's first [`TAIL_KMER_LEN`] bp. 1 (not 2) because a
/// kit-match triggers FASTA substitution — emitting the kit's canonical sequence
/// in place of the in-sample consensus — and 2 mm over 16 bp is too loose a bar
/// for that substitution. A novel adapter sharing 14/16 bp with a kit prefix
/// should not be silently relabeled as that kit. Real sequencing-error fuzzy
/// matches almost always sit at 1 mm given the per-base error rate.
const KIT_NEAR_MATCH_HAMMING: usize = 1;

/// ASCII bases corresponding to slots in [`TailAccumulator::base_counts`].
/// `BASES[i]` is the byte to emit when slot `i` (A=0, C=1, G=2, T=3, N=4)
/// wins the per-position consensus vote.
const BASES: [u8; 5] = [b'A', b'C', b'G', b'T', b'N'];

/// Absolute minimum majority fraction required to emit a base into the
/// consensus. Below this the column is treated as genuinely ambiguous
/// (barcode / soft-clipped region / etc.) and the consensus is truncated
/// regardless of prior columns. 0.50 == "one base must at least be more
/// common than all others combined."
const CONSENSUS_MAJORITY_FLOOR: f64 = 0.50;

/// Drop in per-column majority fraction below the running-minimum baseline
/// that triggers a discontinuity cut in the consensus. The intuition: the
/// kit-stable prefix rocks along at 90+% majority; the first variable-region
/// position (barcode, sample index, primer landing pad) collapses to a much
/// lower majority — even if a single index in the pool happens to dominate,
/// the drop from the prefix baseline is large. Cutting at the discontinuity
/// keeps the FASTA cross-sample-portable rather than emitting one sample's
/// plurality barcode base as if it were canonical.
const CONSENSUS_DROP_TOLERANCE: f64 = 0.10;

/// Identify the adapter sequence(s) present in a FASTQ file.
///
/// Paired-end input uses R1/R2 overlap detection to locate the adapter start in
/// each pair and harvests the post-template bases as candidate adapter k-mers, so
/// novel adapters are discovered even with no prior knowledge of the chemistry.
/// Single-end input scores each read against every built-in kit adapter plus any
/// candidate(s) supplied via `--adapter-sequence` / `--adapter-fasta`.
///
/// In both modes, sampling stops once `--num-detections` usable detections have
/// been observed (or `--max-reads` records have been scanned, whichever comes
/// first). Distinct adapter sequences accounting for at least `--min-fraction` of
/// usable detections are reported. A FASTA file with the winning adapters can be
/// written via `--output-fasta` for direct re-use by `chelae trim --adapter-fasta`.
///
/// # Example
///
/// ```bash
/// # Paired-end discovery
/// chelae detect -i r1.fq.gz r2.fq.gz -o adapters.fa
///
/// # Single-end against the built-in kits
/// chelae detect -i reads.fq.gz
/// ```
#[derive(Parser, Debug)]
#[command(version)]
#[clap(verbatim_doc_comment)]
pub(crate) struct Detect {
    /// One or two input FASTQ files. A single path is single-end; two paths are paired-end.
    /// Inputs may be plain, gzip, or bgzf (auto-detected).
    #[clap(long, short = 'i', required = true, num_args = 1..=2)]
    inputs: Vec<PathBuf>,

    /// Optional FASTA output path. Discovered/winning adapter sequence(s) are written
    /// here, in a format ready to feed back into `chelae trim --adapter-fasta`.
    #[clap(long, short = 'o')]
    output_fasta: Option<PathBuf>,

    /// (SE only) Extra candidate adapter sequence(s) to score against, in addition to
    /// every built-in kit. Pass repeatedly for multiple candidates (`-a AAA -a CCC`).
    /// Errors if supplied on PE input.
    #[clap(long, short = 'a')]
    adapter_sequence: Vec<String>,

    /// (SE only) FASTA file of extra candidate adapter sequences. Each record's name is
    /// preserved in the report. Errors if supplied on PE input.
    #[clap(long, short = 'f')]
    adapter_fasta: Option<PathBuf>,

    /// Stop sampling once this many usable detections (pairs/reads where an adapter
    /// signal was observed with at least `--min-tail-length` bases of evidence) have
    /// accumulated. Higher values give more confident estimates of mixture composition.
    #[clap(long, short = 'n', default_value = "5000")]
    num_detections: u64,

    /// Hard cap on the number of input records scanned even if `--num-detections` is
    /// not reached. Protects against inputs with little or no adapter readthrough.
    #[clap(long, default_value = "1000000")]
    max_reads: u64,

    /// Minimum number of usable detections required to emit a report. If sampling
    /// stops below this floor (either via EOF, `--max-reads`, or both), `detect`
    /// errors out rather than producing a confident-looking summary from a tiny
    /// sample. Lower it on small/sparse libraries where you intentionally want a
    /// best-effort answer; raise it to demand more evidence before reporting.
    #[clap(long, default_value = "20")]
    min_detections_for_report: u64,

    /// Minimum fraction of usable detections (0.0 - 1.0) that an adapter must account
    /// for to be reported. Lower it to surface minor contaminants; raise it to hide
    /// noise.
    #[clap(long, default_value = "0.05")]
    min_fraction: f64,

    /// Minimum length of adapter evidence per detection (bp). For PE, the post-
    /// template tail on each mate must reach this length; for SE, the matched
    /// adapter alignment — `min(read_len - trim_pos, candidate_adapter_len)` —
    /// must reach this length. Shorter evidence is too fragile to reliably
    /// identify the adapter.
    #[clap(long, default_value = "8")]
    min_tail_length: usize,

    /// (PE only) Minimum overlap length (bp) required for PE-overlap detection.
    #[clap(long, default_value = "30")]
    overlap_min_length: usize,

    /// (PE only) Maximum fraction of mismatches permitted in the PE-overlap probe
    /// (0.0 - 1.0).
    #[clap(long, default_value = "0.10")]
    overlap_max_mismatch_rate: f64,

    /// (PE only) Upper bound on the probe length per overlap-length candidate (bp).
    #[clap(long, default_value = "64")]
    overlap_diagnostic_length: usize,

    /// (SE only) Minimum match length (bp) when scoring a candidate adapter against
    /// a read's 3' end. Independent of `--min-tail-length`, which gates how many
    /// matched bases are required before the detection is counted.
    ///
    /// Default 10 (higher than `chelae trim`'s 6) because detect scores every
    /// built-in kit's adapter plus any user-supplied candidate against each read.
    /// At a 6 bp floor the per-position random-match rate (~1/4096 for ACGT) over
    /// ~50 candidate start positions × 8+ candidates leaves a noticeable false-
    /// positive floor; 10 bp drops that ~16-fold.
    #[clap(long, default_value = "10")]
    adapter_min_length: usize,

    /// (SE only) Maximum fraction of mismatches (0.0 - 1.0) permitted when matching
    /// an adapter against a read's 3' end.
    #[clap(long, default_value = "0.125")]
    adapter_mismatch_rate: f64,

    /// 3' poly-G trim minimum run length. Matches `chelae trim`'s default: on by
    /// default to clean Illumina 2-color chemistry artifacts (G is the "no signal"
    /// call) that would otherwise corrupt the overlap probe. Pass `0` to disable.
    #[clap(long, default_value = "10")]
    trim_polyg: usize,

    /// 3' poly-X trim minimum run length (trims homopolymer A/C/T tails). On by
    /// default in detect with a more aggressive default than `chelae trim`'s
    /// opt-in 10: detect's goal is specificity, and short homopolymer tails (e.g.
    /// short polyA contamination) can corrupt the overlap probe without
    /// hurting the discovery of real adapter sequence. Pass `0` to disable.
    #[clap(long, default_value = "5")]
    trim_polyx: usize,

    /// Sliding-window 3' quality trim as `WINDOW:QUAL` with cut-right semantics:
    /// scans 5'→3' and truncates the read at the start of the first window whose
    /// mean Phred quality falls below the threshold. On by default in detect at
    /// `4:20` — a tighter window than `chelae trim`'s opt-in `--quality-trim-3p`
    /// default of `8:20`, so quality-degraded 3' tails are cleaned before the
    /// overlap/adapter probe sees them. Pass `off`, `none`, or `no` (case-
    /// insensitive) to disable.
    #[clap(long, default_value = "4:20")]
    quality_trim: QualityTrimSetting,
}

impl Detect {
    /// Validates CLI inputs, aggregating every problem into one error message.
    fn validate(&self) -> Result<()> {
        let mut errors: Vec<String> = Vec::new();

        for path in &self.inputs {
            if !path.exists() {
                errors.push(format!("Input file {path:?} does not exist."));
            }
        }

        if let Some(path) = &self.output_fasta
            && let Some(parent) = path.parent()
            && !parent.as_os_str().is_empty()
            && !parent.exists()
        {
            errors.push(format!(
                "Output FASTA parent directory {parent:?} does not exist (for {path:?})."
            ));
        }

        // --adapter-sequence / --adapter-fasta are SE-only knobs; warn-then-ignore is
        // worse UX than refusing, since a user passing them on PE almost certainly
        // expects them to do something.
        if self.inputs.len() == 2 {
            if !self.adapter_sequence.is_empty() {
                errors.push(
                    "--adapter-sequence is not used in paired-end mode (PE detection \
                     discovers adapters via overlap)."
                        .to_string(),
                );
            }
            if self.adapter_fasta.is_some() {
                errors.push(
                    "--adapter-fasta is not used in paired-end mode (PE detection \
                     discovers adapters via overlap)."
                        .to_string(),
                );
            }
        }

        for seq in &self.adapter_sequence {
            if seq.is_empty() {
                errors.push("--adapter-sequence values must not be empty.".to_string());
                continue;
            }
            if let Err(msg) = validate_adapter_bases(seq.as_bytes()) {
                errors.push(format!("--adapter-sequence {seq:?}: {msg}"));
            }
            // A candidate shorter than --adapter-min-length can never satisfy
            // find_adapter_3prime's alignment-length guard, so it would
            // silently contribute zero hits — the user would see "no
            // candidate matched" with no clue that their sequence was
            // structurally unmatchable.
            if seq.len() < self.adapter_min_length {
                errors.push(format!(
                    "--adapter-sequence {seq:?} is {} bp, shorter than --adapter-min-length ({}); \
                     it could never match. Lower --adapter-min-length or supply a longer sequence.",
                    seq.len(),
                    self.adapter_min_length,
                ));
            }
        }

        if !(0.0..=1.0).contains(&self.min_fraction) {
            errors.push(format!("--min-fraction must be in 0.0..=1.0, got {}.", self.min_fraction));
        }
        if !(0.0..=1.0).contains(&self.overlap_max_mismatch_rate) {
            errors.push(format!(
                "--overlap-max-mismatch-rate must be in 0.0..=1.0, got {}.",
                self.overlap_max_mismatch_rate
            ));
        }
        if !(0.0..=1.0).contains(&self.adapter_mismatch_rate) {
            errors.push(format!(
                "--adapter-mismatch-rate must be in 0.0..=1.0, got {}.",
                self.adapter_mismatch_rate
            ));
        }
        if self.num_detections == 0 {
            errors.push("--num-detections must be at least 1.".to_string());
        }
        if self.max_reads == 0 {
            errors.push("--max-reads must be at least 1.".to_string());
        }
        if self.min_detections_for_report == 0 {
            errors.push("--min-detections-for-report must be at least 1.".to_string());
        }
        if self.min_detections_for_report > self.num_detections {
            errors.push(format!(
                "--min-detections-for-report ({}) cannot exceed --num-detections ({}) — the \
                 sampler would never reach the floor.",
                self.min_detections_for_report, self.num_detections,
            ));
        }
        if self.min_detections_for_report > self.max_reads {
            errors.push(format!(
                "--min-detections-for-report ({}) cannot exceed --max-reads ({}) — the \
                 sampler would scan all input and still bail under the floor.",
                self.min_detections_for_report, self.max_reads,
            ));
        }
        if self.min_tail_length == 0 {
            errors.push("--min-tail-length must be at least 1.".to_string());
        }
        if self.overlap_min_length == 0 {
            errors.push("--overlap-min-length must be at least 1.".to_string());
        }
        if self.overlap_diagnostic_length == 0 {
            errors.push("--overlap-diagnostic-length must be at least 1.".to_string());
        }
        if self.adapter_min_length == 0 {
            errors.push("--adapter-min-length must be at least 1.".to_string());
        }

        if errors.is_empty() {
            Ok(())
        } else {
            use std::fmt::Write;
            let detail = errors.iter().fold(String::new(), |mut s, e| {
                let _ = writeln!(s, "    - {e}");
                s
            });
            Err(anyhow!("Input validation failed:\n{detail}"))
        }
    }

    /// Paired-end discovery loop. Runs the standard overlap walk with an empty
    /// evidence library, harvests the post-template tail on each mate at the
    /// detected insert size, and bumps per-mate k-mer counters. Sampling stops at
    /// the first of: `--num-detections` usable detections, `--max-reads` records
    /// scanned, or EOF on either input.
    fn run_pe(
        &self,
        mut reader1: FastqReader<Box<dyn BufRead + Send>>,
        mut reader2: FastqReader<Box<dyn BufRead + Send>>,
    ) -> Result<()> {
        // Empty library disables the adapter-evidence side-check in the overlap
        // walk: every probe-accepted shift becomes a candidate detection regardless
        // of whether the post-cut bases look like a known adapter. This is what
        // lets detect discover novel adapters, and the trade-off is that an
        // overlap-acceptance for a low-complexity template that self-aligns is
        // accepted on probe alone. With `--overlap-min-length 30` and the default
        // 10% mismatch rate, false acceptances are rare and the fuzzy-merge step
        // further dilutes them; pathological inputs (heavy poly-A, repeat tracts)
        // could still pollute the harvested k-mer counts.
        let empty_lib = OverlapAdapterLibrary::default();
        // `OverlapStats` is created with no hint and `stats_on=false` is passed to
        // `observe` below. The histogram/unknown counters therefore stay unused —
        // the only reason we maintain this state is so the walk's `center_shift`
        // self-tunes toward the running-mean insert after a warm-up, which keeps
        // per-pair probe iteration counts low on long-insert libraries.
        let mut stats = OverlapStats::new(None);
        let mut rc_scratch: Vec<u8> = Vec::new();

        let mut r1_kmers: HashMap<Vec<u8>, TailAccumulator> = HashMap::new();
        let mut r2_kmers: HashMap<Vec<u8>, TailAccumulator> = HashMap::new();
        let mut detections: u64 = 0;
        let mut reads_scanned: u64 = 0;
        let mut overlap_hits: u64 = 0;

        loop {
            if detections >= self.num_detections || reads_scanned >= self.max_reads {
                break;
            }
            let r1 = match reader1.next() {
                Some(Ok(rec)) => rec,
                Some(Err(e)) => return Err(anyhow!("R1 FASTQ read error: {e}")),
                // R1 EOF: confirm R2 is also at EOF; otherwise the inputs are
                // out of sync and we want to surface that explicitly rather
                // than silently accepting the truncation.
                None => match reader2.next() {
                    None => break,
                    Some(Ok(_)) => {
                        return Err(anyhow!("R1 exhausted before R2 (inputs out of sync)"));
                    }
                    Some(Err(e)) => return Err(anyhow!("R2 FASTQ read error: {e}")),
                },
            };
            let r2 = match reader2.next() {
                Some(Ok(rec)) => rec,
                Some(Err(e)) => return Err(anyhow!("R2 FASTQ read error: {e}")),
                None => return Err(anyhow!("R2 exhausted before R1 (inputs out of sync)")),
            };
            reads_scanned += 1;

            // Apply the 3'-end cleanups (poly-G / poly-X / quality) before the
            // overlap probe sees the reads. detect's bias is specificity over
            // sensitivity — these scrub the noisy tails most likely to corrupt
            // the probe (Illumina 2-color G-runs, short polyA contamination,
            // quality-degraded ends) without touching the adapter region.
            let r1_full = r1.seq();
            let r2_full = r2.seq();
            let r1_end = effective_trimmed_len(
                r1_full,
                r1.qual(),
                self.trim_polyg,
                self.trim_polyx,
                self.quality_trim.as_option(),
            );
            let r2_end = effective_trimmed_len(
                r2_full,
                r2.qual(),
                self.trim_polyg,
                self.trim_polyx,
                self.quality_trim.as_option(),
            );
            let r1_seq = &r1_full[..r1_end];
            let r2_seq = &r2_full[..r2_end];
            let center = stats.center_shift(r2_seq.len());
            let result = detect_pe_overlap(
                r1_seq,
                r2_seq,
                self.overlap_min_length,
                self.overlap_max_mismatch_rate,
                self.overlap_diagnostic_length,
                &empty_lib,
                center,
                false,
                &mut rc_scratch,
            );
            stats.observe(result, false);

            let Some(insert) = result.inferred_insert else { continue };
            overlap_hits += 1;
            // Both mates must have a post-template tail of at least --min-tail-length
            // bp for the pair to count as a detection. `r1_kmers` / `r2_kmers` are
            // pushed together per detection (one atomic increment to the shared
            // `detections` denominator), so an asymmetric push — R2 only, say —
            // would corrupt the per-mate fraction calculations in
            // filter_min_fraction. The two branches below therefore both use
            // logical-OR: any one failing drops the whole pair.
            if insert >= r1_seq.len() || insert >= r2_seq.len() {
                continue;
            }
            let r1_tail = &r1_seq[insert..];
            let r2_tail = &r2_seq[insert..];
            if r1_tail.len() < self.min_tail_length || r2_tail.len() < self.min_tail_length {
                continue;
            }

            push_tail_kmer(&mut r1_kmers, r1_tail);
            push_tail_kmer(&mut r2_kmers, r2_tail);
            detections += 1;
        }

        if detections == 0 {
            return Err(anyhow!(
                "No usable PE adapter detections in {} pair(s) scanned ({} overlap hits, but \
                 none with a post-template tail of >= {} bp on both mates). The library may \
                 have inserts longer than the read length on every pair (no adapter \
                 readthrough), or the library may have very low readthrough — try increasing \
                 `--max-reads`, lowering `--overlap-min-length`, or (if overlap hits > 0) \
                 lowering `--min-tail-length`.",
                fmt_count(reads_scanned),
                fmt_count(overlap_hits),
                self.min_tail_length,
            ));
        }
        if detections < self.min_detections_for_report {
            return Err(anyhow!(
                "Only {} usable detections after scanning {} pair(s); below the \
                 `--min-detections-for-report` floor of {}. Reporting fractions from this \
                 small a sample would mislead. Increase `--max-reads`, lower \
                 `--min-detections-for-report` if you accept a noisier estimate, or rerun on \
                 a larger input.",
                fmt_count(detections),
                fmt_count(reads_scanned),
                fmt_count(self.min_detections_for_report),
            ));
        }
        if detections < self.num_detections {
            warn!(
                "Reached EOF / --max-reads with only {} of the requested {} detections; \
                 reported fractions may be noisier than expected.",
                fmt_count(detections),
                fmt_count(self.num_detections),
            );
        }

        let r1_hits = aggregate_kmers(r1_kmers, FUZZY_MERGE_HAMMING);
        let r2_hits = aggregate_kmers(r2_kmers, FUZZY_MERGE_HAMMING);
        let r1_reported = filter_min_fraction(&r1_hits, detections, self.min_fraction);
        let r2_reported = filter_min_fraction(&r2_hits, detections, self.min_fraction);
        let r1_annot = annotate_pe_hits(&r1_reported, KitMate::R1);
        let r2_annot = annotate_pe_hits(&r2_reported, KitMate::R2);

        emit_pe_report(detections, reads_scanned, overlap_hits, &r1_annot, &r2_annot);

        if let Some(path) = &self.output_fasta {
            // The console report shows what we found (or didn't); a FASTA
            // that would be silently half-empty (R1 hits, no R2) or fully
            // empty (nothing above --min-fraction on either mate) is worse
            // than no FASTA at all — downstream `chelae trim --adapter-fasta`
            // would run against an incomplete/absent adapter list and
            // silently under-trim. Fail loudly so the user can either lower
            // --min-fraction or drop --output-fasta and inspect the report.
            let missing = match (r1_reported.is_empty(), r2_reported.is_empty()) {
                (true, true) => Some("both R1 and R2"),
                (true, false) => Some("R1"),
                (false, true) => Some("R2"),
                (false, false) => None,
            };
            if let Some(which) = missing {
                return Err(anyhow!(
                    "No adapter reached --min-fraction ({}) on {which}. The console report \
                     above shows what was found. Lower --min-fraction and rerun to write a \
                     FASTA, or drop --output-fasta if you only want the diagnostic.",
                    self.min_fraction,
                ));
            }
            write_fasta(path, &pe_fasta_records(&r1_annot, &r2_annot))?;
            info!("Wrote discovered adapter FASTA to {path:?}");
        }

        Ok(())
    }

    /// Single-end scoring loop. Each read is run against every kit + user candidate
    /// via [`find_adapter_3prime`]; the candidate with the longest matched overhang
    /// (smallest trim position) wins that read's vote. A read with no candidate
    /// match contributes nothing.
    fn run_se(&self, mut reader: FastqReader<Box<dyn BufRead + Send>>) -> Result<()> {
        let candidates = build_se_candidates(
            &self.adapter_sequence,
            &self.adapter_fasta,
            self.adapter_min_length,
        )?;
        // `build_se_candidates` always seeds with every entry in `ALL_KITS`, which
        // is non-empty by construction, so this branch is unreachable at runtime.
        debug_assert!(!candidates.is_empty(), "candidate list is empty; ALL_KITS broken?");

        let mut counts = vec![0u64; candidates.len()];
        let mut detections: u64 = 0;
        let mut reads_scanned: u64 = 0;

        loop {
            if detections >= self.num_detections || reads_scanned >= self.max_reads {
                break;
            }
            let rec = match reader.next() {
                Some(Ok(r)) => r,
                Some(Err(e)) => return Err(anyhow!("FASTQ read error: {e}")),
                None => break,
            };
            reads_scanned += 1;
            // 3'-end cleanups before the candidate scan. Same rationale as PE:
            // a noisy 3' tail (poly-G, polyA, quality dropout) can produce
            // spurious "winning" matches to a candidate's prefix; we'd rather
            // see no detection than a wrong one given detect's specificity
            // bias.
            let full = rec.seq();
            let end = effective_trimmed_len(
                full,
                rec.qual(),
                self.trim_polyg,
                self.trim_polyx,
                self.quality_trim.as_option(),
            );
            let seq = &full[..end];
            // Pick the candidate with the longest matched overhang (smallest trim
            // position `k`). On exact ties (`k_new == k_best`) we keep the prior
            // best — i.e. the candidate that appears earliest in `candidates`,
            // which is kit-insertion order. This is deterministic but introduces a
            // mild bias toward earlier kits; in practice the kit adapters in
            // ALL_KITS have distinct first-bases so ties at the same `k` are rare.
            let mut best: Option<(usize, usize)> = None;
            for (idx, cand) in candidates.iter().enumerate() {
                if let Some(k) = find_adapter_3prime(
                    seq,
                    &cand.adapter,
                    self.adapter_min_length,
                    self.adapter_mismatch_rate,
                    None,
                ) && best.is_none_or(|(_, bk)| k < bk)
                {
                    best = Some((idx, k));
                }
            }
            if let Some((idx, k)) = best {
                // Matched alignment length is capped by both the available read tail
                // and the candidate adapter's own length — `seq.len() - k` alone
                // overstates evidence when the read tail extends past the adapter's
                // end and the trailing bases happen to be unrelated.
                let matched_len = (seq.len() - k).min(candidates[idx].adapter.bytes.len());
                if matched_len >= self.min_tail_length {
                    counts[idx] += 1;
                    detections += 1;
                }
            }
        }

        if detections == 0 {
            return Err(anyhow!(
                "No candidate adapter matched any read in {} record(s) scanned. The library \
                 may be free of read-through, or the configured candidates may not include \
                 the adapter actually present — try supplying it via `--adapter-sequence` or \
                 `--adapter-fasta`, or (if you have paired-end reads) rerun in PE mode where \
                 detect discovers novel adapters via overlap.",
                fmt_count(reads_scanned),
            ));
        }
        if detections < self.min_detections_for_report {
            return Err(anyhow!(
                "Only {} candidate-matched reads after scanning {} record(s); below the \
                 `--min-detections-for-report` floor of {}. Reporting fractions from this \
                 small a sample would mislead. Increase `--max-reads`, lower \
                 `--min-detections-for-report` if you accept a noisier estimate, or rerun on \
                 a larger input.",
                fmt_count(detections),
                fmt_count(reads_scanned),
                fmt_count(self.min_detections_for_report),
            ));
        }
        if detections < self.num_detections {
            warn!(
                "Reached EOF / --max-reads with only {} of the requested {} detections.",
                fmt_count(detections),
                fmt_count(self.num_detections),
            );
        }

        // Build a Hit list directly from the per-candidate counters (no fuzzy merge
        // needed — candidates are already distinct named sequences).
        let mut hits: Vec<Hit> = candidates
            .iter()
            .zip(counts.iter())
            .filter(|&(_, &c)| c > 0)
            .map(|(c, &count)| Hit {
                name: Some(c.name.clone()),
                seq: c.adapter.bytes.clone(),
                count,
            })
            .collect();
        hits.sort_by_key(|h| std::cmp::Reverse(h.count));

        let reported = filter_min_fraction(&hits, detections, self.min_fraction);
        let annotated = annotate_se_hits(&reported);
        emit_se_report(detections, reads_scanned, &annotated);

        if let Some(path) = &self.output_fasta {
            // Fail loud rather than write an empty FASTA — same reasoning as
            // the PE branch: an empty adapters.fa fed back into
            // `chelae trim --adapter-fasta` silently under-trims.
            if reported.is_empty() {
                return Err(anyhow!(
                    "No candidate adapter reached --min-fraction ({}). The console report \
                     above shows the per-candidate counts. Lower --min-fraction and rerun to \
                     write a FASTA, or drop --output-fasta if you only want the diagnostic.",
                    self.min_fraction,
                ));
            }
            // SE candidates already carry the canonical kit (or user-supplied)
            // sequence — no per-sample consensus to substitute. Use the
            // candidate's name and bytes directly.
            let records: Vec<(String, &[u8])> = reported
                .iter()
                .map(|h| {
                    (h.name.clone().unwrap_or_else(|| "adapter".to_string()), h.seq.as_slice())
                })
                .collect();
            write_fasta(path, &records)?;
            info!("Wrote winning adapter FASTA to {path:?}");
        }

        Ok(())
    }
}

impl Command for Detect {
    /// Validates inputs, opens the FASTQ reader(s), and dispatches to the PE or SE
    /// detection loop. Both paths share the aggregation, reporting, and FASTA-output
    /// shell via free helpers in this module.
    fn execute(&self) -> Result<()> {
        self.validate()?;
        info!(
            "Detecting adapters in {} input file(s) (target {} detections, hard cap {} reads)",
            self.inputs.len(),
            fmt_count(self.num_detections),
            fmt_count(self.max_reads),
        );
        let mut readers = open_fastq_inputs(&self.inputs)?;
        match readers.len() {
            1 => self.run_se(readers.pop().unwrap()),
            2 => {
                let r2 = readers.pop().unwrap();
                let r1 = readers.pop().unwrap();
                self.run_pe(r1, r2)
            }
            // clap's `num_args = 1..=2` already enforces this, but be defensive.
            n => Err(anyhow!("Expected 1 or 2 inputs; got {n}.")),
        }
    }
}

/// `--quality-trim` setting. Distinguishes "off" (explicit opt-out) from "on
/// with these `WINDOW:QUAL` parameters". A dedicated enum (rather than
/// `Option<QualityTrim>` with a custom value-parser) sidesteps clap's
/// reflection limitations around `Option<T>` fields.
#[derive(Debug, Clone, Copy)]
enum QualityTrimSetting {
    /// Quality trimming disabled — read end is left at its full length.
    Off,
    /// Quality trimming enabled with the supplied window/threshold.
    On(QualityTrim),
}

impl QualityTrimSetting {
    /// Returns the wrapped `QualityTrim` parameters when enabled, `None` otherwise.
    /// Used to bridge the CLI-level [`QualityTrimSetting`] into the
    /// `Option<QualityTrim>` parameter that [`effective_trimmed_len`] expects.
    fn as_option(self) -> Option<QualityTrim> {
        match self {
            Self::Off => None,
            Self::On(qt) => Some(qt),
        }
    }
}

impl FromStr for QualityTrimSetting {
    type Err = String;

    /// Parses the `--quality-trim` flag value: accepts any of "off" / "none" /
    /// "no" (case-insensitive) as the explicit-disable form, otherwise delegates
    /// to [`QualityTrim`]'s `WINDOW:QUAL` parser.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let lower = s.to_ascii_lowercase();
        if matches!(lower.as_str(), "off" | "none" | "no") {
            Ok(Self::Off)
        } else {
            QualityTrim::from_str(s).map(Self::On)
        }
    }
}

/// One reported adapter hit after aggregation. `count` is the sum of detections
/// folded into this entry (including any near-identical sequences merged in via
/// [`FUZZY_MERGE_HAMMING`]). `name` is set for SE (kit / user adapter names) and
/// unset for PE (where the entry is a discovered k-mer with no a priori name).
#[derive(Debug, Clone)]
struct Hit {
    /// Display name for the hit. Set by the SE path to the matching candidate's
    /// name (kit name, `user_N`, or the FASTA record's own `>name`); `None` in the
    /// PE path where the hit is a sequence discovered de novo with no a priori
    /// label.
    name: Option<String>,
    /// Sequence used as the hit's identity. For PE this is the (uppercased) tail
    /// k-mer; for SE it is the full candidate adapter sequence.
    seq: Vec<u8>,
    /// Number of detections folded into this hit, including any near-identical
    /// sequences merged in during [`aggregate_kmers`].
    count: u64,
}

/// One candidate adapter for the SE scoring loop: a display name plus the
/// `Adapter` value the per-read scanner expects.
struct Candidate {
    /// Display name surfaced in the stdout report and FASTA output: a kit name,
    /// `<kit>_r2`, a synthetic `user_N` for `--adapter-sequence` entries, or the
    /// FASTA record's own `>name` (with synthetic `record_N` fallback) for entries
    /// loaded via `--adapter-fasta`.
    name: String,
    /// Wrapped adapter sequence the per-read scanner consumes; the `pure_acgt`
    /// flag is precomputed so [`find_adapter_3prime`] can pick the SIMD fast path.
    adapter: Adapter,
}

/// Identifies the relationship of a discovered or scored adapter to a known kit.
/// Used in the stdout report so the user can tell at a glance whether the result
/// hits a known kit or is novel, and in the FASTA writer so a kit-matched
/// discovery emits the kit's published sequence (the cross-sample-portable
/// artifact) rather than the in-sample consensus (which extends into the i7/i5
/// barcode region and is therefore sample-specific).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum KitMatch {
    /// The first [`TAIL_KMER_LEN`] bases of the discovered sequence are byte-for-
    /// byte equal to the named kit's adapter on the queried mate.
    Exact {
        /// Name of the matched kit (from `KitAdapter::name`).
        kit: &'static str,
        /// The kit's full published adapter sequence on the queried mate. The
        /// cross-sample-portable artifact: emitted into FASTA verbatim when this
        /// hit is a kit match.
        kit_seq: &'static [u8],
    },
    /// The discovered sequence is within [`KIT_NEAR_MATCH_HAMMING`] mismatches of
    /// the named kit's adapter on the queried mate, but not exact. Typical cause
    /// is sequencing error in the harvested tail.
    Fuzzy {
        /// Name of the closest kit.
        kit: &'static str,
        /// The kit's full published adapter sequence on the queried mate (same
        /// role as in [`KitMatch::Exact`]). Fuzzy matches still emit this
        /// verbatim in FASTA on the premise that mismatches are sequencing
        /// errors against the canonical, not real per-library variation.
        kit_seq: &'static [u8],
        /// Hamming-distance to the kit adapter's first [`TAIL_KMER_LEN`] bases.
        mismatches: usize,
    },
}

impl KitMatch {
    /// Returns the matched kit's name regardless of variant.
    fn kit(&self) -> &'static str {
        match self {
            KitMatch::Exact { kit, .. } | KitMatch::Fuzzy { kit, .. } => kit,
        }
    }

    /// Returns the matched kit's published adapter sequence regardless of variant.
    fn kit_seq(&self) -> &'static [u8] {
        match self {
            KitMatch::Exact { kit_seq, .. } | KitMatch::Fuzzy { kit_seq, .. } => kit_seq,
        }
    }
}

/// Which mate-side adapter list to compare against in [`classify_against_kits`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum KitMate {
    /// Compare against each kit's `seq_r1` adapter.
    R1,
    /// Compare against each kit's `seq_r2` adapter; kits whose `seq_r2` is `None`
    /// (e.g. small-RNA) are skipped.
    R2,
}

/// A reported hit bundled with its mate label (for PE; `None` in SE) and the
/// pre-computed kit-match annotation. Built once via [`annotate_pe_hits`] /
/// [`annotate_se_hits`] so the section-1 / section-2 console emitters and the
/// FASTA writer all share the same kit-match decision (and so the FASTA writer
/// can substitute the kit's published sequence consistently with what section 1
/// of the report announces).
struct AnnotatedHit<'a> {
    /// Underlying hit (consensus sequence + count + optional name).
    hit: &'a Hit,
    /// Mate label for the report ("R1" / "R2") in PE mode; `None` for SE.
    mate: Option<&'static str>,
    /// Kit-match annotation, or `None` if the hit doesn't resemble any kit
    /// adapter within `KIT_NEAR_MATCH_HAMMING` of the queried mate.
    kit_match: Option<KitMatch>,
}

/// Per-cluster accumulator for the PE-discovery path: a count of detections
/// folded into the bucket plus per-position base counts across every tail
/// observed in that bucket. The first [`TAIL_KMER_LEN`] positions are guaranteed
/// to be present (that's what the hashmap is keyed on); later positions are
/// populated only by reads whose tail extended that far, so coverage tails off
/// with read-tail length distribution. Consensus extraction stops where coverage
/// falls below the adaptive floor — see [`Self::consensus`].
#[derive(Default, Debug, Clone)]
struct TailAccumulator {
    /// Number of detections folded into this bucket.
    count: u64,
    /// Per-position base counts; outer index is the position in the tail
    /// (0-based from the adapter start), inner index follows [`BASE_INDEX`]
    /// (`A`=0, `C`=1, `G`=2, `T`=3, `N`/other=4).
    base_counts: Vec<[u64; 5]>,
}

impl TailAccumulator {
    /// Map an ASCII base byte to the [`BASES`] slot index. Any non-ACGT byte
    /// (including IUPAC ambiguity codes and lowercase letters that don't
    /// survive case-fold) lands in the `N` slot, so the consensus call never
    /// emits a surprise base.
    fn base_to_slot(b: u8) -> usize {
        match b.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        }
    }

    /// Increments the count and folds `tail`'s base distribution into
    /// `base_counts`, extending the per-position vectors as needed.
    fn observe(&mut self, tail: &[u8]) {
        self.count += 1;
        if self.base_counts.len() < tail.len() {
            self.base_counts.resize(tail.len(), [0; 5]);
        }
        for (i, &b) in tail.iter().enumerate() {
            self.base_counts[i][Self::base_to_slot(b)] += 1;
        }
    }

    /// Merges `other` into `self`, summing counts and base counts element-wise.
    /// Used by [`aggregate_kmers`] when folding a near-identical bucket into a
    /// denser primary.
    fn merge(&mut self, other: &TailAccumulator) {
        self.count += other.count;
        if self.base_counts.len() < other.base_counts.len() {
            self.base_counts.resize(other.base_counts.len(), [0; 5]);
        }
        for (i, oc) in other.base_counts.iter().enumerate() {
            for (s, &v) in oc.iter().enumerate() {
                self.base_counts[i][s] += v;
            }
        }
    }

    /// Builds the consensus sequence by walking positions and stopping at the
    /// first column that fails any of three checks:
    ///   1. column coverage below `min_coverage` (the adaptive floor
    ///      [`aggregate_kmers`] computes from cluster size);
    ///   2. majority fraction below [`CONSENSUS_MAJORITY_FLOOR`] — a
    ///      genuinely ambiguous column (barcode / soft-clip / etc.);
    ///   3. majority fraction below the running-min baseline by more than
    ///      [`CONSENSUS_DROP_TOLERANCE`] — a discontinuity, typically the
    ///      first position of a variable region behind a stable adapter
    ///      prefix (e.g. i7 barcode after a kit-canonical 16 bp).
    ///
    /// The baseline is a running minimum of the majority fractions of all
    /// columns emitted so far, so the check adapts to gently-noisy libraries
    /// (90 → 88 → 86 → 85 → 80 all pass) but catches an abrupt drop
    /// (95 → 40 cuts at the 40).
    ///
    /// Returns uppercase ASCII bytes. The caller ([`aggregate_kmers`]) is
    /// responsible for rejecting results shorter than a caller-supplied
    /// minimum — a consensus that cuts short of [`TAIL_KMER_LEN`] indicates
    /// the cluster itself is heterogeneous within the k-mer window and
    /// should not be reported as an adapter.
    fn consensus(&self, min_coverage: u64) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.base_counts.len());
        let mut baseline: Option<f64> = None;
        for col in &self.base_counts {
            let total: u64 = col.iter().sum();
            if total < min_coverage {
                break;
            }
            let (best, &best_count) = col.iter().enumerate().max_by_key(|&(_, &v)| v).unwrap();
            let majority = best_count as f64 / total as f64;
            if majority < CONSENSUS_MAJORITY_FLOOR {
                break;
            }
            if let Some(b) = baseline
                && majority < b - CONSENSUS_DROP_TOLERANCE
            {
                break;
            }
            out.push(BASES[best]);
            baseline = Some(match baseline {
                Some(b) => b.min(majority),
                None => majority,
            });
        }
        out
    }
}

/// Computes the effective 3' end of a read after applying the user-configured
/// 3'-end cleanups in order: poly-G trim, then poly-X trim (max across A/C/T),
/// then cut-right quality trim. Returns the new effective length — callers
/// reslice as `&seq[..effective_len]` rather than mutating an OwnedRecord.
///
/// The order matters: poly-G first removes 2-color "no signal" tails so the
/// poly-X scan sees the underlying bases; quality trim runs last so it accounts
/// for any bases the homopolymer scans already removed.
fn effective_trimmed_len(
    seq: &[u8],
    qual: &[u8],
    polyg_min_run: usize,
    polyx_min_run: usize,
    quality_trim: Option<QualityTrim>,
) -> usize {
    let mut end = seq.len();
    if polyg_min_run > 0 {
        let tail = find_polyx_tail_len(&seq[..end], b'G');
        if tail >= polyg_min_run {
            end -= tail;
        }
    }
    if polyx_min_run > 0 {
        let tail = [b'A', b'C', b'T']
            .iter()
            .map(|&x| find_polyx_tail_len(&seq[..end], x))
            .max()
            .unwrap_or(0);
        if tail >= polyx_min_run {
            end -= tail;
        }
    }
    if let Some(qt) = quality_trim
        && let Some(cut_at) = cut_right_quality_position(&qual[..end], qt.window, qt.threshold)
    {
        end = end.min(cut_at);
    }
    end
}

/// Increments the bucket keyed by the first [`TAIL_KMER_LEN`] bases of `tail`
/// (case-folded to uppercase) and folds the full tail's base distribution into
/// that bucket's [`TailAccumulator`]. Per-position counts beyond the k-mer key
/// length are populated only by reads whose tail extended that far — which is
/// what drives the consensus length tracking the library's typical adapter
/// readthrough rather than a fixed cap.
///
/// The key is owned because it becomes the HashMap key on a miss; the full
/// tail bytes are passed to `observe` as-is (no pre-uppercase clone) because
/// [`TailAccumulator::base_to_slot`] already case-folds each byte.
fn push_tail_kmer(buckets: &mut HashMap<Vec<u8>, TailAccumulator>, tail: &[u8]) {
    let n = tail.len().min(TAIL_KMER_LEN);
    let key = tail[..n].to_ascii_uppercase();
    buckets.entry(key).or_default().observe(tail);
}

/// Greedy Hamming-distance aggregation of the tail buckets. Buckets are sorted
/// by count (desc), then key length (desc, full 16-bp keys before truncated
/// short-tail keys), then key lexicographically (for determinism). Each bucket
/// folds into the first prior primary whose common prefix matches within
/// `max_hamming`; otherwise it becomes a new primary.
///
/// The returned [`Hit`]s carry the consensus sequence for each primary with the
/// discontinuity-aware column stop from [`TailAccumulator::consensus`]. Hits
/// whose consensus terminates below [`TAIL_KMER_LEN`] are dropped — a cluster
/// whose consensus can't survive the k-mer window's worth of columns is
/// too heterogeneous to be reported as a real adapter.
fn aggregate_kmers(buckets: HashMap<Vec<u8>, TailAccumulator>, max_hamming: usize) -> Vec<Hit> {
    let mut entries: Vec<(Vec<u8>, TailAccumulator)> = buckets.into_iter().collect();
    entries.sort_by(|a, b| {
        b.1.count
            .cmp(&a.1.count)
            .then_with(|| b.0.len().cmp(&a.0.len()))
            .then_with(|| a.0.cmp(&b.0))
    });

    // Cluster: each new bucket either folds into the densest matching primary or
    // becomes a fresh one. Comparisons use the bucket's 16-bp k-mer key, not the
    // full consensus, so clustering decisions don't depend on later-position
    // noise. `primary_count` tracks the primary's OWN count separately from the
    // running merged total, so we can cap the coverage floor at the primary's
    // originating support and prevent very-dense short-key buckets folding into
    // a longer primary from inflating the floor past what the long-tail columns
    // can actually support.
    let mut primaries: Vec<(Vec<u8>, TailAccumulator, u64)> = Vec::new();
    for (kmer, acc) in entries {
        let mut merged = false;
        for (prim_kmer, prim_acc, _) in primaries.iter_mut() {
            let common = prim_kmer.len().min(kmer.len());
            if common == 0 {
                continue;
            }
            let mm =
                count_mismatches_ci_bounded(&prim_kmer[..common], &kmer[..common], max_hamming);
            if mm <= max_hamming {
                prim_acc.merge(&acc);
                merged = true;
                break;
            }
        }
        if !merged {
            let own = acc.count;
            primaries.push((kmer, acc, own));
        }
    }

    primaries.sort_by(|a, b| b.1.count.cmp(&a.1.count).then_with(|| a.0.cmp(&b.0)));
    primaries
        .into_iter()
        .filter_map(|(_kmer, acc, primary_count)| {
            // Adaptive floor: 5% of the ORIGINATING primary's count, never
            // below 5 reads. Using `primary_count` rather than `acc.count`
            // matters when a very-dense short-key bucket folded in: the
            // long-tail positions only have coverage from the original
            // primary, and inflating the floor with the short bucket's
            // count would truncate an otherwise-informative consensus.
            let floor = (primary_count / 20).max(5);
            let seq = acc.consensus(floor);
            if seq.len() < TAIL_KMER_LEN {
                return None;
            }
            Some(Hit { name: None, seq, count: acc.count })
        })
        .collect()
}

/// Keeps only hits whose count is at least `ceil(min_fraction * total)` of the
/// running total of detections (with an implicit floor of 1). Operates on a
/// pre-built `Vec<Hit>` so both the PE and SE paths can share it.
fn filter_min_fraction(hits: &[Hit], total: u64, min_fraction: f64) -> Vec<Hit> {
    if total == 0 {
        return Vec::new();
    }
    let cutoff = (min_fraction * total as f64).ceil() as u64;
    hits.iter().filter(|h| h.count >= cutoff.max(1)).cloned().collect()
}

/// Compares `seq` (uppercased) against every kit adapter's first
/// [`TAIL_KMER_LEN`] bases on the selected mate, returning the closest match
/// within [`KIT_NEAR_MATCH_HAMMING`] mismatches (exact match wins ties).
fn classify_against_kits(seq: &[u8], mate: KitMate) -> Option<KitMatch> {
    let n = seq.len().min(TAIL_KMER_LEN);
    if n == 0 {
        return None;
    }
    let mut best: Option<(usize, &'static str, &'static [u8])> = None;
    for kit in ALL_KITS {
        let kit_seq: &'static [u8] = match mate {
            KitMate::R1 => kit.seq_r1,
            // Many kits have no published R2 adapter (small-RNA). Skip those for R2 lookup.
            KitMate::R2 => match kit.seq_r2 {
                Some(s) => s,
                None => continue,
            },
        };
        let k = n.min(kit_seq.len());
        if k == 0 {
            continue;
        }
        let limit = best.map(|(m, _, _)| m).unwrap_or(KIT_NEAR_MATCH_HAMMING);
        let mm = count_mismatches_ci_bounded(&seq[..k], &kit_seq[..k], limit);
        if mm <= KIT_NEAR_MATCH_HAMMING && best.is_none_or(|(m, _, _)| mm < m) {
            best = Some((mm, kit.name, kit_seq));
        }
    }
    best.map(|(mm, name, kit_seq)| {
        if mm == 0 {
            KitMatch::Exact { kit: name, kit_seq }
        } else {
            KitMatch::Fuzzy { kit: name, kit_seq, mismatches: mm }
        }
    })
}

/// Builds the SE candidate list: every adapter in [`ALL_KITS`] (R1 sequence, and
/// R2 sequence when present, named `<kit>` and `<kit>_r2`) plus any user
/// `--adapter-sequence` (named `user_1`, …) and FASTA records (named after their
/// `>` header). Duplicates are not pruned — distinct names with identical
/// sequences will each get their own counter, which is the right shape for the
/// per-candidate report.
fn build_se_candidates(
    adapter_sequence: &[String],
    adapter_fasta: &Option<PathBuf>,
    adapter_min_length: usize,
) -> Result<Vec<Candidate>> {
    let mut out: Vec<Candidate> = Vec::new();
    for kit in ALL_KITS {
        out.push(Candidate {
            name: kit.name.to_string(),
            adapter: Adapter::new(kit.seq_r1.to_ascii_uppercase()),
        });
        if let Some(s2) = kit.seq_r2 {
            let r2_name = format!("{}_r2", kit.name);
            // Skip if the R2 byte sequence equals R1 (e.g. Nextera) — listing twice
            // would double-count any read that hits this kit on both candidates.
            if s2 != kit.seq_r1 {
                out.push(Candidate {
                    name: r2_name,
                    adapter: Adapter::new(s2.to_ascii_uppercase()),
                });
            }
        }
    }
    for (i, s) in adapter_sequence.iter().enumerate() {
        out.push(Candidate {
            name: format!("user_{}", i + 1),
            adapter: Adapter::new(s.as_bytes().to_ascii_uppercase()),
        });
    }
    if let Some(path) = adapter_fasta {
        let entries = load_adapter_fasta_with_names(path)?;
        for (name, seq) in entries {
            // Reject too-short FASTA candidates for the same reason as
            // --adapter-sequence: they can never satisfy find_adapter_3prime's
            // alignment-length guard and would silently contribute nothing.
            if seq.len() < adapter_min_length {
                return Err(anyhow!(
                    "--adapter-fasta entry {name:?} ({path:?}) is {} bp, shorter than \
                     --adapter-min-length ({adapter_min_length}); it could never match. Lower \
                     --adapter-min-length or edit the FASTA.",
                    seq.len(),
                ));
            }
            out.push(Candidate { name, adapter: Adapter::new(seq.to_ascii_uppercase()) });
        }
    }
    Ok(out)
}

/// Annotates each PE-mate hit with its kit-match decision against the requested
/// mate-side adapter list. Used by both the console report and the FASTA writer.
fn annotate_pe_hits<'a>(hits: &'a [Hit], mate: KitMate) -> Vec<AnnotatedHit<'a>> {
    let label = match mate {
        KitMate::R1 => "R1",
        KitMate::R2 => "R2",
    };
    hits.iter()
        .map(|h| AnnotatedHit {
            hit: h,
            mate: Some(label),
            kit_match: classify_against_kits(&h.seq, mate),
        })
        .collect()
}

/// Annotates each SE hit. SE candidates are already named (kit names, `user_N`,
/// or FASTA header names), but we still run a kit-similarity check on the
/// candidate's bytes so user-supplied sequences that happen to be kit prefixes
/// get the same annotation surface as built-in kits. Checks the R1-side first,
/// then R2-side as a fallback (some SE libraries use the R2 chemistry).
fn annotate_se_hits<'a>(hits: &'a [Hit]) -> Vec<AnnotatedHit<'a>> {
    hits.iter()
        .map(|h| {
            let kit_match = classify_against_kits(&h.seq, KitMate::R1)
                .or_else(|| classify_against_kits(&h.seq, KitMate::R2));
            AnnotatedHit { hit: h, mate: None, kit_match }
        })
        .collect()
}

/// Formats `consensus` with the kit-matched prefix uppercased and any extension
/// past the kit's published length lowercased. Used in the section-2 lines so
/// the reader can see at a glance where the kit-stable region ends and where
/// sample-specific bases (e.g. an i7 barcode tail) begin.
///
/// When `kit_seq` is `None` (no kit match), the entire consensus is uppercased
/// — there's no boundary to mark.
fn format_consensus_with_kit_marker(consensus: &[u8], kit_seq: Option<&[u8]>) -> String {
    let kit_len = kit_seq.map(|k| k.len()).unwrap_or(consensus.len());
    let split = consensus.len().min(kit_len);
    let mut out = String::with_capacity(consensus.len());
    for &b in &consensus[..split] {
        out.push(b.to_ascii_uppercase() as char);
    }
    for &b in &consensus[split..] {
        out.push(b.to_ascii_lowercase() as char);
    }
    out
}

/// Emits the human-readable PE report on stdout via `info!`. Two sections:
///   1. Matched kit(s): one line per mate per kit match, with the kit's full
///      published adapter sequence.
///   2. Full-length discovered consensus: per-mate, the consensus as observed
///      in the data, with uppercase marking the kit-matched prefix and
///      lowercase marking the extension past it (typically the i7 barcode).
fn emit_pe_report(
    detections: u64,
    reads_scanned: u64,
    overlap_hits: u64,
    r1: &[AnnotatedHit<'_>],
    r2: &[AnnotatedHit<'_>],
) {
    info!("chelae detect complete (paired-end):");
    info!(
        "  scanned:    {} pair(s) ({} overlap hits, {} usable detections)",
        fmt_count(reads_scanned),
        fmt_count(overlap_hits),
        fmt_count(detections),
    );
    emit_matched_kits_section(r1.iter().chain(r2.iter()));
    info!("  Full-length discovered consensus (uppercase = matches kit, lowercase = extension):");
    emit_full_length_rows(r1, detections);
    emit_full_length_rows(r2, detections);
}

/// Emits the human-readable SE report on stdout via `info!`. Same two-section
/// layout as PE, with the section-2 sequences rendered through the same
/// uppercase/lowercase formatter (no-op for SE because the candidate's
/// sequence and the matched kit's sequence are typically identical).
fn emit_se_report(detections: u64, reads_scanned: u64, reported: &[AnnotatedHit<'_>]) {
    info!("chelae detect complete (single-end):");
    info!(
        "  scanned:    {} read(s); {} matched a candidate adapter",
        fmt_count(reads_scanned),
        fmt_count(detections),
    );
    emit_matched_kits_section(reported.iter());
    info!("  Candidate adapter(s) above --min-fraction:");
    emit_full_length_rows(reported, detections);
}

/// Section-1 emitter: collects unique `(mate, kit, kit_seq)` triples across all
/// annotated hits and prints each on its own line with the kit's published
/// sequence. Skips the section entirely when no hit matched any kit.
fn emit_matched_kits_section<'a, I>(hits: I)
where
    I: IntoIterator<Item = &'a AnnotatedHit<'a>>,
{
    let mut rows: Vec<(Option<&'static str>, &'static str, &'static [u8], bool)> = Vec::new();
    for ah in hits {
        if let Some(km) = ah.kit_match {
            let exact = matches!(km, KitMatch::Exact { .. });
            let row = (ah.mate, km.kit(), km.kit_seq(), exact);
            if !rows.iter().any(|r| r.0 == row.0 && r.1 == row.1) {
                rows.push(row);
            }
        }
    }
    if rows.is_empty() {
        info!("  Matched kit(s): (none — discovered adapter does not resemble any known kit)");
        return;
    }
    info!("  Matched kit(s):");
    for (mate, kit, kit_seq, exact) in rows {
        let mate_prefix = match mate {
            Some(m) => format!("{m}: "),
            None => String::new(),
        };
        let exactness = if exact { "exact" } else { "fuzzy" };
        info!("    {mate_prefix}{kit} ({exactness}) — {}", String::from_utf8_lossy(kit_seq),);
    }
}

/// Section-2 emitter: per-hit lines showing the full-length consensus (or
/// candidate sequence) with uppercase/lowercase kit-boundary marking, count,
/// fraction, and the kit-match annotation.
fn emit_full_length_rows(hits: &[AnnotatedHit<'_>], total: u64) {
    if hits.is_empty() {
        info!("    (no candidate above the --min-fraction cutoff)");
        return;
    }
    for (i, ah) in hits.iter().enumerate() {
        let pct = (ah.hit.count as f64 * 100.0) / total.max(1) as f64;
        let kit_seq = ah.kit_match.as_ref().map(|km| km.kit_seq());
        let display = format_consensus_with_kit_marker(&ah.hit.seq, kit_seq);
        let kit_note = match ah.kit_match {
            Some(KitMatch::Exact { kit, .. }) => format!("  matches kit {kit}"),
            Some(KitMatch::Fuzzy { kit, mismatches, .. }) => {
                format!("  near kit {kit} ({mismatches} mm in first {TAIL_KMER_LEN} bp)")
            }
            None => "  no kit match".to_string(),
        };
        let mate_label = match (ah.mate, ah.hit.name.as_deref()) {
            (Some(m), _) => format!("{m} [{}]: ", i + 1),
            (None, Some(name)) => format!("[{}] {name}: ", i + 1),
            (None, None) => format!("[{}] ", i + 1),
        };
        info!(
            "    {mate_label}{display}  count={}  ({:.2}%){kit_note}",
            fmt_count(ah.hit.count),
            pct,
        );
    }
}

/// Writes the supplied `(name, sequence)` records to `path` in 2-line FASTA
/// format (`>name\nseq\n`). All FASTA outputs from `detect` route through here
/// — PE uses synthetic names like `r1_adapter` / `r1_adapter_1`, SE uses the
/// candidate's display name.
///
/// Routes through [`Io::new_writer`] so a `.gz`-extensioned path is
/// transparently gzip-compressed — matches trim.rs's writer pattern and
/// keeps the FASTA usable directly with `chelae trim --adapter-fasta` whose
/// reader also handles gzip-by-extension. The first argument to `Io::new`
/// is the gzip compression level; level 5 is the middle ground trim.rs
/// uses for its analogous writers.
fn write_fasta(path: &Path, records: &[(String, &[u8])]) -> Result<()> {
    let mut w = Io::new(5, BUFFER_SIZE)
        .new_writer(path)
        .map_err(|e| anyhow!("Failed to create {path:?}: {e}"))?;
    for (name, seq) in records {
        writeln!(w, ">{name}").map_err(|e| anyhow!("Failed to write {path:?}: {e}"))?;
        w.write_all(seq).map_err(|e| anyhow!("Failed to write {path:?}: {e}"))?;
        writeln!(w).map_err(|e| anyhow!("Failed to write {path:?}: {e}"))?;
    }
    w.flush().map_err(|e| anyhow!("Failed to flush {path:?}: {e}"))
}

/// Builds the PE FASTA record list. Per the design: when a hit is kit-matched
/// (exact OR fuzzy), the FASTA entry uses the **kit's published adapter
/// sequence** rather than the in-sample consensus — the in-sample sequence
/// extends past the kit-stable region into the i7/i5 barcode, which is
/// sample-specific and would break round-trip use of the FASTA on a different
/// sample in the same batch. For novel adapters (no kit match), there's no
/// kit-stable portion to fall back to and the full consensus is emitted.
fn pe_fasta_records<'a>(
    r1: &'a [AnnotatedHit<'a>],
    r2: &'a [AnnotatedHit<'a>],
) -> Vec<(String, &'a [u8])> {
    let mut out: Vec<(String, &[u8])> = Vec::with_capacity(r1.len() + r2.len());
    push_named(&mut out, r1, "r1_adapter");
    push_named(&mut out, r2, "r2_adapter");
    out
}

/// Appends every annotated hit in `hits` to `out` with synthetic names:
/// `{prefix}` for a single entry, `{prefix}_{N}` (1-based) otherwise. Each
/// entry's sequence is the kit's published adapter when a kit-match exists,
/// otherwise the in-sample consensus.
fn push_named<'a>(out: &mut Vec<(String, &'a [u8])>, hits: &'a [AnnotatedHit<'a>], prefix: &str) {
    match hits {
        [] => {}
        [only] => out.push((prefix.to_string(), fasta_seq_for(only))),
        _ => {
            for (i, ah) in hits.iter().enumerate() {
                out.push((format!("{prefix}_{}", i + 1), fasta_seq_for(ah)));
            }
        }
    }
}

/// Returns the sequence to embed in the FASTA for one annotated hit. Kit-matched
/// hits return the kit's full published adapter; novel hits return the
/// consensus.
fn fasta_seq_for<'a>(ah: &'a AnnotatedHit<'a>) -> &'a [u8] {
    match ah.kit_match {
        Some(km) => km.kit_seq(),
        None => &ah.hit.seq,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chelae_lib::adapter_db::{NEXTERA, TRUSEQ};
    use fgoxide::io::Io;
    use std::fs;
    use tempfile::TempDir;

    /// Builds a bucket map from `(sequence, count)` pairs by observing the sequence
    /// `count` times into a fresh [`TailAccumulator`]. The bucket key is the first
    /// `min(seq.len(), TAIL_KMER_LEN)` bases of `seq`, matching the production
    /// keying in [`push_tail_kmer`]. Used by the aggregate_kmers unit tests to
    /// avoid hand-rolling per-position base counts.
    fn bucket_map(entries: &[(&[u8], u64)]) -> HashMap<Vec<u8>, TailAccumulator> {
        let mut out: HashMap<Vec<u8>, TailAccumulator> = HashMap::new();
        for &(seq, count) in entries {
            let n = seq.len().min(TAIL_KMER_LEN);
            let key = seq[..n].to_vec();
            let acc = out.entry(key).or_default();
            for _ in 0..count {
                acc.observe(seq);
            }
        }
        out
    }

    /// Minimal ACGT reverse-complement for synthesizing PE reads in tests. Production
    /// code uses the SIMD `reverse_complement_acgt_into` over in `trim.rs`; tests can
    /// afford the scalar form since they generate at most a few hundred bp.
    fn rc(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            })
            .collect()
    }

    /// Builds a deterministic, non-periodic template of length `n` via a small
    /// xorshift PRBG. Periodic templates (e.g. `ACGTACGT…`) self-align at multiple
    /// shifts and break PE-overlap probing for tests; a PRBG sequence is reproducible
    /// without that pathology.
    fn template_of_len(n: usize) -> Vec<u8> {
        let bases = b"ACGT";
        let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
        (0..n)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                bases[(state as usize) & 0x3]
            })
            .collect()
    }

    /// Builds a paired-end FASTQ record body (no header) for a single pair carrying
    /// the supplied template plus `tail_r1` / `tail_r2` past the template. No tail
    /// padding — the result is `template.len() + tail.len()` long, identical between
    /// pairs so the production reader sees properly synchronized records.
    fn pe_pair(template: &[u8], tail_r1: &[u8], tail_r2: &[u8]) -> (Vec<u8>, Vec<u8>) {
        let mut r1 = template.to_vec();
        r1.extend_from_slice(tail_r1);
        let mut r2 = rc(template);
        r2.extend_from_slice(tail_r2);
        (r1, r2)
    }

    /// Writes a 4-lines-per-record FASTQ at `tmp/{name}.fq`.
    fn write_fq(tmp: &TempDir, name: &str, records: &[(String, Vec<u8>)]) -> PathBuf {
        let path = tmp.path().join(format!("{name}.fq"));
        let mut lines: Vec<String> = Vec::with_capacity(records.len() * 4);
        for (id, seq) in records {
            lines.push(format!("@{id}"));
            lines.push(String::from_utf8(seq.clone()).unwrap());
            lines.push("+".to_string());
            lines.push("I".repeat(seq.len()));
        }
        Io::default().write_lines(&path, &lines).unwrap();
        path
    }

    /// Single test-helper constructor used for both SE (one input) and PE (two
    /// inputs) — the per-mode dispatch happens inside `execute()` based on the
    /// length of `inputs`. Defaults are picked for fast tests (50 detections, 10k
    /// max reads) rather than production realism.
    fn make_detect(inputs: Vec<PathBuf>, output_fasta: Option<PathBuf>) -> Detect {
        Detect {
            inputs,
            output_fasta,
            adapter_sequence: vec![],
            adapter_fasta: None,
            num_detections: 50,
            max_reads: 10_000,
            min_fraction: 0.05,
            min_tail_length: 8,
            min_detections_for_report: 5,
            overlap_min_length: 30,
            overlap_max_mismatch_rate: 0.10,
            overlap_diagnostic_length: 64,
            adapter_min_length: 10,
            adapter_mismatch_rate: 0.125,
            trim_polyg: 10,
            trim_polyx: 5,
            quality_trim: QualityTrimSetting::On(QualityTrim { window: 4, threshold: 20 }),
        }
    }

    #[test]
    fn aggregate_kmers_merges_within_hamming() {
        // Two near-identical 16-mers (1 mismatch at position 0) and an unrelated one.
        let input = bucket_map(&[
            (b"AGATCGGAAGAGCACA".as_slice(), 100),
            (b"TGATCGGAAGAGCACA".as_slice(), 5), // 1 mm from primary
            (b"CTGTCTCTTATACACA".as_slice(), 50),
        ]);
        let out = aggregate_kmers(input, 2);
        assert_eq!(out.len(), 2, "near-identical entry should merge into primary");
        // The primary's position-0 consensus stays A because 100 As beat 5 Ts.
        assert_eq!(out[0].seq, b"AGATCGGAAGAGCACA");
        assert_eq!(out[0].count, 105);
        assert_eq!(out[1].seq, b"CTGTCTCTTATACACA");
        assert_eq!(out[1].count, 50);
    }

    #[test]
    fn aggregate_kmers_keeps_unrelated_entries_distinct() {
        let input = bucket_map(&[
            (b"AGATCGGAAGAGCACA".as_slice(), 100),
            (b"CTGTCTCTTATACACA".as_slice(), 50),
            (b"GGGGGGGGGGGGGGGG".as_slice(), 25),
        ]);
        let out = aggregate_kmers(input, 2);
        assert_eq!(out.len(), 3);
        let counts: Vec<u64> = out.iter().map(|h| h.count).collect();
        assert_eq!(counts, vec![100, 50, 25]);
    }

    #[test]
    fn aggregate_kmers_merges_at_boundary() {
        // Exactly FUZZY_MERGE_HAMMING mismatches should still merge.
        let input = bucket_map(&[
            (b"AAAAAAAAAAAAAAAA".as_slice(), 100),
            (b"CCAAAAAAAAAAAAAA".as_slice(), 10),
        ]);
        let out = aggregate_kmers(input, 2);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].count, 110);
        // After the merge position 0 has A=100, C=10 — consensus picks A.
        assert_eq!(out[0].seq, b"AAAAAAAAAAAAAAAA");
    }

    #[test]
    fn aggregate_kmers_folds_shorter_prefix_into_full_primary() {
        // A shorter tail (truncated by a short adapter overhang) that is a strict
        // prefix of the dense primary should fold into the primary; the primary
        // keeps its (longer) consensus sequence.
        let input = bucket_map(&[
            (b"AGATCGGAAGAGCACA".as_slice(), 100), // 16 bp
            (b"AGATCGGAAGAGC".as_slice(), 7),      // 13 bp prefix of primary
            (b"AGATCGGAAGAGCAC".as_slice(), 4),    // 15 bp prefix of primary
        ]);
        let out = aggregate_kmers(input, 2);
        assert_eq!(out.len(), 1, "shorter prefixes should fold into the longer primary");
        assert_eq!(out[0].seq, b"AGATCGGAAGAGCACA");
        assert_eq!(out[0].count, 111);
    }

    #[test]
    fn aggregate_kmers_keeps_beyond_boundary() {
        // 3 mismatches > FUZZY_MERGE_HAMMING should not merge.
        let input = bucket_map(&[
            (b"AAAAAAAAAAAAAAAA".as_slice(), 100),
            (b"CCCAAAAAAAAAAAAA".as_slice(), 10),
        ]);
        let out = aggregate_kmers(input, 2);
        assert_eq!(out.len(), 2);
    }

    #[test]
    fn aggregate_kmers_extends_consensus_past_kmer_when_supported() {
        // Two reads contribute a 25 bp tail past the primary's 16-bp key — the
        // consensus should run all 25 bp because every position has high coverage.
        let tail_long = b"AGATCGGAAGAGCACACGTCTGAACT";
        let mut input: HashMap<Vec<u8>, TailAccumulator> = HashMap::new();
        let mut acc = TailAccumulator::default();
        for _ in 0..100 {
            acc.observe(tail_long);
        }
        input.insert(tail_long[..16].to_vec(), acc);
        let out = aggregate_kmers(input, 2);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].seq, tail_long);
        assert_eq!(out[0].count, 100);
    }

    #[test]
    fn aggregate_kmers_truncates_consensus_where_coverage_drops() {
        // Bucket where only a small minority of reads have a tail extending
        // past position 16. The consensus should stop where coverage falls
        // below the adaptive floor (max(5, count/20)).
        let mut acc = TailAccumulator::default();
        // 100 reads cover positions 0..16
        for _ in 0..100 {
            acc.observe(b"AGATCGGAAGAGCACA");
        }
        // 3 reads cover positions 0..25 (below the floor of 5)
        for _ in 0..3 {
            acc.observe(b"AGATCGGAAGAGCACATAILBASES");
        }
        let mut input: HashMap<Vec<u8>, TailAccumulator> = HashMap::new();
        input.insert(b"AGATCGGAAGAGCACA".to_vec(), acc);
        let out = aggregate_kmers(input, 2);
        assert_eq!(out.len(), 1);
        // Coverage at positions 16+ is 3, floor is max(5, 103/20) = 5 — truncated.
        assert_eq!(out[0].seq, b"AGATCGGAAGAGCACA");
        assert_eq!(out[0].count, 103);
    }

    #[test]
    fn filter_min_fraction_drops_below_cutoff() {
        let hits = vec![
            Hit { name: None, seq: b"A".to_vec(), count: 600 },
            Hit { name: None, seq: b"B".to_vec(), count: 300 },
            Hit { name: None, seq: b"C".to_vec(), count: 50 }, // 5% exactly
            Hit { name: None, seq: b"D".to_vec(), count: 49 }, // < 5%
        ];
        // total = 999 detections, cutoff = ceil(0.05 * 999) = 50
        let out = filter_min_fraction(&hits, 999, 0.05);
        assert_eq!(out.len(), 3);
        assert_eq!(out[2].seq, b"C");
    }

    #[test]
    fn classify_against_kits_exact_truseq_r1() {
        // First 16 bp of TruSeq R1 adapter, byte-identical.
        let seq = &TRUSEQ.seq_r1[..16];
        let res = classify_against_kits(seq, KitMate::R1).unwrap();
        assert_eq!(res.kit(), "truseq");
        assert!(matches!(res, KitMatch::Exact { .. }));
        // The kit_seq field carries the kit's *full* published adapter sequence.
        assert_eq!(res.kit_seq(), TRUSEQ.seq_r1);
    }

    #[test]
    fn classify_against_kits_fuzzy_truseq_r2() {
        // First 16 bp of TruSeq R2 with one base swapped.
        let r2 = TRUSEQ.seq_r2.unwrap();
        let mut seq = r2[..16].to_vec();
        seq[0] = b'C'; // perturbed from A
        let res = classify_against_kits(&seq, KitMate::R2).unwrap();
        assert_eq!(res.kit(), "truseq");
        assert!(matches!(res, KitMatch::Fuzzy { mismatches: 1, .. }));
        assert_eq!(res.kit_seq(), TRUSEQ.seq_r2.unwrap());
    }

    #[test]
    fn classify_against_kits_no_match_for_random_seq() {
        // 16 bp that doesn't match any kit's first 16 within KIT_NEAR_MATCH_HAMMING.
        let seq = b"GGGGGGGGGGGGGGGG";
        let res = classify_against_kits(seq, KitMate::R1);
        assert_eq!(res, None);
    }

    #[test]
    fn classify_against_kits_nextera_symmetric() {
        // Nextera uses the same adapter on R1 and R2; both lookups should hit.
        let seq = &NEXTERA.seq_r1[..16];
        let r1_match = classify_against_kits(seq, KitMate::R1).unwrap();
        assert_eq!(r1_match.kit(), "nextera");
        assert!(matches!(r1_match, KitMatch::Exact { .. }));
        let r2_match = classify_against_kits(seq, KitMate::R2).unwrap();
        assert_eq!(r2_match.kit(), "nextera");
        assert!(matches!(r2_match, KitMatch::Exact { .. }));
    }

    #[test]
    fn build_se_candidates_includes_every_kit() {
        let cands = build_se_candidates(&[], &None, 10).unwrap();
        let names: Vec<&str> = cands.iter().map(|c| c.name.as_str()).collect();
        // Every kit appears at least once by name.
        for kit in ALL_KITS {
            assert!(names.contains(&kit.name), "kit {} missing from SE candidate list", kit.name);
        }
        // Nextera (R1 == R2) shouldn't appear as both "nextera" and "nextera_r2".
        assert!(!names.contains(&"nextera_r2"));
        // TruSeq (R1 != R2) should have a _r2 entry.
        assert!(names.contains(&"truseq_r2"));
    }

    #[test]
    fn build_se_candidates_adds_user_sequences() {
        let user = vec!["AAAAAAAAAAAA".to_string()];
        let cands = build_se_candidates(&user, &None, 10).unwrap();
        assert!(cands.iter().any(|c| c.name == "user_1"));
    }

    #[test]
    fn push_tail_kmer_uppercases_and_truncates_key_but_keeps_full_tail() {
        let mut h: HashMap<Vec<u8>, TailAccumulator> = HashMap::new();
        push_tail_kmer(&mut h, b"agatcggaagagcacaTAIL"); // 20 bp, mixed case
        let key: Vec<u8> = h.keys().next().cloned().unwrap();
        assert_eq!(key, b"AGATCGGAAGAGCACA");
        let acc = &h[&key];
        assert_eq!(acc.count, 1);
        // The full 20-base tail (uppercased) was observed — positions 16..19 have one
        // count each. So the consensus when coverage permits would extend to 20 bp.
        assert_eq!(acc.base_counts.len(), 20);
        assert_eq!(acc.base_counts[16][TailAccumulator::base_to_slot(b'T')], 1);
        assert_eq!(acc.base_counts[19][TailAccumulator::base_to_slot(b'L')], 1); // L → N slot
    }

    #[test]
    fn push_tail_kmer_short_tail_uses_full_length() {
        let mut h: HashMap<Vec<u8>, TailAccumulator> = HashMap::new();
        push_tail_kmer(&mut h, b"AGAT"); // 4 bp, below TAIL_KMER_LEN
        let acc = &h[b"AGAT".as_ref()];
        assert_eq!(acc.count, 1);
        assert_eq!(acc.base_counts.len(), 4);
    }

    #[test]
    fn pe_end_to_end_identifies_truseq() {
        let tmp = TempDir::new().unwrap();
        let template = template_of_len(80);
        let tail_r1 = &TRUSEQ.seq_r1[..20];
        let tail_r2 = &TRUSEQ.seq_r2.unwrap()[..20];

        let mut r1_recs: Vec<(String, Vec<u8>)> = Vec::new();
        let mut r2_recs: Vec<(String, Vec<u8>)> = Vec::new();
        for i in 0..200 {
            let (r1, r2) = pe_pair(&template, tail_r1, tail_r2);
            r1_recs.push((format!("pair_{i}/1"), r1));
            r2_recs.push((format!("pair_{i}/2"), r2));
        }
        let r1_path = write_fq(&tmp, "r1", &r1_recs);
        let r2_path = write_fq(&tmp, "r2", &r2_recs);
        let fasta_path = tmp.path().join("out.fa");

        make_detect(vec![r1_path, r2_path], Some(fasta_path.clone())).execute().unwrap();

        let fasta = fs::read_to_string(&fasta_path).unwrap();
        // The hit is kit-matched, so the FASTA emits the kit's *published* 33 bp
        // adapter (not the in-sample 20 bp consensus). This is the cross-sample-
        // portable artifact: when the user feeds adapters.fa to another sample
        // in the batch the per-sample i7 barcode tail doesn't poison the match.
        let truseq_r1: &str = std::str::from_utf8(TRUSEQ.seq_r1).unwrap();
        let truseq_r2: &str = std::str::from_utf8(TRUSEQ.seq_r2.unwrap()).unwrap();
        assert!(
            fasta.contains(">r1_adapter\n"),
            "FASTA should label single-winner R1 hit as >r1_adapter; got:\n{fasta}"
        );
        assert!(
            fasta.contains(truseq_r1),
            "FASTA missing kit-published TruSeq R1 {truseq_r1:?}; got:\n{fasta}"
        );
        assert!(
            fasta.contains(">r2_adapter\n"),
            "FASTA should label single-winner R2 hit as >r2_adapter; got:\n{fasta}"
        );
        assert!(
            fasta.contains(truseq_r2),
            "FASTA missing kit-published TruSeq R2 {truseq_r2:?}; got:\n{fasta}"
        );
    }

    #[test]
    fn se_end_to_end_identifies_truseq() {
        let tmp = TempDir::new().unwrap();
        let template = template_of_len(70);
        let tail = &TRUSEQ.seq_r1[..25];

        let mut recs: Vec<(String, Vec<u8>)> = Vec::new();
        for i in 0..200 {
            let mut read = template.clone();
            read.extend_from_slice(tail);
            recs.push((format!("read_{i}"), read));
        }
        let in_path = write_fq(&tmp, "r1", &recs);
        let fasta_path = tmp.path().join("out.fa");

        make_detect(vec![in_path], Some(fasta_path.clone())).execute().unwrap();

        let fasta = fs::read_to_string(&fasta_path).unwrap();
        assert!(
            fasta.contains(">truseq\n"),
            "SE FASTA should label winning candidate as >truseq; got:\n{fasta}"
        );
        let truseq_r1: &str = std::str::from_utf8(TRUSEQ.seq_r1).unwrap();
        assert!(
            fasta.contains(truseq_r1),
            "SE FASTA missing TruSeq R1 full sequence; got:\n{fasta}"
        );
    }

    #[test]
    fn se_does_not_match_unrelated_kits() {
        // Reads that contain only the TruSeq adapter should not register Nextera-shape
        // hits — guards against accidental cross-kit collisions in the SE scorer.
        let tmp = TempDir::new().unwrap();
        let template = template_of_len(70);
        let tail = &TRUSEQ.seq_r1[..25];

        let mut recs: Vec<(String, Vec<u8>)> = Vec::new();
        for i in 0..200 {
            let mut read = template.clone();
            read.extend_from_slice(tail);
            recs.push((format!("read_{i}"), read));
        }
        let in_path = write_fq(&tmp, "r1", &recs);
        let fasta_path = tmp.path().join("out.fa");
        make_detect(vec![in_path], Some(fasta_path.clone())).execute().unwrap();
        let fasta = fs::read_to_string(&fasta_path).unwrap();
        assert!(!fasta.contains(&format!(">{}\n", NEXTERA.name)));
    }

    #[test]
    fn pe_validation_rejects_se_only_flags() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fq(&tmp, "r1", &[("x".to_string(), b"ACGT".to_vec())]);
        let r2 = write_fq(&tmp, "r2", &[("x".to_string(), b"ACGT".to_vec())]);
        let mut cmd = make_detect(vec![r1, r2], None);
        cmd.adapter_sequence = vec!["AAAA".to_string()];
        let err = cmd.execute().unwrap_err().to_string();
        assert!(err.contains("--adapter-sequence is not used in paired-end mode"));
    }

    #[test]
    fn aggregate_kmers_empty_input_returns_empty() {
        let out = aggregate_kmers(HashMap::<Vec<u8>, TailAccumulator>::new(), 2);
        assert!(out.is_empty());
    }

    #[test]
    fn filter_min_fraction_zero_total_returns_empty() {
        // Defensive: a divide-by-zero on the cutoff calculation would be a real bug.
        let hits = vec![Hit { name: None, seq: b"A".to_vec(), count: 0 }];
        assert!(filter_min_fraction(&hits, 0, 0.5).is_empty());
    }

    #[test]
    fn push_tail_kmer_increments_existing_entry() {
        let mut h: HashMap<Vec<u8>, TailAccumulator> = HashMap::new();
        push_tail_kmer(&mut h, b"AGATCGGAAGAGCACA");
        push_tail_kmer(&mut h, b"AGATCGGAAGAGCACA");
        push_tail_kmer(&mut h, b"agatcggaagagcaca"); // case-insensitive
        assert_eq!(h[b"AGATCGGAAGAGCACA".as_ref()].count, 3);
    }

    #[test]
    fn classify_against_kits_short_sequence_still_classifies() {
        // Tails shorter than TAIL_KMER_LEN are common when the adapter overhang is
        // small; we still want to identify the kit from whatever we have.
        let short = &TRUSEQ.seq_r1[..10];
        let res = classify_against_kits(short, KitMate::R1).unwrap();
        assert_eq!(res.kit(), "truseq");
        assert!(matches!(res, KitMatch::Exact { .. }));
    }

    #[test]
    fn classify_against_kits_r2_skips_small_rna() {
        // small-RNA has no R2 adapter; an R2 lookup for a sequence resembling its
        // R1 must not accidentally return "small-rna" from cross-side matching.
        let small_rna_r1_prefix: &[u8] = &chelae_lib::adapter_db::SMALL_RNA.seq_r1[..16];
        let res = classify_against_kits(small_rna_r1_prefix, KitMate::R2);
        let kit_name = res.map(|m| m.kit());
        assert!(
            kit_name != Some("small-rna"),
            "R2 lookup should skip the R1-only small-rna kit; got {res:?}"
        );
    }

    #[test]
    fn build_se_candidates_preserves_fasta_record_names() {
        let tmp = TempDir::new().unwrap();
        let fa = tmp.path().join("extra.fa");
        std::fs::write(&fa, ">foo description here\nACGTACGTACGT\n>bar\nTTTTAAAATTTT\n").unwrap();
        let cands = build_se_candidates(&[], &Some(fa), 10).unwrap();
        let names: Vec<&str> = cands.iter().map(|c| c.name.as_str()).collect();
        assert!(names.contains(&"foo"), "expected FASTA name 'foo' in {names:?}");
        assert!(names.contains(&"bar"), "expected FASTA name 'bar' in {names:?}");
        let foo = cands.iter().find(|c| c.name == "foo").unwrap();
        assert_eq!(foo.adapter.bytes, b"ACGTACGTACGT");
    }

    #[test]
    fn build_se_candidates_fasta_without_header_gets_synthetic_name() {
        let tmp = TempDir::new().unwrap();
        let fa = tmp.path().join("noheader.fa");
        // FASTA with a header-less leading record (malformed but recoverable);
        // the loader should give it `record_1`.
        std::fs::write(&fa, "ACGTACGTACGT\n>named\nTTTTAAAATTTT\n").unwrap();
        let cands = build_se_candidates(&[], &Some(fa), 10).unwrap();
        let names: Vec<&str> = cands.iter().map(|c| c.name.as_str()).collect();
        assert!(names.contains(&"record_1"));
        assert!(names.contains(&"named"));
    }

    #[test]
    fn pe_r2_short_produces_clear_error() {
        // R2 with fewer records than R1 should produce a synchronization error,
        // not silently drop the extra R1s or panic.
        let tmp = TempDir::new().unwrap();
        let r1_recs = vec![
            ("a/1".to_string(), b"ACGTACGTACGT".to_vec()),
            ("b/1".to_string(), b"ACGTACGTACGT".to_vec()),
        ];
        let r2_recs = vec![("a/2".to_string(), b"ACGTACGTACGT".to_vec())];
        let r1_path = write_fq(&tmp, "r1", &r1_recs);
        let r2_path = write_fq(&tmp, "r2", &r2_recs);
        let err = make_detect(vec![r1_path, r2_path], None).execute().unwrap_err().to_string();
        assert!(err.contains("R2 exhausted before R1"), "expected out-of-sync error; got: {err}");
    }

    #[test]
    fn pe_r1_short_produces_clear_error() {
        // Symmetric to the R2-short case: when R1 hits EOF but R2 still has
        // records, we should report it rather than silently truncating.
        let tmp = TempDir::new().unwrap();
        let r1_recs = vec![("a/1".to_string(), b"ACGTACGTACGT".to_vec())];
        let r2_recs = vec![
            ("a/2".to_string(), b"ACGTACGTACGT".to_vec()),
            ("b/2".to_string(), b"ACGTACGTACGT".to_vec()),
        ];
        let r1_path = write_fq(&tmp, "r1", &r1_recs);
        let r2_path = write_fq(&tmp, "r2", &r2_recs);
        let err = make_detect(vec![r1_path, r2_path], None).execute().unwrap_err().to_string();
        assert!(err.contains("R1 exhausted before R2"), "expected out-of-sync error; got: {err}");
    }

    #[test]
    fn pe_fasta_records_uses_unsuffixed_name_for_single_winner() {
        // Single-winner naming is part of the user-facing contract — a `>r1_adapter`
        // entry feeds back into `chelae trim` more naturally than `>r1_adapter_1`.
        let r1 = [Hit { name: None, seq: b"AAAA".to_vec(), count: 10 }];
        let r2 = [Hit { name: None, seq: b"CCCC".to_vec(), count: 10 }];
        let r1_annot: Vec<AnnotatedHit<'_>> =
            r1.iter().map(|h| AnnotatedHit { hit: h, mate: Some("R1"), kit_match: None }).collect();
        let r2_annot: Vec<AnnotatedHit<'_>> =
            r2.iter().map(|h| AnnotatedHit { hit: h, mate: Some("R2"), kit_match: None }).collect();
        let recs = pe_fasta_records(&r1_annot, &r2_annot);
        let names: Vec<&str> = recs.iter().map(|(n, _)| n.as_str()).collect();
        assert_eq!(names, ["r1_adapter", "r2_adapter"]);
    }

    #[test]
    fn pe_fasta_records_suffixes_multi_winner_names() {
        let r1 = [
            Hit { name: None, seq: b"AAAA".to_vec(), count: 10 },
            Hit { name: None, seq: b"GGGG".to_vec(), count: 5 },
        ];
        let r1_annot: Vec<AnnotatedHit<'_>> =
            r1.iter().map(|h| AnnotatedHit { hit: h, mate: Some("R1"), kit_match: None }).collect();
        let recs = pe_fasta_records(&r1_annot, &[]);
        let names: Vec<&str> = recs.iter().map(|(n, _)| n.as_str()).collect();
        assert_eq!(names, ["r1_adapter_1", "r1_adapter_2"]);
    }

    #[test]
    fn pe_fasta_records_emits_kit_seq_when_kit_matched() {
        // Consensus is "kit + extra bases" (simulating barcode readthrough); the
        // FASTA must emit the kit's canonical published sequence, not the consensus.
        let consensus = {
            let mut v = TRUSEQ.seq_r1.to_vec();
            v.extend_from_slice(b"GGGGGG"); // 6 bp of "barcode"
            v
        };
        let r1 = [Hit { name: None, seq: consensus.clone(), count: 100 }];
        let r1_annot: Vec<AnnotatedHit<'_>> = r1
            .iter()
            .map(|h| AnnotatedHit {
                hit: h,
                mate: Some("R1"),
                kit_match: Some(KitMatch::Exact { kit: "truseq", kit_seq: TRUSEQ.seq_r1 }),
            })
            .collect();
        let recs = pe_fasta_records(&r1_annot, &[]);
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].0, "r1_adapter");
        assert_eq!(recs[0].1, TRUSEQ.seq_r1, "kit-matched FASTA should emit kit's published seq");
        // Crucially: the 6 bp "barcode" extension is NOT in the FASTA output.
        assert!(
            !recs[0].1.windows(6).any(|w| w == b"GGGGGG"),
            "FASTA must not contain the in-sample barcode extension"
        );
    }

    #[test]
    fn pe_fasta_records_emits_consensus_when_no_kit_match() {
        // Novel adapter: no kit_seq to fall back to, so the FASTA contains the
        // full discovered consensus.
        let novel = b"GCGCGCGCGCGCGCGCGCGCGCGCGC".to_vec();
        let r1 = [Hit { name: None, seq: novel.clone(), count: 100 }];
        let r1_annot: Vec<AnnotatedHit<'_>> =
            r1.iter().map(|h| AnnotatedHit { hit: h, mate: Some("R1"), kit_match: None }).collect();
        let recs = pe_fasta_records(&r1_annot, &[]);
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].1, novel);
    }

    #[test]
    fn format_consensus_marks_kit_prefix_uppercase_extension_lowercase() {
        let consensus = b"AGATCGGAAGAGCACATAILBASES";
        let kit_seq = b"AGATCGGAAGAGCACA"; // 16 bp
        let out = format_consensus_with_kit_marker(consensus, Some(kit_seq));
        assert_eq!(out, "AGATCGGAAGAGCACAtailbases");
    }

    #[test]
    fn format_consensus_all_upper_when_no_kit() {
        let consensus = b"AGATCGGAAGAGCACA";
        let out = format_consensus_with_kit_marker(consensus, None);
        assert_eq!(out, "AGATCGGAAGAGCACA");
    }

    #[test]
    fn format_consensus_handles_consensus_shorter_than_kit() {
        // Consensus didn't extend the full kit length — uppercase everything we
        // have (the entire string is within the kit-stable region).
        let consensus = b"AGATCGGAA";
        let kit_seq = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        let out = format_consensus_with_kit_marker(consensus, Some(kit_seq));
        assert_eq!(out, "AGATCGGAA");
    }

    #[test]
    fn pe_validation_reports_no_detections_clearly() {
        // A library with every insert larger than the read length yields zero adapter
        // readthrough, so the PE path should bail with a specific error rather than
        // returning empty results silently.
        let tmp = TempDir::new().unwrap();
        let template = template_of_len(150);
        let mut r1_recs: Vec<(String, Vec<u8>)> = Vec::new();
        let mut r2_recs: Vec<(String, Vec<u8>)> = Vec::new();
        for i in 0..50 {
            // Insert > read length: just take 100 bp of the template on each mate,
            // no adapter tail at all.
            let r1 = template[..100].to_vec();
            let r2 = rc(&template[50..150]);
            r1_recs.push((format!("pair_{i}/1"), r1));
            r2_recs.push((format!("pair_{i}/2"), r2));
        }
        let r1_path = write_fq(&tmp, "r1", &r1_recs);
        let r2_path = write_fq(&tmp, "r2", &r2_recs);
        let err = make_detect(vec![r1_path, r2_path], None).execute().unwrap_err().to_string();
        assert!(
            err.contains("No usable PE adapter detections"),
            "expected zero-detection error message; got: {err}"
        );
    }

    // ---- validation tests (one per documented constraint) ----

    /// Returns a baseline valid `Detect` configured for one (SE) input that
    /// validate() accepts; tests perturb a single field and check the right error
    /// surfaces. The file content doesn't matter for `validate()` — only that the
    /// path exists.
    fn valid_se_baseline(tmp: &TempDir) -> Detect {
        let p = tmp.path().join("ok.fq");
        std::fs::write(&p, b"@x\nACGT\n+\nIIII\n").unwrap();
        make_detect(vec![p], None)
    }

    fn assert_validate_err_contains(cmd: &Detect, substr: &str) {
        let err = cmd.validate().unwrap_err().to_string();
        assert!(err.contains(substr), "expected error containing {substr:?}; got: {err}");
    }

    #[test]
    fn validate_rejects_missing_input() {
        let tmp = TempDir::new().unwrap();
        let mut cmd = valid_se_baseline(&tmp);
        cmd.inputs = vec![tmp.path().join("nope.fq")];
        assert_validate_err_contains(&cmd, "does not exist");
    }

    #[test]
    fn validate_rejects_missing_output_parent() {
        let tmp = TempDir::new().unwrap();
        let mut cmd = valid_se_baseline(&tmp);
        cmd.output_fasta = Some(tmp.path().join("no_such_dir/out.fa"));
        assert_validate_err_contains(&cmd, "parent directory");
    }

    #[test]
    fn validate_rejects_empty_adapter_sequence() {
        let tmp = TempDir::new().unwrap();
        let mut cmd = valid_se_baseline(&tmp);
        cmd.adapter_sequence = vec!["".to_string()];
        assert_validate_err_contains(&cmd, "must not be empty");
    }

    #[test]
    fn validate_rejects_invalid_adapter_base() {
        let tmp = TempDir::new().unwrap();
        let mut cmd = valid_se_baseline(&tmp);
        cmd.adapter_sequence = vec!["ACGZ".to_string()];
        assert_validate_err_contains(&cmd, "invalid base");
    }

    #[test]
    fn validate_rejects_out_of_range_min_fraction() {
        let tmp = TempDir::new().unwrap();
        let mut cmd = valid_se_baseline(&tmp);
        cmd.min_fraction = -0.1;
        assert_validate_err_contains(&cmd, "--min-fraction");
    }

    #[test]
    fn validate_rejects_out_of_range_overlap_max_mismatch_rate() {
        let tmp = TempDir::new().unwrap();
        let mut cmd = valid_se_baseline(&tmp);
        cmd.overlap_max_mismatch_rate = 1.5;
        assert_validate_err_contains(&cmd, "--overlap-max-mismatch-rate");
    }

    #[test]
    fn validate_rejects_out_of_range_adapter_mismatch_rate() {
        let tmp = TempDir::new().unwrap();
        let mut cmd = valid_se_baseline(&tmp);
        cmd.adapter_mismatch_rate = 2.0;
        assert_validate_err_contains(&cmd, "--adapter-mismatch-rate");
    }

    #[test]
    fn validate_rejects_zero_num_detections() {
        let tmp = TempDir::new().unwrap();
        let mut cmd = valid_se_baseline(&tmp);
        cmd.num_detections = 0;
        assert_validate_err_contains(&cmd, "--num-detections");
    }

    #[test]
    fn validate_rejects_zero_max_reads() {
        let tmp = TempDir::new().unwrap();
        let mut cmd = valid_se_baseline(&tmp);
        cmd.max_reads = 0;
        assert_validate_err_contains(&cmd, "--max-reads");
    }

    #[test]
    fn validate_rejects_zero_min_detections_for_report() {
        let tmp = TempDir::new().unwrap();
        let mut cmd = valid_se_baseline(&tmp);
        cmd.min_detections_for_report = 0;
        assert_validate_err_contains(&cmd, "--min-detections-for-report must be at least 1");
    }

    #[test]
    fn validate_rejects_floor_above_target_detections() {
        let tmp = TempDir::new().unwrap();
        let mut cmd = valid_se_baseline(&tmp);
        cmd.min_detections_for_report = 100;
        cmd.num_detections = 50;
        assert_validate_err_contains(&cmd, "cannot exceed --num-detections");
    }

    #[test]
    fn validate_rejects_floor_above_max_reads() {
        let tmp = TempDir::new().unwrap();
        let mut cmd = valid_se_baseline(&tmp);
        cmd.num_detections = 200_000; // big enough to pass floor-vs-target check
        cmd.max_reads = 100;
        cmd.min_detections_for_report = 1000;
        assert_validate_err_contains(&cmd, "cannot exceed --max-reads");
    }

    #[test]
    fn validate_rejects_zero_min_tail_length() {
        let tmp = TempDir::new().unwrap();
        let mut cmd = valid_se_baseline(&tmp);
        cmd.min_tail_length = 0;
        assert_validate_err_contains(&cmd, "--min-tail-length");
    }

    #[test]
    fn validate_rejects_zero_overlap_min_length() {
        let tmp = TempDir::new().unwrap();
        let mut cmd = valid_se_baseline(&tmp);
        cmd.overlap_min_length = 0;
        assert_validate_err_contains(&cmd, "--overlap-min-length");
    }

    #[test]
    fn validate_rejects_zero_overlap_diagnostic_length() {
        let tmp = TempDir::new().unwrap();
        let mut cmd = valid_se_baseline(&tmp);
        cmd.overlap_diagnostic_length = 0;
        assert_validate_err_contains(&cmd, "--overlap-diagnostic-length");
    }

    #[test]
    fn validate_rejects_zero_adapter_min_length() {
        let tmp = TempDir::new().unwrap();
        let mut cmd = valid_se_baseline(&tmp);
        cmd.adapter_min_length = 0;
        assert_validate_err_contains(&cmd, "--adapter-min-length");
    }

    #[test]
    fn validate_pe_rejects_adapter_fasta() {
        let tmp = TempDir::new().unwrap();
        let r1 = tmp.path().join("r1.fq");
        let r2 = tmp.path().join("r2.fq");
        std::fs::write(&r1, b"@x\nACGT\n+\nIIII\n").unwrap();
        std::fs::write(&r2, b"@x\nACGT\n+\nIIII\n").unwrap();
        let mut cmd = make_detect(vec![r1, r2], None);
        cmd.adapter_fasta = Some(tmp.path().join("x.fa"));
        assert_validate_err_contains(&cmd, "--adapter-fasta is not used in paired-end mode");
    }

    // ---- min-fraction boundary ----

    #[test]
    fn filter_min_fraction_passes_singleton_at_zero_fraction() {
        // The .max(1) floor means even --min-fraction 0.0 still requires count >= 1,
        // i.e. you always get every hit that was actually seen.
        let hits = vec![Hit { name: None, seq: b"A".to_vec(), count: 1 }];
        let out = filter_min_fraction(&hits, 1000, 0.0);
        assert_eq!(out.len(), 1);
    }

    #[test]
    fn filter_min_fraction_at_one_keeps_only_full_winners() {
        // --min-fraction 1.0 -> cutoff = total; only hits with count == total pass.
        let hits = vec![
            Hit { name: None, seq: b"A".to_vec(), count: 50 },
            Hit { name: None, seq: b"B".to_vec(), count: 50 },
        ];
        // Neither alone makes 100% of total=100; both should be dropped.
        let out = filter_min_fraction(&hits, 100, 1.0);
        assert!(out.is_empty());
    }

    // ---- --output-fasta runtime errors ----

    #[test]
    fn se_missing_adapter_fasta_errors_cleanly() {
        let tmp = TempDir::new().unwrap();
        let in_path = write_fq(
            &tmp,
            "r1",
            &[("x".to_string(), b"ACGT".repeat(25))].into_iter().collect::<Vec<_>>(),
        );
        let mut cmd = make_detect(vec![in_path], None);
        cmd.adapter_fasta = Some(tmp.path().join("does_not_exist.fa"));
        let err = cmd.execute().unwrap_err().to_string();
        assert!(
            err.contains("Failed to open adapter FASTA") || err.contains("does_not_exist"),
            "expected FASTA-open error; got: {err}"
        );
    }

    #[test]
    fn se_empty_adapter_fasta_loads_to_no_extra_candidates() {
        let tmp = TempDir::new().unwrap();
        let fa = tmp.path().join("empty.fa");
        std::fs::write(&fa, b"").unwrap();
        let cands = build_se_candidates(&[], &Some(fa), 10).unwrap();
        // Should be exactly the kit candidates — no extras from the empty FASTA.
        let baseline = build_se_candidates(&[], &None, 10).unwrap();
        assert_eq!(cands.len(), baseline.len());
    }

    // ---- detection floor ----

    #[test]
    fn pe_below_detection_floor_errors_instead_of_reporting() {
        // Build a library where only a few pairs have detectable readthrough — fewer
        // than the floor. The PE path should refuse to produce a report rather than
        // emit confident-looking percentages from a tiny sample.
        let tmp = TempDir::new().unwrap();
        let mut r1_recs: Vec<(String, Vec<u8>)> = Vec::new();
        let mut r2_recs: Vec<(String, Vec<u8>)> = Vec::new();
        let tail_r1 = &TRUSEQ.seq_r1[..20];
        let tail_r2 = &TRUSEQ.seq_r2.unwrap()[..20];
        for i in 0..3 {
            let template = template_of_len(80);
            let (r1, r2) = pe_pair(&template, tail_r1, tail_r2);
            r1_recs.push((format!("p_{i}/1"), r1));
            r2_recs.push((format!("p_{i}/2"), r2));
        }
        let r1_path = write_fq(&tmp, "r1", &r1_recs);
        let r2_path = write_fq(&tmp, "r2", &r2_recs);
        let mut cmd = make_detect(vec![r1_path, r2_path], None);
        cmd.min_detections_for_report = 10;
        let err = cmd.execute().unwrap_err().to_string();
        assert!(
            err.contains("below the `--min-detections-for-report` floor"),
            "expected floor-violation error; got: {err}"
        );
    }

    // ---- varied-insert PE smoke test (exercises OverlapStats center update) ----

    #[test]
    fn pe_varied_inserts_still_identify_truseq() {
        // Each pair gets its own deterministic template AND insert size jitter, so
        // OverlapStats.expected_insert actually drifts during the run. If the
        // center-update logic regressed, this test would catch the resulting
        // detection drop.
        let tmp = TempDir::new().unwrap();
        let tail_r1 = &TRUSEQ.seq_r1[..20];
        let tail_r2 = &TRUSEQ.seq_r2.unwrap()[..20];
        let mut r1_recs: Vec<(String, Vec<u8>)> = Vec::new();
        let mut r2_recs: Vec<(String, Vec<u8>)> = Vec::new();
        // Vary insert from 60 to 90 bp across 200 pairs.
        for i in 0..200u64 {
            let ins = 60 + (i.wrapping_mul(31) % 31) as usize; // 60..=90
            let template = template_of_len_seeded(ins, 0xABCD ^ i);
            let (r1, r2) = pe_pair(&template, tail_r1, tail_r2);
            r1_recs.push((format!("p_{i}/1"), r1));
            r2_recs.push((format!("p_{i}/2"), r2));
        }
        let r1_path = write_fq(&tmp, "r1", &r1_recs);
        let r2_path = write_fq(&tmp, "r2", &r2_recs);
        let fasta_path = tmp.path().join("out.fa");
        make_detect(vec![r1_path, r2_path], Some(fasta_path.clone())).execute().unwrap();
        let fasta = std::fs::read_to_string(&fasta_path).unwrap();
        let truseq_r1: &str = std::str::from_utf8(&TRUSEQ.seq_r1[..TAIL_KMER_LEN]).unwrap();
        assert!(
            fasta.contains(truseq_r1),
            "varied-insert PE should still identify TruSeq R1; got:\n{fasta}"
        );
    }

    /// Variant of `template_of_len` that takes a seed so callers can vary the
    /// template per pair.
    fn template_of_len_seeded(n: usize, seed: u64) -> Vec<u8> {
        let bases = b"ACGT";
        let mut state: u64 = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        (0..n)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                bases[(state as usize) & 0x3]
            })
            .collect()
    }

    // ---- 3'-end trim integration ----

    #[test]
    fn effective_trimmed_len_strips_polyg_tail() {
        // 20 bp of normal bases + 15 G's; default polyg_min_run=10 should remove all 15.
        let seq = b"AGATCGGAAGAGCACAGTGTGGGGGGGGGGGGGGG";
        let qual = vec![b'I'; seq.len()];
        let end = effective_trimmed_len(seq, &qual, 10, 0, None);
        assert_eq!(end, 20);
    }

    #[test]
    fn effective_trimmed_len_strips_polyx_a_tail() {
        let seq = b"AGATCGGAAGAGCACAGTGTAAAAAAAAAAAAAAA"; // 20 + 15 A's
        let qual = vec![b'I'; seq.len()];
        // poly-X with min_run=5 should strip the 15-base A tail.
        let end = effective_trimmed_len(seq, &qual, 0, 5, None);
        assert_eq!(end, 20);
    }

    #[test]
    fn effective_trimmed_len_polyg_runs_before_polyx() {
        // GGGGG (5) + AAAAA (5) at the 3' end. polyg trims the G-run, exposing the
        // poly-A which the polyx pass then trims. Without the ordered application
        // the polyx pass would see G's at the tail and find no A-run to trim.
        let seq = b"AGATCGGAAGAGCACAGTGTAAAAAGGGGG"; // 20 + 5 A + 5 G
        let qual = vec![b'I'; seq.len()];
        let end = effective_trimmed_len(seq, &qual, 5, 5, None);
        assert_eq!(end, 20);
    }

    #[test]
    fn effective_trimmed_len_quality_cut_after_homopolymers() {
        // 20 high-Q bases + 10 low-Q bases. Default 4:20 quality trim with cut-right
        // semantics cuts at the start of the first failing window: window starting
        // at position 19 covers qual[19..23] = high,low,low,low → mean 10 < 20 →
        // fail. The cut therefore lands at 19 (keeps 19 bases), not 20 — the last
        // high-Q base is "shared" with the failing window and dropped with it.
        let seq = b"AGATCGGAAGAGCACAGTGTACGTACGTAC"; // 30 bp
        let mut qual = vec![b'I'; 20];
        qual.extend(vec![b'!'; 10]); // Phred 0 — fails any reasonable threshold
        let end =
            effective_trimmed_len(seq, &qual, 0, 0, Some(QualityTrim { window: 4, threshold: 20 }));
        assert_eq!(end, 19);
    }

    #[test]
    fn effective_trimmed_len_all_disabled_is_noop() {
        let seq = b"AGATCGGAAGAGCACAGTGT";
        let qual = vec![b'I'; seq.len()];
        assert_eq!(effective_trimmed_len(seq, &qual, 0, 0, None), seq.len());
    }

    #[test]
    fn effective_trimmed_len_min_run_gate_respected() {
        // 3 trailing G's, below the default min_run of 10. No trim.
        let seq = b"AGATCGGAAGAGCACAGTGTGGG";
        let qual = vec![b'I'; seq.len()];
        let end = effective_trimmed_len(seq, &qual, 10, 0, None);
        assert_eq!(end, seq.len());
    }

    #[test]
    fn quality_trim_setting_accepts_off_aliases() {
        // The CLI accepts "off"/"none"/"no" (case-insensitive) to disable quality trim.
        assert!(QualityTrimSetting::from_str("off").unwrap().as_option().is_none());
        assert!(QualityTrimSetting::from_str("OFF").unwrap().as_option().is_none());
        assert!(QualityTrimSetting::from_str("None").unwrap().as_option().is_none());
        assert!(QualityTrimSetting::from_str("no").unwrap().as_option().is_none());
        // And a real W:Q parses through.
        let qt = QualityTrimSetting::from_str("8:25").unwrap().as_option().unwrap();
        assert_eq!(qt.window, 8);
        assert_eq!(qt.threshold, 25);
        // Garbage rejected.
        assert!(QualityTrimSetting::from_str("garbage").is_err());
    }

    // ---- Consensus discontinuity cut ----

    /// Helper: build a `TailAccumulator` for `n_positions` where each position has
    /// `majority` of `base` and the rest distributed evenly across the other three
    /// ACGT slots. Coverage per column is `total`. Used by the discontinuity tests
    /// to construct controlled majority curves.
    fn tail_accumulator_with_majorities(
        majorities: &[f64],
        base: u8,
        total: u64,
    ) -> TailAccumulator {
        let base_slot = TailAccumulator::base_to_slot(base);
        let mut acc = TailAccumulator { count: total, base_counts: Vec::new() };
        for &frac in majorities {
            let top = (total as f64 * frac).round() as u64;
            let other = (total - top) / 3;
            let mut col = [0u64; 5];
            col[base_slot] = top;
            for (slot, c) in col.iter_mut().enumerate().take(4) {
                if slot != base_slot {
                    *c = other;
                }
            }
            // Any accounting slack from integer rounding lands in N so column sum stays == total.
            let assigned: u64 = col.iter().sum();
            col[4] += total - assigned;
            acc.base_counts.push(col);
        }
        acc
    }

    #[test]
    fn consensus_extends_through_steady_high_majority() {
        // 20 columns all at ~95% majority: consensus should emit all 20.
        let acc = tail_accumulator_with_majorities(&[0.95; 20], b'A', 100);
        let seq = acc.consensus(5);
        assert_eq!(seq.len(), 20);
        assert!(seq.iter().all(|&b| b == b'A'));
    }

    #[test]
    fn consensus_extends_through_gently_declining_majority() {
        // 96, 94, 92, 90, 88, 86, 84, 82, 80 — each drop is 2 pp so the running-min
        // baseline drifts down smoothly and no column exceeds CONSENSUS_DROP_TOLERANCE.
        let majorities = [0.96, 0.94, 0.92, 0.90, 0.88, 0.86, 0.84, 0.82, 0.80];
        let acc = tail_accumulator_with_majorities(&majorities, b'A', 100);
        let seq = acc.consensus(5);
        assert_eq!(seq.len(), majorities.len());
    }

    #[test]
    fn consensus_cuts_at_sharp_discontinuity() {
        // 16 stable positions at 95%, then a sudden drop to 82% (>10 pp below baseline).
        let mut majorities = vec![0.95; 16];
        majorities.push(0.82);
        majorities.push(0.85);
        let acc = tail_accumulator_with_majorities(&majorities, b'A', 100);
        let seq = acc.consensus(5);
        assert_eq!(seq.len(), 16, "should cut at the 0.82 discontinuity");
    }

    #[test]
    fn consensus_cuts_at_absolute_floor_when_pool_dominates() {
        // Variable-region column: user's imbalanced 10-index pool where one index
        // dominates at 40%. 40% is below CONSENSUS_MAJORITY_FLOOR (50%), so we
        // cut regardless of baseline drift.
        let mut majorities = vec![0.95; 16];
        majorities.push(0.40);
        let acc = tail_accumulator_with_majorities(&majorities, b'A', 100);
        let seq = acc.consensus(5);
        assert_eq!(seq.len(), 16);
    }

    #[test]
    fn consensus_stops_when_coverage_below_min() {
        // Coverage 100 for positions 0-9, then drops to 4 (below min_coverage of 5).
        let mut acc = tail_accumulator_with_majorities(&[0.95; 10], b'A', 100);
        acc.base_counts.push([4, 0, 0, 0, 0]); // A=4, total=4 < 5
        let seq = acc.consensus(5);
        assert_eq!(seq.len(), 10);
    }

    #[test]
    fn aggregate_kmers_drops_primaries_with_short_consensus() {
        // A cluster whose 16-bp k-mer prefix itself is heterogeneous (position 3
        // splits ~50/50) — consensus falls below CONSENSUS_MAJORITY_FLOOR before
        // reaching TAIL_KMER_LEN. aggregate_kmers should drop it.
        let mut input: HashMap<Vec<u8>, TailAccumulator> = HashMap::new();
        // Two 16-bp buckets identical except at position 3 (T vs A), Hamming 1 => merge.
        let mut acc_a = TailAccumulator::default();
        let mut acc_b = TailAccumulator::default();
        for _ in 0..50 {
            acc_a.observe(b"ACGTACGTACGTACGT");
            acc_b.observe(b"ACGAACGTACGTACGT");
        }
        input.insert(b"ACGTACGTACGTACGT".to_vec(), acc_a);
        input.insert(b"ACGAACGTACGTACGT".to_vec(), acc_b);
        let out = aggregate_kmers(input, 2);
        // The merged cluster's position 3 is 50/50 => below the majority floor at
        // position 3 (<16), so the primary is dropped.
        assert!(out.is_empty(), "merged 50/50-at-pos-3 cluster should be dropped; got {out:?}");
    }

    // ---- Validate() boundary equality tests ----

    #[test]
    fn validate_accepts_floor_equal_to_num_detections() {
        let tmp = TempDir::new().unwrap();
        let fq = write_fq(&tmp, "r", &[("x".to_string(), b"ACGT".to_vec())]);
        let mut cmd = make_detect(vec![fq], None);
        cmd.num_detections = 50;
        cmd.min_detections_for_report = 50; // equality — sampler can just reach the floor
        assert!(cmd.validate().is_ok());
    }

    #[test]
    fn validate_accepts_floor_equal_to_max_reads() {
        let tmp = TempDir::new().unwrap();
        let fq = write_fq(&tmp, "r", &[("x".to_string(), b"ACGT".to_vec())]);
        let mut cmd = make_detect(vec![fq], None);
        cmd.max_reads = 50;
        cmd.num_detections = 50;
        cmd.min_detections_for_report = 50;
        assert!(cmd.validate().is_ok());
    }

    #[test]
    fn validate_rejects_user_candidate_shorter_than_adapter_min_length() {
        let tmp = TempDir::new().unwrap();
        let fq = write_fq(&tmp, "r", &[("x".to_string(), b"ACGTACGTAC".to_vec())]);
        let mut cmd = make_detect(vec![fq], None);
        cmd.adapter_sequence = vec!["ACGTAC".to_string()]; // 6 bp < adapter_min_length=10
        let err = cmd.validate().unwrap_err().to_string();
        assert!(
            err.contains("shorter than --adapter-min-length"),
            "expected too-short error; got: {err}"
        );
    }

    // ---- SE min-detections-for-report floor ----

    #[test]
    fn se_below_detection_floor_errors_instead_of_reporting() {
        // Mirror of the PE floor test: only a handful of SE reads carry a real
        // TruSeq tail; scanning stops well below the floor.
        let tmp = TempDir::new().unwrap();
        let tail = &TRUSEQ.seq_r1[..25];
        let mut recs: Vec<(String, Vec<u8>)> = Vec::new();
        for i in 0..3 {
            let mut read = template_of_len(70);
            read.extend_from_slice(tail);
            recs.push((format!("r_{i}"), read));
        }
        let path = write_fq(&tmp, "r1", &recs);
        let mut cmd = make_detect(vec![path], None);
        cmd.min_detections_for_report = 10;
        let err = cmd.execute().unwrap_err().to_string();
        assert!(
            err.contains("below the `--min-detections-for-report` floor"),
            "expected SE floor-violation error; got: {err}"
        );
    }

    // ---- Hard-fail on incomplete/empty FASTA output ----

    #[test]
    fn pe_output_fasta_errors_when_min_fraction_excludes_all_hits() {
        // Split the library across three distinct adapter tails so no cluster
        // holds a large fraction of detections; a tight --min-fraction then
        // excludes them all. --output-fasta should hard-fail rather than write
        // an empty file that downstream trim would silently accept.
        let tmp = TempDir::new().unwrap();
        // Slice to a length shorter than the shortest kit adapter (Nextera == 19 bp).
        let cap = 18;
        let truseq_r1 = &TRUSEQ.seq_r1[..cap];
        let truseq_r2 = &TRUSEQ.seq_r2.unwrap()[..cap];
        let nextera_r1 = &NEXTERA.seq_r1[..cap];
        let nextera_r2 = &NEXTERA.seq_r2.unwrap()[..cap];
        let novel_r1: &[u8] = b"CCCCGGGGAAAATTTTAAAA";
        let novel_r2: &[u8] = b"TTTTAAAAGGGGCCCCGGGG";
        let mut r1_recs: Vec<(String, Vec<u8>)> = Vec::new();
        let mut r2_recs: Vec<(String, Vec<u8>)> = Vec::new();
        for i in 0..30 {
            let template = template_of_len(80);
            let (r1, r2) = pe_pair(&template, truseq_r1, truseq_r2);
            r1_recs.push((format!("a_{i}/1"), r1));
            r2_recs.push((format!("a_{i}/2"), r2));
            let template = template_of_len(80);
            let (r1, r2) = pe_pair(&template, nextera_r1, nextera_r2);
            r1_recs.push((format!("b_{i}/1"), r1));
            r2_recs.push((format!("b_{i}/2"), r2));
            let template = template_of_len(80);
            let (r1, r2) = pe_pair(&template, novel_r1, novel_r2);
            r1_recs.push((format!("c_{i}/1"), r1));
            r2_recs.push((format!("c_{i}/2"), r2));
        }
        let r1_path = write_fq(&tmp, "r1", &r1_recs);
        let r2_path = write_fq(&tmp, "r2", &r2_recs);
        let fasta_path = tmp.path().join("out.fa");
        let mut cmd = make_detect(vec![r1_path, r2_path], Some(fasta_path.clone()));
        // Each of the three clusters accounts for ~1/3 of detections; a 0.90 cutoff
        // excludes every one.
        cmd.min_fraction = 0.90;
        let err = cmd.execute().unwrap_err().to_string();
        assert!(
            err.contains("No adapter reached --min-fraction"),
            "expected empty-FASTA hard-fail; got: {err}"
        );
        assert!(!fasta_path.exists(), "empty FASTA should not be created; found file");
    }

    #[test]
    fn se_output_fasta_errors_when_min_fraction_excludes_all_hits() {
        // Interleave TruSeq and Nextera tails so each candidate holds only ~half
        // the detections; a --min-fraction well above 0.5 then excludes both.
        let tmp = TempDir::new().unwrap();
        let cap = 18;
        let truseq_tail = &TRUSEQ.seq_r1[..cap];
        let nextera_tail = &NEXTERA.seq_r1[..cap];
        let mut recs: Vec<(String, Vec<u8>)> = Vec::new();
        for i in 0..50 {
            let mut a = template_of_len(70);
            a.extend_from_slice(truseq_tail);
            recs.push((format!("t_{i}"), a));
            let mut b = template_of_len(70);
            b.extend_from_slice(nextera_tail);
            recs.push((format!("n_{i}"), b));
        }
        let path = write_fq(&tmp, "r1", &recs);
        let fasta_path = tmp.path().join("out.fa");
        let mut cmd = make_detect(vec![path], Some(fasta_path.clone()));
        cmd.min_fraction = 0.90;
        let err = cmd.execute().unwrap_err().to_string();
        assert!(
            err.contains("No candidate adapter reached --min-fraction"),
            "expected SE empty-FASTA hard-fail; got: {err}"
        );
        assert!(!fasta_path.exists());
    }

    // ---- FASTA loader: empty-body header must not be silently dropped ----

    #[test]
    fn load_adapter_fasta_errors_on_header_without_body() {
        let tmp = TempDir::new().unwrap();
        let fa = tmp.path().join("bad.fa");
        // `>foo` immediately followed by `>bar` — foo has no body. Old behavior
        // silently dropped foo; new behavior errors so the user notices.
        std::fs::write(&fa, ">foo\n>bar\nACGTACGTACGT\n").unwrap();
        let err = match build_se_candidates(&[], &Some(fa), 10) {
            Ok(_) => panic!("expected empty-header error; loader accepted the malformed FASTA"),
            Err(e) => e.to_string(),
        };
        assert!(
            err.contains("empty sequence") || err.contains("header without a body"),
            "expected empty-header error; got: {err}"
        );
    }
}
