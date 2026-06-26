//! Small utilities shared across subcommands.

use anyhow::{Result, anyhow};
use fgoxide::io::Io;
use seq_io::fastq::Reader as FastqReader;
use std::io::BufRead;
use std::path::Path;

/// BufReader / BufWriter capacity used by the FASTQ I/O paths in every
/// subcommand. 512 KiB chosen via a 32k→4MiB sweep on Graviton4 (Neoverse-V2) and
/// x86 Granite Rapids (c8i): wall time and cycle count are flat across 256k–2MiB
/// on both architectures, with instruction count showing a shallow U-shape that
/// bottoms out at 512k–1024k. 512k is at the floor everywhere tested and halves
/// resident memory per reader/writer vs 1 MiB.
pub(crate) const BUFFER_SIZE: usize = 512 * 1024;

/// Opens every input path as a [`FastqReader`] with `BUFFER_SIZE` capacity.
/// fgoxide's `Io` pool runs gzip decompression on a background thread per
/// reader, sized to match the input count (benchmarks showed more threads here
/// don't help — decompression isn't the bottleneck).
pub(crate) fn open_fastq_inputs(
    paths: &[std::path::PathBuf],
) -> Result<Vec<FastqReader<Box<dyn BufRead + Send>>>> {
    let fgio = Io::new(paths.len().max(1) as u32, BUFFER_SIZE);
    paths
        .iter()
        .map(|p: &std::path::PathBuf| {
            let p_ref: &Path = p.as_ref();
            fgio.new_reader(p_ref)
                .map(|r| FastqReader::with_capacity(r, BUFFER_SIZE))
                .map_err(|e| anyhow!("Failed to open input {p_ref:?}: {e}"))
        })
        .collect()
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
}
