//! Crate-level library surface shared by `chelae` subcommands.
//!
//! Currently exposes:
//! - [`IUPAC_MASKS`] — a 256-entry byte LUT mapping ASCII base characters (and the
//!   subset of IUPAC ambiguity codes) onto 4-bit masks. Used by the adapter matcher
//!   to evaluate IUPAC-aware base compatibility: two bases are compatible iff
//!   `IUPAC_MASKS[a] & IUPAC_MASKS[b] != 0`.
//! - [`adapter_db`] — built-in 3' adapter presets for common sequencing kits,
//!   selectable via `--kit`.

pub mod adapter_db;

/// ASCII-indexed IUPAC ambiguity masks. Bit 0 = A, bit 1 = C, bit 2 = G, bit 3 = T
/// (U is aliased to T). Case-insensitive: both uppercase and lowercase entries are
/// populated. Non-IUPAC bytes map to 0 — callers compute compatibility via
/// `MASKS[a] & MASKS[b] != 0`, and a 0 mask always yields "incompatible".
pub const IUPAC_MASKS: [u8; 256] = {
    let mut masks = [0u8; 256];
    let (a, c, g, t) = (1, 2, 4, 8);
    let pairs: &[(u8, u8)] = &[
        (b'A', a),
        (b'C', c),
        (b'G', g),
        (b'T', t),
        (b'U', t),
        (b'M', a | c),
        (b'R', a | g),
        (b'W', a | t),
        (b'S', c | g),
        (b'Y', c | t),
        (b'K', g | t),
        (b'V', a | c | g),
        (b'H', a | c | t),
        (b'D', a | g | t),
        (b'B', c | g | t),
        (b'N', a | c | g | t),
    ];
    let mut i = 0;
    while i < pairs.len() {
        let (upper, mask) = pairs[i];
        masks[upper as usize] = mask;
        // Case-fold to lowercase (ASCII: `| 0x20`).
        masks[(upper | 0x20) as usize] = mask;
        i += 1;
    }
    masks
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn acgt_masks_are_one_hot() {
        assert_eq!(IUPAC_MASKS[b'A' as usize], 0b0001);
        assert_eq!(IUPAC_MASKS[b'C' as usize], 0b0010);
        assert_eq!(IUPAC_MASKS[b'G' as usize], 0b0100);
        assert_eq!(IUPAC_MASKS[b'T' as usize], 0b1000);
        assert_eq!(IUPAC_MASKS[b'U' as usize], 0b1000, "U aliases T");
    }

    #[test]
    fn lowercase_matches_uppercase() {
        for b in b"ACGTUMRWSYKVHDBN" {
            assert_eq!(
                IUPAC_MASKS[*b as usize],
                IUPAC_MASKS[(*b | 0x20) as usize],
                "case mismatch for {:?}",
                *b as char,
            );
        }
    }

    #[test]
    fn n_matches_everything() {
        for b in b"ACGT" {
            assert_ne!(IUPAC_MASKS[*b as usize] & IUPAC_MASKS[b'N' as usize], 0);
        }
    }

    #[test]
    fn non_iupac_bytes_are_zero() {
        assert_eq!(IUPAC_MASKS[b'Z' as usize], 0);
        assert_eq!(IUPAC_MASKS[b'.' as usize], 0);
        assert_eq!(IUPAC_MASKS[0_usize], 0);
    }

    #[test]
    fn r_is_superset_of_a_and_g_only() {
        let r = IUPAC_MASKS[b'R' as usize];
        assert_ne!(r & IUPAC_MASKS[b'A' as usize], 0);
        assert_ne!(r & IUPAC_MASKS[b'G' as usize], 0);
        assert_eq!(r & IUPAC_MASKS[b'C' as usize], 0);
        assert_eq!(r & IUPAC_MASKS[b'T' as usize], 0);
    }
}
