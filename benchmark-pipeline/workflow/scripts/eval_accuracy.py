#!/usr/bin/env python3
"""Compute the adapter-trim accuracy matrix for one (sample, tool, threads, rep)
run by comparing an original (holodeck-simulated) FASTQ against a trimmed one.

Holodeck encodes ground truth in each read name as colon-separated fields. The
format uses `::` as the delimiter between the top-level fields so that contig
names containing `:` (HLA alleles, etc.) parse cleanly:

    @holodeck::READ_NUM::FRAG_LEN::CONTIG::POS1+STRAND::POS2+STRAND::HAP::ERRS1::ERRS2

For a read of length L with template length F:
  - F >= L  ⇒ the read contains no adapter; expected_trim_len = 0
  - F <  L  ⇒ bases [F..L) are adapter/padding; expected_trim_len = L - F

Observed trim length is L - len(trimmed.seq). Reads present in the original
FASTQ but absent from the trimmed FASTQ were filtered out by the tool; those
count as `dropped` rows (observed_trim_len = L, dropped = True).

Scoring walks the sim and trim FASTQs in lockstep — every tool we benchmark
preserves input read order (bbduk needs `ordered=t`; others do so natively),
so dropped reads are detected when the trim stream's next read_num exceeds
the sim stream's. That keeps memory O(1) regardless of dataset size.

Limitation: assumes all trimming is 3'. A 5'-trimming tool would mis-score.
"""

import argparse
import gzip
import sys
from collections import Counter


class _MissingFragLen(Exception):
    """Raised when the holodeck read-name schema doesn't carry FRAG_LEN.
    The pipeline catches this and produces an empty matrix with a note
    rather than failing the whole run."""


def _open_fastq(path: str):
    return gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")


def _parse_name_fields(header: bytes) -> tuple[int, int]:
    """Return (read_num, frag_len) from a holodeck-encoded header line."""
    name = header[1:].split(None, 1)[0]
    parts = name.split(b"::")
    if len(parts) < 3 or parts[0] != b"holodeck":
        raise _MissingFragLen(
            f"header does not match holodeck `::` schema (need FRAG_LEN at "
            f"field 2 via `::`): {header!r}"
        )
    return int(parts[1]), int(parts[2])


def _iter_fastq(path: str):
    """Yield (read_num, frag_len, seq_len) per record."""
    fh = _open_fastq(path)
    try:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline().rstrip(b"\n")
            fh.readline()  # '+'
            fh.readline()  # qual
            rn, fl = _parse_name_fields(header.rstrip(b"\n"))
            yield rn, fl, len(seq)
    finally:
        fh.close()


def score_mate(sim_path: str, trim_path: str, read_len: int, mate: str,
               buckets: Counter) -> None:
    """Walk sim and trim FASTQs in lockstep, tallying accuracy buckets.
    Assumes both files are in the same read_num order and the trim file
    contains a subset of sim's reads (tools drop reads; they never invent)."""
    sim_iter = _iter_fastq(sim_path)
    trim_iter = _iter_fastq(trim_path)

    trim_rec = next(trim_iter, None)
    for sim_rn, sim_fl, _sim_seq_len in sim_iter:
        expected = max(0, read_len - sim_fl)
        if trim_rec is not None and trim_rec[0] == sim_rn:
            observed = max(0, read_len - trim_rec[2])
            buckets[(mate, expected, observed, False)] += 1
            trim_rec = next(trim_iter, None)
        else:
            # Either no trim record left, or trim is ahead of sim (which means
            # the tool dropped the current sim read).
            if trim_rec is not None and trim_rec[0] < sim_rn:
                raise RuntimeError(
                    f"trim stream out of order: sim read_num={sim_rn}, "
                    f"trim read_num={trim_rec[0]} (tool must preserve input "
                    f"order; bbduk needs `ordered=t`)"
                )
            buckets[(mate, expected, read_len, True)] += 1

    if trim_rec is not None:
        raise RuntimeError(
            f"trim stream has records beyond sim: first extra read_num="
            f"{trim_rec[0]}"
        )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sim-r1", required=True)
    p.add_argument("--sim-r2", required=True)
    p.add_argument("--trim-r1", required=True)
    p.add_argument("--trim-r2", required=True)
    p.add_argument("--paired", required=True)
    p.add_argument("--r1-len", type=int, required=True)
    p.add_argument("--r2-len", type=int, required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()
    paired = str(args.paired).lower() == "true"

    buckets: Counter = Counter()
    note = ""
    try:
        score_mate(args.sim_r1, args.trim_r1, args.r1_len, "r1", buckets)
        if paired:
            score_mate(args.sim_r2, args.trim_r2, args.r2_len, "r2", buckets)
    except _MissingFragLen as e:
        note = f"# SKIPPED: {e}\n"
        sys.stderr.write(note)

    with open(args.out, "w") as fh:
        if note:
            fh.write(note)
        fh.write("mate\texpected_trim_len\tobserved_trim_len\tdropped\tcount\n")
        for (mate, exp, obs, dropped), n in sorted(buckets.items()):
            fh.write(f"{mate}\t{exp}\t{obs}\t{int(dropped)}\t{n}\n")


if __name__ == "__main__":
    main()
