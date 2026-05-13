#!/usr/bin/env python3
"""Count reads and bases across one or more gzipped FASTQ files. Output is a
two-line stats file: `reads\tN` and `bases\tM`."""

import gzip
import sys


def count(path: str) -> tuple[int, int]:
    reads = 0
    bases = 0
    with gzip.open(path, "rb") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip(b"\n")
            fh.readline()  # +
            fh.readline()  # qual
            reads += 1
            bases += len(seq)
    return reads, bases


def main():
    total_reads = 0
    total_bases = 0
    for path in sys.argv[1:]:
        r, b = count(path)
        total_reads += r
        total_bases += b
    print(f"reads\t{total_reads}")
    print(f"bases\t{total_bases}")


if __name__ == "__main__":
    main()
