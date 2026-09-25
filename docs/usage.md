# Using chelae

The full reference for `chelae trim` and `chelae detect`: how input and output work, more examples, and every option. `chelae <command> --help` has the longest explanation of each option and is the authoritative reference.

- [Input and output](#input-and-output)
  - [Files, stdin and stdout](#files-stdin-and-stdout)
  - [Paired-end layout](#paired-end-layout)
  - [How a single input is sniffed](#how-a-single-input-is-sniffed)
  - [How pairs are checked](#how-pairs-are-checked)
  - [Output compression](#output-compression)
  - [Pipes that close early](#pipes-that-close-early)
- [`chelae trim`](#chelae-trim)
- [`chelae detect`](#chelae-detect)

## Input and output

Both subcommands read FASTQ the same way; `chelae trim` also writes it.

### Files, stdin and stdout

- `-i`/`--inputs` takes one or two paths, and `chelae trim`'s `-o`/`--outputs` takes one or two. Both default to `-`, meaning stdin or stdout, and `-` may appear at most once in each.
- Inputs may be plain text, gzip or BGZF. The format is detected from the content, not the file name, so compressed stdin works too.
- Two inputs must be different files (after resolving symlinks), two outputs must be different files, and no output may be one of the inputs.
- chelae won't read FASTQ from an interactive terminal (`stdin is a terminal; pass --inputs or pipe data in`), but it will write to one, so `chelae trim -i r1.fq.gz r2.fq.gz | head` works.
- `--metrics` and `--json` need file paths; they don't accept `-`.
- Progress and the end-of-run summary are logged to stderr, so stdout carries only FASTQ.

### Paired-end layout

The layout comes from the number of inputs and outputs; there is no interleaving flag.

| Inputs              | Outputs | Meaning                          |
|---------------------|---------|----------------------------------|
| 2 files             | 2       | split PE in → split PE out       |
| 2 files             | 1       | split PE in → interleaved out    |
| 1 file (sniffed PE) | 2       | interleaved in → split out       |
| 1 file (sniffed PE) | 1       | interleaved in → interleaved out |
| 1 file (sniffed SE) | 1       | single-end                       |
| 1 file (sniffed SE) | 2       | error                            |

- Two inputs are always R1 then R2, by position. They are never sniffed.
- One input is single-end unless its first records sniff as an interleaved pair (below).
- A completely empty single input carries no evidence either way, so the layout follows the number of outputs (or of `--read-structures` or `--adapter-sequence` values), and chelae writes valid empty output.
- Two separate input files are the fastest layout: R1 and R2 are decompressed in parallel on their own threads.

### How a single input is sniffed

chelae peeks at up to the first four records of a single input; nothing is lost to the peek. The input is interleaved paired-end if one of the naming conventions below explains records 1 and 2 as mate 1 then mate 2 and, when there are at least four records, explains records 3 and 4 as a pair under the same convention. Otherwise it is single-end. The decision and the convention chosen are logged.

The conventions, tried in this order:

| Convention      | Mate 1 name         | Mate 2 name         | Typical source                                  |
|-----------------|---------------------|---------------------|-------------------------------------------------|
| `/1` / `/2`     | `@frag/1`           | `@frag/2`           | older Illumina pipelines, `samtools fastq`      |
| `.1` / `.2`     | `@SRR123.7.1`       | `@SRR123.7.2`       | SRA `fastq-dump -I --split-spot`                |
| `_1` / `_2`     | `@frag_1`           | `@frag_2`           | some pipelines                                  |
| Identical names | `@frag 1:N:0:ACGT`  | `@frag 2:N:0:ACGT`  | Casava 1.8+, SRA and ENA defaults, bare names   |

A suffix convention compares the first whitespace-delimited word of each name and requires the part before the suffix to match. For identical names that first word must match exactly, and if both comments carry a mate number, mate 1's must be 1 and mate 2's must be 2. A mate number is a Casava 1.8+ field starting `1:` or `2:`, or a first comment word ending `/1` or `/2`, as in ENA's `@ERR000589.1 EAS139_45:5:1:2:111/1`.

Records 3 and 4 are checked because a single-end SRA file (`@SRR390728.1`, `@SRR390728.2`, `@SRR390728.3`, …) looks like a `.1`/`.2` pair in its first two records; its records 3 and 4 don't pair, so it is correctly read as single-end. An input with only two or three records is decided on records 1 and 2 alone.

If records 1 and 2, or records 3 and 4 under the convention records 1 and 2 chose, are a pair in mate 2, mate 1 order, chelae stops with an error rather than trim reversed interleaved input as single-end.

### How pairs are checked

Paired-end input is checked pair by pair for the whole run:

- **Interleaved input:** every pair must satisfy the convention chosen when sniffing. A pair that doesn't, or a stream that ends mid-pair (an odd number of records), is an error naming the pair and its record numbers in the file.
- **Two input files:** the convention is chosen from the first pair and every later pair must satisfy it, and the two files must end together. For the first pair:
  - names in mate 2, mate 1 order are an error that suggests the inputs were given the wrong way round;
  - names that both carry mate markers but don't correspond are an error;
  - names that match no convention and aren't both mate-marked get a warning, and the records are paired by position only; the files must still end together.

### Output compression

With the default `--output-compression auto`, outputs whose path ends in `.gz` or `.bgz` (in any case) are written as BGZF, and all others as plain text, including stdout. `--output-compression bgzf` or `none` forces one encoding for every output. `-c`/`--compression-level` (1–12, default 1) sets the BGZF level; higher levels cost much more CPU for slightly smaller files.

### Pipes that close early

If a downstream reader closes stdout before chelae has finished (e.g. `chelae trim | head`), chelae stops reading input and exits successfully with whatever it had written, without an error or a `SIGPIPE` death. If `--metrics` or `--json` is set, it warns at that point that their counts include records it processed but never wrote.

## `chelae trim`

`chelae trim` runs these steps in one pass:

1. Poly-G 3' trim (on by default)
2. Adapter trimming: paired-end overlap detection, confirmed against the adapter sequences given with `--kit`, `--adapter-sequence` or `--adapter-fasta`, which are also searched for directly in single-end reads and in inserts too short to overlap
3. [Read-structure](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures) based hard-trim and UMI extraction, after adapter trimming so that tail-skip segments act on the cleaned template
4. Optional poly-X 3' trim (`--trim-polyx`)
5. Optional 5'→3' and/or 3'→5' sliding-window quality trim
6. Length filter (`--filter-length MIN[:MAX]`)
7. Optional N-base, mean-quality and low-quality-fraction filters

If a read is shorter than its read structure's fixed-length segments require, the pair is dropped and counted as filtered on length. A fastp-compatible JSON report (`--json`) feeds MultiQC's fastp module unchanged.

### Examples

#### Trim single-end reads with the Nextera adapters

```bash
chelae trim -i sample.fq.gz -o trimmed.fq.gz --kit nextera
```

#### Add a 3' quality trim with an 8 bp sliding window at Q20

```bash
chelae trim \
    -i sample.r1.fq.gz sample.r2.fq.gz \
    -o trimmed.r1.fq.gz trimmed.r2.fq.gz \
    --kit truseq \
    --quality-trim-3p 8:20
```

#### Write a report for MultiQC

```bash
chelae trim \
    -i sample.r1.fq.gz sample.r2.fq.gz \
    -o trimmed.r1.fq.gz trimmed.r2.fq.gz \
    --kit truseq \
    --json sample.chelae.json
```

#### Trim an interleaved file to split R1/R2 files

```bash
chelae trim -i interleaved.fq.gz -o trimmed.r1.fq.gz trimmed.r2.fq.gz --kit truseq
```

### Options

#### Inputs, outputs, and runtime

| Option                          | Description                                                                                                  | Default |
|---------------------------------|--------------------------------------------------------------------------------------------------------------|---------|
| `-i, --inputs <PATHS>...`       | One or two FASTQ paths; `-` means stdin. Two files are split R1/R2; one is SE unless sniffed as interleaved PE. See [Input and output](#input-and-output) | `-`     |
| `-o, --outputs <PATHS>...`      | One or two output FASTQ paths; `-` means stdout. One output interleaves both mates; two write split R1/R2   | `-`     |
| `--output-compression <MODE>`   | `auto` (BGZF for `.gz`/`.bgz` paths, case-insensitive; plain text otherwise), `bgzf`, or `none` — forces the encoding for every output | `auto`  |
| `-t, --threads <N>`             | Number of threads to use                                                                                     | `4`     |
| `-c, --compression-level <1-12>`| Compression level for BGZF outputs; ignored for plain-text outputs                                           | `1`     |
| `-m, --metrics <PATH>`          | Optional path for the trimming metrics TSV (does not accept `-`); a summary is always logged to stderr       | —       |
| `-j, --json <PATH>`             | Optional fastp-shape JSON report (does not accept `-`); consumed by MultiQC's `fastp` module unchanged        | —       |

#### Read-structure (hard-trim + UMI extraction)

| Option                                | Description                                                                                                | Default |
|---------------------------------------|------------------------------------------------------------------------------------------------------------|---------|
| `-r, --read-structures <RS>...`       | Optional [read-structures](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures) per input; supports `T` (template), `M` (UMI → read name), `S` (skip); applied after adapter trim | —       |
| `--discard-unsupported-segments`      | Treat `B` (sample barcode) and `C` (cellular barcode) segments as `S` (skip) instead of erroring          | off     |

#### Adapter trimming

| Option                              | Description                                                                                                  | Default |
|-------------------------------------|--------------------------------------------------------------------------------------------------------------|---------|
| `-k, --kit <NAME>...`               | Built-in kit preset; repeatable. Known: `truseq`, `nextera`, `small-rna`, `aviti`, `mgi` (alias `dnbseq`), `all` | —       |
| `-a, --adapter-sequence <SEQ>...`   | 3' adapter sequence(s); 1 for SE, 1 or 2 for PE (R1, R2); ACGT or IUPAC                                      | —       |
| `-f, --adapter-fasta <PATH>`        | FASTA of adapter sequences; best match is trimmed                                                            | —       |
| `--adapter-min-length <N>`          | Minimum match length when searching the 3' end for an adapter sequence (SE mode; PE mode only for inserts < `overlap-min-length`) | `6`     |
| `--adapter-mismatch-rate <0..1>`    | Max fraction of mismatches when matching adapter against the 3' end (default ≈ 1 mismatch / 8 bases)         | `0.125` |

#### Paired-end overlap detection

| Option                                | Description                                                                                                | Default |
|---------------------------------------|------------------------------------------------------------------------------------------------------------|---------|
| `--no-overlap-detection`              | Disable PE-overlap trim-point detection (ignored for SE); rely on sequence matching alone                  | on (PE) |
| `--overlap-min-length <N>`            | Minimum overlap (bp) required to declare R1/R2 overlap                                                     | `30`    |
| `--overlap-max-mismatch-rate <0..1>`  | Max fraction of mismatches in the overlap probe window                                                     | `0.10`  |
| `--overlap-diagnostic-length <N>`     | When evaluating PE overlap, only examine this many overlapping bases. Multiples of 16 ideal.               | `64`    |
| `--expected-insert-size <BP>`         | Hint for typical insert size; seeds the overlap candidate-walk order so the right overlap is found sooner  | —       |
| `--insert-size-stats`                 | Emit a fastp-shape per-pair insert-size histogram under `insert_size` in the JSON (extends overlap probing to I > R configurations) | off     |

#### Poly-G / poly-X trimming

| Option                  | Description                                                                                                          | Default |
|-------------------------|----------------------------------------------------------------------------------------------------------------------|---------|
| `--trim-polyg <N>`      | 3' poly-G trim minimum run length; pass `0` to disable                                                               | `10`    |
| `--trim-polyx [<N>]`    | Enable 3' poly-X trim (A/C/T homopolymer tails, e.g. poly-A from RNA-seq) with the given minimum run length         | off     |

#### Quality trimming

Both quality-trim modes shorten the read at the 3' end; the `-3p` / `-5p` suffix is the scan direction, not the trim location. `-3p` is conservative (keeps everything up to the last good window from the 3' end); `-5p` is aggressive (cuts at the first bad window encountered from the 5' end).

| Option                          | Description                                                                                                  | Default      |
|---------------------------------|--------------------------------------------------------------------------------------------------------------|--------------|
| `--quality-trim-3p [<W:Q>]`     | Scan 3'→5'; trim trailing bases until a window of size `W` has mean quality ≥ `Q` (fastp `--cut_tail`)        | off (`8:20`) |
| `--quality-trim-5p [<W:Q>]`     | Scan 5'→3'; truncate at the first window of size `W` with mean quality < `Q` (fastp `--cut_right`)            | off (`8:20`) |

#### Filters (applied after trimming; pair dropped if either mate fails)

| Option                            | Description                                                                                              | Default |
|-----------------------------------|----------------------------------------------------------------------------------------------------------|---------|
| `-l, --filter-length <MIN[:MAX]>` | Drop reads/pairs with post-trim length below `MIN` (or above `MAX`)                                      | `15`    |
| `--filter-max-ns <N>`             | Drop reads/pairs whose per-mate count of ambiguous (N) bases exceeds `N`                                 | off     |
| `--filter-mean-qual <Q>`          | Drop reads/pairs whose post-trim mean Phred quality is below `Q` (runs last, after every trim stage)     | off     |
| `--filter-low-qual <Q:F>`         | Drop reads/pairs where the fraction of bases below quality `Q` exceeds `F` (e.g. `15:0.4`)               | off     |

## `chelae detect`

`chelae detect` samples reads from one or two FASTQ files and reports the adapter sequences present. Paired-end input discovers adapters from the reads that run past the end of their insert into adapter, found by R1/R2 overlap, so it needs no knowledge of kits. Single-end input scores reads against every built-in kit plus any candidates you supply. The discovered or winning sequences can be written as FASTA for `chelae trim --adapter-fasta`. Input follows the same rules as `chelae trim` ([Input and output](#input-and-output)), and `-o -` writes the FASTA to stdout.

### Examples

#### Score a single-end library against extra candidates as well as the built-in kits

```bash
chelae detect \
    -i sample.fq.gz \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -f custom-adapters.fa
```

### Options

Most options apply to paired-end or single-end input only, as each row says. Run `chelae detect --help` for the full rationale.

| Option                                | Description                                                                                                | Default      |
|---------------------------------------|------------------------------------------------------------------------------------------------------------|--------------|
| `-i, --inputs <PATHS>...`             | One or two FASTQ paths; `-` means stdin. Two files are split R1/R2; one is SE unless sniffed as interleaved PE | `-`          |
| `-o, --output-fasta <PATH>`           | Optional FASTA output of discovered/winning adapter(s); `-` writes to stdout; ready to feed back into `chelae trim --adapter-fasta` | —            |
| `-a, --adapter-sequence <SEQ>...`     | (SE only) Extra adapter candidate(s) to score against, in addition to every built-in kit                    | —            |
| `-f, --adapter-fasta <PATH>`          | (SE only) FASTA of extra adapter candidates; record names are preserved in the report                       | —            |
| `-n, --num-detections <N>`            | Target number of usable detections before stopping. Higher = more confident composition estimate            | `5000`       |
| `--max-reads <N>`                     | Hard cap on records scanned even if `--num-detections` isn't reached                                        | `1000000`    |
| `--min-detections-for-report <N>`     | Refuse to report if final detection count falls below this floor (avoids confident-looking tiny samples)    | `20`         |
| `--min-fraction <0..1>`               | Minimum share of detections an adapter must account for to be reported                                      | `0.05`       |
| `--min-tail-length <N>`               | Minimum length of adapter evidence (bp) per detection (PE post-template tail; SE matched alignment)         | `8`          |
| `--overlap-min-length <N>`            | (PE) Minimum overlap (bp) required for PE-overlap detection                                                 | `30`         |
| `--overlap-max-mismatch-rate <0..1>`  | (PE) Max fraction of mismatches in the overlap probe                                                        | `0.10`       |
| `--overlap-diagnostic-length <N>`     | (PE) Upper bound on the probe length per overlap-length candidate (bp)                                      | `64`         |
| `--adapter-min-length <N>`            | (SE) Minimum match length (bp) when scoring a candidate against a read's 3' end                             | `10`         |
| `--adapter-mismatch-rate <0..1>`      | (SE) Max fraction of mismatches when matching a candidate against a read's 3' end                           | `0.125`      |
| `--trim-polyg <N>`                    | 3' poly-G trim min run length applied before the probe (cleans 2-color "no signal" tails); `0` disables     | `10`         |
| `--trim-polyx <N>`                    | 3' poly-X (A/C/T) trim min run length applied before the probe (more aggressive than trim's `10`); `0` disables | `5`        |
| `--quality-trim <W:Q>`                | 3' cut-right quality trim applied before the probe; pass `off`/`none`/`no` to disable                       | `4:20`       |
