# chelae benchmark pipeline

A reproducible Snakemake pipeline that benchmarks `chelae` against other FASTQ
trimming tools on simulated (and later, real) short-read data, recording
throughput, resource use, and — on configs that disable quality-based trimming —
per-read accuracy.

Designed to run on a single machine with the box to itself: benchmarks are
serialized via a Snakemake resource so only one trim job executes at a time,
even though data-prep (simulation, reference fetch, merging) runs in parallel.

## What it measures

For each `(sample, trim_config, tool, threads, replicate)` combination the
pipeline records:

- **Performance** — wall time, user/sys CPU seconds, max RSS (GNU `time -v`);
  optionally `perf stat` counters (cycles, instructions, task-clock) where
  available.
- **Throughput** — reads/s and bases/s (derived from wall time and input
  record counts). Thread-normalized throughput (reads/s/thread) is computed at
  aggregation time.
- **Accuracy** *(for `adapter_only` configs only)* — a melted
  `expected_trim_len × observed_trim_len` histogram per read, emitted as a TSV
  by `chelae bench-eval`. Includes edge buckets for reads dropped by the tool.
- **Provenance** — tool name/version, host CPU/arch/mem, kernel, date.

Results land in a single wide TSV per run (`results/bench.tsv`); merging runs
across hosts is just `cat` after deduping headers. All plotting reads from that
merged TSV.

## Design at a glance

```
config/
  config.yaml              pipeline-level knobs (thread counts, replicates, paths)
  samples.tsv              sample sheet; one row per (name, depth, fragment, read_len, error_rate)
  trim_configs/            semantic "how to trim" presets
    adapter_only.yaml      adapter trim only — for accuracy scoring
    wgs.yaml               adapter + quality + length filter — for performance
  tools.yaml               tool inventory with pinned versions
  adapters.yaml            named adapter sets used by trim_configs
tools/
  chelae/render.py         each tool has a render.py that reads the
  fastp/render.py          semantic trim_config + runtime knobs and emits
  cutadapt/render.py       the tool-specific argv. These are the canonical
  ...                      documentation of how each tool was configured.
workflow/
  rules/                   Snakemake modules: simulate, trim, eval, aggregate
  scripts/                 host_info, gnu-time parser, merge, plot.R
```

**Separation of accuracy and performance runs.** Quality-based trimming is
tool-specific and would conflate adapter-trim accuracy with quality-trim
heuristics. So the pipeline runs each sample twice:

1. Against a `adapter_only` config (adapter detection/trim only; no quality
   trim, no length filter, no N filter). Output is scored for accuracy.
2. Against a realistic config (e.g. `wgs`, which enables adapter + quality +
   length filter). Output is measured for performance but not scored.

Both runs produce timing numbers; only the first feeds the accuracy eval.

## Tool matrix

Pinned to match [nf-core/modules](https://github.com/nf-core/modules) at the
time of writing (April 2026). These are the trimmers most commonly used in
real nf-core pipelines today; adding more is a matter of adding a
`tools/<name>/render.py` and an entry in `tools.yaml`.

| Tool            | Version | Source        | Notes                                  |
|-----------------|---------|---------------|----------------------------------------|
| chelae          | HEAD    | this workspace| built from `../` via `cargo install`   |
| fastp           | 1.1.0   | bioconda      | nf-core default for most pipelines     |
| cutadapt        | 5.2     | bioconda      | adapter-focused; underlies Trim Galore |
| trim-galore     | 0.6.10  | bioconda      | wraps cutadapt; nf-core rnaseq default |
| trimmomatic     | 0.39    | bioconda      | legacy but still in nf-core/modules    |
| bbmap (bbduk)   | 39.18   | bioconda      | nf-core alternative in mag/taxprofiler |
| adapterremoval  | 2.3.4   | bioconda      | nf-core alternative in mag/taxprofiler |
| atria           | current | bioconda      | newer, Julia, claims SIMD speedups     |

## Requirements

- **pixi** ≥ 0.40 — `install.sh` installs it for you if missing
- a reference FASTA — either a pre-indexed local path, or any http/https/ftp
  URL. URL references get downloaded to `results/reference/` on first run;
  `.gz` URLs are decompressed; `.fai` is fetched alongside the URL if
  available, else built with `samtools faidx` from the pixi env. Set in
  `config.yaml` (or overridden per-run via `--config reference=...`).

Everything else — pixi, the rust toolchain (sandboxed under `.rust/`),
`cargo-multivers`, every trimmer, GNU `time`, samtools — is installed by
`install.sh`. holodeck is `cargo install`-ed into the same sandbox to
guarantee a specific version on every platform (the bioconda 0.2.0
upload is missing linux-64; 0.2.1 has the IUPAC-only-bases fix we want
anyway). The host's `~/.cargo` and `~/.rustup` are not touched.

## Three configs: smoke, accuracy, performance

The pipeline ships three complementary configurations:

- **`config/smoke.config.yaml`** — every tool, 2 samples at 0.1x WGS,
  4 threads, 1 replicate. ~15–20 minute end-to-end run that exercises
  every code path. **Run this before any bigger sweep** to validate
  render.py / pixi env / DAG changes.
- **`config/accuracy.config.yaml`** *(default for `run.sh`)* — every tool,
  10 scenarios spanning insert/read-length geometry and error rate, 2x WGS
  per scenario, 8 threads, 1 replicate. Wide tool stratification + per-read
  accuracy scoring on the `adapter_only` rows.
- **`config/performance.config.yaml`** — fast-tier tools only (chelae,
  fastp 1.3.2, fastp 1.1.0, cutadapt, adapterremoval), 3 scenarios at 30x
  WGS, full thread sweep `[1, 2, 4, 8]`, 3 replicates. Throughput +
  threading curves at realistic data scale.

The `accuracy` and `performance` configs target an 8-core host — threads
are capped at 8 so a tool that silently fans out to all available cores
can't game the comparison.

## Quick start

```bash
cd benchmark-pipeline
./install.sh                                    # one-time setup
# edit config/accuracy.config.yaml: point reference: at an indexed FASTA
./run.sh config/smoke.config.yaml               # smoke test (~15-20 min)
./run.sh --dry-run                              # preview the accuracy job graph
./run.sh                                        # accuracy run (default)
./run.sh config/performance.config.yaml         # performance run
./run.sh path/to/cfg.yaml path/to/samples.tsv   # alt inputs
# Plot once a bench.tsv exists
pixi run Rscript workflow/scripts/plot.R results/bench.tsv results/plots
```

`install.sh` flags:
- `--system-rust` — use the cargo on PATH instead of installing rustup into `.rust/`
- `--skip-build` — skip the chelae build (re-run for env-only changes)

`run.sh` accepts `--dry-run`, `--cores N`, and forwards any args after `--`
straight to snakemake.

## Output

- `results/sim/<sample>/` — holodeck outputs (`r1.fastq.gz`, `r2.fastq.gz`)
- `results/trim/<sample>/<trim_config>/<tool>/t<threads>/rep<N>/` — trimmed FASTQ + tool log + `time.txt` + (optional) `perf.txt`
- `results/eval/<sample>/<tool>/t<threads>/rep<N>/matrix.tsv` — melted accuracy matrix (adapter-only configs only)
- `results/bench.tsv` — one row per execution: tool, version, threads, wall, user, sys, rss, reads/s, bases/s, host fields
- `results/accuracy.tsv` — one row per `(run, expected_len, observed_len)` bucket
- `results/bench_summary.tsv` — bench.tsv with replicates collapsed (median/min/max/IQR per cell)
- `results/accuracy_summary.tsv` — headline accuracy stats per `(sample, trim_config, tool, threads, rep, mate)`: exact-match rate, false-positive/negative counts, MAE, RMSE
- `results/plots/` — PDF bundle: throughput / wall / speedup / max-RSS vs threads; accuracy heatmap

## Limits / caveats

- Accuracy scoring assumes **all trimming is from the 3' end**. Tools that 5'-trim (quality or fixed) mix positions, and the matrix can't distinguish. Use adapter-only configs when evaluating accuracy.
- Holodeck adapter boundary is encoded as `FRAG_LEN` in read names: for a read of length `L`, bases `[FRAG_LEN .. L)` are adapter (possibly N-padded). When `FRAG_LEN >= L`, the read has no adapter and any trimming is a false positive.
- Dropped reads (filtered by min-length, etc.) are scored as "observed = read_length" — everything trimmed. This is the right call for adapter-only configs because the only reason to drop a read there is the render script bug or tool behavior worth flagging.
- `perf stat` requires `kernel.perf_event_paranoid ≤ 2` (default on most Linux hosts). If unavailable the pipeline skips it and continues.
- macOS lacks GNU time by default; pixi installs the `time` package which provides it. Likewise `perf` is Linux-only.
