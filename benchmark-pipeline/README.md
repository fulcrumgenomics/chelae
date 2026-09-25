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
  by `workflow/scripts/eval_accuracy.py`. Includes edge buckets for reads dropped by the tool.
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

Every tool is at the latest version on bioconda (September 2026). Adding more is a matter of adding a `tools/<name>/render.py` and an entry in `tools.yaml`.

| Tool            | Version | Source         | Notes                                        |
|-----------------|---------|----------------|----------------------------------------------|
| chelae          | HEAD    | this workspace | built from `../` by `install.sh`             |
| fastp           | 1.3.7   | bioconda       |                                              |
| cutadapt        | 5.2     | bioconda       | adapter-focused; underlies Perl Trim Galore  |
| trim-galore-rs  | 2.3.0   | bioconda       | Trim Galore's Rust rewrite; own pixi env     |
| trimmomatic     | 0.41    | bioconda       | no output compression-level option           |
| bbmap (bbduk)   | 40.02   | bioconda       |                                              |
| adapterremoval  | 3.0.2   | bioconda       | `adapterremoval3`                            |

The Perl Trim Galore (0.6.x) is not included: it runs cutadapt, so its accuracy is cutadapt's, and it was the slowest tool by far. atria has no bioconda package.

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

- **`config/smoke.config.yaml`**: every tool, 2 samples at 0.1x WGS, 4 threads, 1 replicate. A short end-to-end run that exercises every code path. **Run this before any bigger sweep** to validate render.py / pixi env / DAG changes.
- **`config/accuracy.config.yaml`** *(default for `run.sh`)*, tier 1: every tool on 11 scenarios spanning insert/read-length geometry, adapter kit and error rate, at 0.1x WGS each (~1M pairs for 2x150), 8 threads, 1 replicate. Per-read accuracy scoring on the `adapter_only` rows; the timings are only a screen for which tools go on to tier 2.
- **`config/performance.config.yaml`**, tier 2: the tools that are competitive on tier 1's paired-end runs, on two paired-end scenarios at ~50M pairs each, `wgs` trim config, 8 threads, 3 replicates. These are the runtime numbers.

Both tiers target an 8-core host, so a tool that fans out to every available core gets no more than the others. Run them into separate `results/` directories (move `results/` aside between them).

### RAM-backed scratch

Set `scratch_dir` to a tmpfs mount and no timed run touches the disk: each sample's inputs are copied there before its trims, and every trimmed FASTQ is written there. Samples are staged one at a time in sample-sheet order, so scratch needs room for one sample's inputs plus one run's outputs (~30 GB for a 50M-pair 2x150 sample at compression level 1). Everything else, including simulated inputs, stays under `results/`.

## Quick start

```bash
cd benchmark-pipeline
./install.sh                                    # one-time setup
# edit config/accuracy.config.yaml: point reference: at an indexed FASTA
./run.sh config/smoke.config.yaml               # smoke test
./run.sh --dry-run                              # preview the accuracy job graph
./run.sh                                        # accuracy run (default)
./run.sh config/performance.config.yaml         # performance run
./run.sh path/to/cfg.yaml path/to/samples.tsv   # alt inputs
./run.sh -- --config scratch_dir=/mnt/ram       # stage inputs/outputs on a tmpfs
# Plot once a bench.tsv exists
pixi run Rscript workflow/scripts/plot.R results/bench.tsv results/plots
```

`install.sh` flags:
- `--system-rust` — use the cargo on PATH instead of installing rustup into `.rust/`
- `--skip-build` — skip the chelae build (re-run for env-only changes)

`run.sh` accepts `--dry-run`, `--cores N`, and forwards any args after `--` straight to snakemake. It runs with `--keep-going`, and each trim is killed and fails its job after `trim_timeout_minutes`, so one hung tool doesn't stall the sweep. With `scratch_dir` set, though, a sample whose trims fail blocks the samples after it, since each waits on the previous sample's trims.

## Output

- `results/sim/<sample>/` — holodeck outputs (`r1.fastq.gz`, `r2.fastq.gz`)
- `results/trim/<sample>/<trim_config>/<tool>/t<threads>/rep<N>/` — trimmed FASTQ (under `scratch_dir` instead when set) + tool log + `time.txt` + (optional) `perf.txt`
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
