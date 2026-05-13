#!/usr/bin/env bash
# Run the chelae benchmark pipeline.
#
# Defaults to config/accuracy.config.yaml (which itself points at
# config/accuracy.samples.tsv via the `samples_tsv` key). Pass an
# alternate configfile / samples sheet positionally, or use the
# performance config:
#
#   ./run.sh                                       # accuracy run (default)
#   ./run.sh config/smoke.config.yaml              # smoke test (~15-20 min)
#   ./run.sh config/performance.config.yaml        # performance run
#   ./run.sh path/to/cfg.yaml path/to/samples.tsv  # alt config + samples
#   ./run.sh --dry-run                             # preview job graph
#   ./run.sh --cores 16                            # cap cores (default: all)
#
# Anything after `--` is forwarded verbatim to snakemake.
set -euo pipefail

cd "$(dirname "$0")"
PIPELINE_DIR="$PWD"

# If install.sh built an isolated rust toolchain, source its env so any
# rebuild triggered by snakemake (or a developer running the script) uses
# it instead of the host's cargo.
if [[ -d "$PIPELINE_DIR/.rust/cargo" ]]; then
  export RUSTUP_HOME="$PIPELINE_DIR/.rust/rustup"
  export CARGO_HOME="$PIPELINE_DIR/.rust/cargo"
  export PATH="$CARGO_HOME/bin:$PATH"
fi

CONFIG_FILE=""
SAMPLES_FILE=""
DRY_RUN=0
CORES="all"
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run|-n) DRY_RUN=1; shift ;;
    --cores)      CORES="$2"; shift 2 ;;
    --cores=*)    CORES="${1#--cores=}"; shift ;;
    --) shift; EXTRA_ARGS+=("$@"); break ;;
    -h|--help)
      sed -n '2,15p' "$0"
      exit 0
      ;;
    -*)
      echo "unknown flag: $1" >&2
      exit 2
      ;;
    *)
      if [[ -z "$CONFIG_FILE" ]]; then
        CONFIG_FILE="$1"
      elif [[ -z "$SAMPLES_FILE" ]]; then
        SAMPLES_FILE="$1"
      else
        echo "too many positional args: $1" >&2
        exit 2
      fi
      shift
      ;;
  esac
done

CONFIG_FILE="${CONFIG_FILE:-config/accuracy.config.yaml}"
[[ -f "$CONFIG_FILE" ]] || { echo "config not found: $CONFIG_FILE" >&2; exit 1; }
if [[ -n "$SAMPLES_FILE" ]]; then
  [[ -f "$SAMPLES_FILE" ]] || { echo "samples not found: $SAMPLES_FILE" >&2; exit 1; }
fi

SNAKE_ARGS=(
  --configfile "$CONFIG_FILE"
  --cores "$CORES"
  # bench=100 sizes the exclusive-lock pool: trim_one consumes all 100 so
  # nothing else can run during a timed measurement; non-trim rules
  # consume 1 each so they can otherwise run in parallel.
  --resources bench=100
  --rerun-incomplete
)
if [[ -n "$SAMPLES_FILE" ]]; then
  # Use absolute path so snakemake can resolve it regardless of where it
  # was launched from (the Snakefile resolves relative paths against
  # benchmark-pipeline/, but absolute is unambiguous).
  ABS_SAMPLES="$(cd "$(dirname "$SAMPLES_FILE")" && pwd)/$(basename "$SAMPLES_FILE")"
  SNAKE_ARGS+=(--config "samples_tsv=$ABS_SAMPLES")
fi
if [[ "$DRY_RUN" -eq 1 ]]; then
  SNAKE_ARGS+=(-n -p)
fi
# Bash 3.2 (macOS default) errors under `set -u` on empty-array expansion;
# the `${name[@]+...}` form expands only if name was set.
SNAKE_ARGS+=(${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"})

echo "==> pixi run snakemake ${SNAKE_ARGS[*]}"
exec pixi run snakemake "${SNAKE_ARGS[@]}"
