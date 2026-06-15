#!/usr/bin/env bash
# One-shot setup for the chelae benchmark pipeline.
#
# Idempotent: re-running is safe and only re-does steps that are out of date.
# Sets up an isolated environment under benchmark-pipeline/.rust/ so the
# host's ~/.cargo and ~/.rustup are never touched.
#
#   ./install.sh                  # default: install everything
#   ./install.sh --system-rust    # use whatever cargo is on PATH instead of
#                                 # installing rustup into .rust/
#   ./install.sh --skip-build     # set up envs but don't (re)build chelae
#
# After this script succeeds, run benchmarks with `./run.sh`.
set -euo pipefail

cd "$(dirname "$0")"
PIPELINE_DIR="$PWD"
REPO_ROOT="$(cd .. && pwd)"

USE_ISOLATED_RUST=1
DO_BUILD=1
for arg in "$@"; do
  case "$arg" in
    --system-rust) USE_ISOLATED_RUST=0 ;;
    --skip-build)  DO_BUILD=0 ;;
    -h|--help)
      sed -n '2,14p' "$0"
      exit 0
      ;;
    *)
      echo "unknown flag: $arg" >&2
      exit 2
      ;;
  esac
done

log() { printf '\033[1;34m==>\033[0m %s\n' "$*"; }

# ---- Build toolchain ------------------------------------------------------
# Slim Linux AMIs (AL2023 minimal etc.) ship without a C/C++ toolchain or
# cmake, which breaks rustc's build scripts and any Rust dep that drives a
# C-library build. We don't want to silently sudo-install system packages
# behind the user's back, so check up-front and bail loudly with concrete
# install instructions for the platform.
missing=()
{ command -v cc    >/dev/null 2>&1 || command -v gcc >/dev/null 2>&1; } || missing+=("cc/gcc")
{ command -v c++   >/dev/null 2>&1 || command -v g++ >/dev/null 2>&1; } || missing+=("c++/g++")
command -v cmake >/dev/null 2>&1 || missing+=("cmake")

if (( ${#missing[@]} > 0 )); then
  echo
  echo "ERROR: required build tools missing: ${missing[*]}" >&2
  echo "These are needed by rustc, cargo-multivers, and holodeck's build scripts." >&2
  echo >&2
  case "$(uname -s)" in
    Linux)
      if command -v dnf >/dev/null 2>&1; then
        echo "  Install with: sudo dnf install -y gcc gcc-c++ cmake" >&2
      elif command -v apt-get >/dev/null 2>&1; then
        echo "  Install with: sudo apt-get install -y build-essential cmake" >&2
      elif command -v yum >/dev/null 2>&1; then
        echo "  Install with: sudo yum install -y gcc gcc-c++ cmake" >&2
      elif command -v pacman >/dev/null 2>&1; then
        echo "  Install with: sudo pacman -S --needed base-devel cmake" >&2
      else
        echo "  Install gcc, g++, and cmake via your distribution's package manager." >&2
      fi
      ;;
    Darwin)
      echo "  Install with: xcode-select --install && brew install cmake" >&2
      ;;
    *)
      echo "  Install gcc, g++, and cmake using your platform's standard mechanism." >&2
      ;;
  esac
  echo >&2
  echo "Re-run ./install.sh after installing." >&2
  exit 1
fi

# ---- pixi -----------------------------------------------------------------
if ! command -v pixi >/dev/null 2>&1; then
  log "Installing pixi to ~/.pixi"
  curl -fsSL https://pixi.sh/install.sh | bash
  export PATH="$HOME/.pixi/bin:$PATH"
else
  log "pixi found: $(pixi --version)"
fi

# ---- pixi environments ----------------------------------------------------
log "Materializing pixi environments (default, run, plot, fastp-nfcore, trim-galore-rs)"
pixi install                              # default
pixi install --environment run
pixi install --environment plot
pixi install --environment fastp-nfcore
# trim-galore-rs is only built for linux-64 / linux-aarch64 / osx-arm64 on
# bioconda; on osx-64 the feature's `platforms` restriction makes the env
# unsolvable. Skip the install there with a warning rather than aborting.
if [[ "$(uname -s)-$(uname -m)" == "Darwin-x86_64" ]]; then
  log "Skipping pixi env 'trim-galore-rs' on osx-64 (no bioconda 2.x build)"
else
  pixi install --environment trim-galore-rs
fi

# ---- rust toolchain -------------------------------------------------------
RUST_DIR="$PIPELINE_DIR/.rust"
if [[ "$USE_ISOLATED_RUST" -eq 1 ]]; then
  export RUSTUP_HOME="$RUST_DIR/rustup"
  export CARGO_HOME="$RUST_DIR/cargo"
  export PATH="$CARGO_HOME/bin:$PATH"

  if [[ ! -x "$CARGO_HOME/bin/rustup" ]]; then
    log "Installing rustup into $RUST_DIR (host ~/.cargo and ~/.rustup are not touched)"
    mkdir -p "$RUST_DIR"
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
      | sh -s -- -y --no-modify-path --default-toolchain none --profile minimal
  else
    log "Isolated rustup present at $CARGO_HOME"
  fi

  # rust-toolchain.toml at the repo root pins channel = 1.95; the first
  # cargo invocation auto-installs it. Force it now so build progress is
  # visible up front rather than mixed with the chelae build.
  ( cd "$REPO_ROOT" && rustup show active-toolchain >/dev/null )
else
  if ! command -v cargo >/dev/null 2>&1; then
    echo "ERROR: --system-rust set but cargo not on PATH. Install rustup or drop --system-rust." >&2
    exit 1
  fi
  log "Using system cargo: $(cargo --version)"
fi

# ---- cargo-multivers (x86_64 only) ----------------------------------------
ARCH="$(uname -m)"
if [[ "$ARCH" == "x86_64" || "$ARCH" == "amd64" ]]; then
  if ! command -v cargo-multivers >/dev/null 2>&1; then
    log "Installing cargo-multivers"
    # Floor of 0.12.0: that release carries the fexecve/memfd fix required for
    # the launcher to run (and preserve argv[0]) under binfmt emulation.
    cargo install cargo-multivers --locked --version '>=0.12.0'
  else
    log "cargo-multivers found: $(cargo multivers --version 2>/dev/null || echo 'present')"
  fi
fi

# ---- holodeck (read simulator, from crates.io) ----------------------------
# Pinned, all-platform, sidesteps the bioconda linux-64 gap on 0.2.0.
# `--root` forces the install location regardless of CARGO_HOME so config.yaml
# can rely on `.rust/cargo/bin/holodeck` whether we're using isolated or
# system rust. cargo install is a no-op when the version already matches.
HOLODECK_VERSION="0.2.1"
HOLODECK_ROOT="$PIPELINE_DIR/.rust/cargo"
if [[ -x "$HOLODECK_ROOT/bin/holodeck" ]] \
   && "$HOLODECK_ROOT/bin/holodeck" --version 2>/dev/null | grep -qF "$HOLODECK_VERSION"; then
  log "holodeck $HOLODECK_VERSION already installed at $HOLODECK_ROOT/bin/holodeck"
else
  log "Installing holodeck $HOLODECK_VERSION via cargo into $HOLODECK_ROOT"
  cargo install --root "$HOLODECK_ROOT" --version "$HOLODECK_VERSION" holodeck
fi

# ---- build chelae ---------------------------------------------------------
if [[ "$DO_BUILD" -eq 1 ]]; then
  cd "$REPO_ROOT"
  case "$ARCH" in
    x86_64|amd64)
      log "Building chelae with cargo-multivers (x86-64 v1/v2/v4 variants)"
      cargo multivers --profile dist
      # cargo-multivers writes the launcher under target/<target-triple>/<profile>/.
      # Symlink to a stable, arch-independent path so config.yaml can use
      # `../target/dist/chelae` regardless of host.
      MV_LAUNCHER="$(find target -type f -path '*/dist/chelae' \
                       -not -path '*/deps/*' -not -path '*/build/*' \
                       -not -path 'target/dist/chelae' \
                       -print -quit)"
      if [[ -z "${MV_LAUNCHER:-}" ]]; then
        echo "ERROR: cargo multivers reported success but no launcher binary was found under target/" >&2
        exit 1
      fi
      mkdir -p target/dist
      ln -sf "../../$MV_LAUNCHER" target/dist/chelae
      ;;
    aarch64|arm64)
      log "Building chelae with cargo build --profile dist (single aarch64 binary)"
      cargo build --profile dist
      ;;
    *)
      echo "ERROR: unsupported architecture: $ARCH" >&2
      exit 1
      ;;
  esac
  cd "$PIPELINE_DIR"
  log "chelae binary: $(realpath ../target/dist/chelae)"
  ../target/dist/chelae --version || true
fi

cat <<EOF

Setup complete.

Next:
  1. Edit config/accuracy.config.yaml: point 'reference:' at an indexed FASTA
     (or leave the URL form to fetch one on first run).
  2. Run the smoke first to validate the pipeline (~15-20 min):
       ./run.sh config/smoke.config.yaml
  3. Then the bigger sweeps:
       ./run.sh                                  # accuracy (default; ~7-8 hr)
       ./run.sh config/performance.config.yaml   # performance (~8-12 hr)
     ./run.sh --dry-run                          # preview the job graph
EOF
