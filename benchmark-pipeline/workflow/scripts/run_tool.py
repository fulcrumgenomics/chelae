#!/usr/bin/env python3
"""Runner shared across all tools. Loads tools/<tool>/render.py, asks it to
translate the semantic trim config into argv for that tool, runs the command
under GNU time (and optionally perf stat), and moves the tool's output to the
canonical path Snakemake expects."""

import argparse
import importlib.util
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import yaml


def load_render(tool: str, render_dir: str | None):
    """Load the render module. `render_dir` overrides the default lookup of
    tools/<tool>/render.py, so version variants (e.g. fastp-nfcore → fastp)
    can share a single render script."""
    subdir = render_dir or tool
    path = Path("tools") / subdir / "render.py"
    if not path.exists():
        sys.exit(f"No render.py for tool {tool!r} at {path}")
    spec = importlib.util.spec_from_file_location(f"render_{tool}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.render


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tool", required=True)
    p.add_argument("--trim-config", required=True)
    p.add_argument("--adapters", required=True)
    p.add_argument("--adapter-set", required=True)
    p.add_argument("--threads", type=int, required=True)
    p.add_argument("--paired", required=True)
    p.add_argument("--in-r1", required=True)
    p.add_argument("--in-r2", required=True)
    p.add_argument("--out-r1", required=True)
    p.add_argument("--out-r2", required=True)
    p.add_argument("--time-txt", required=True)
    p.add_argument("--cmdline", required=True)
    p.add_argument("--log", required=True)
    p.add_argument("--chelae-bin", required=True)
    p.add_argument("--compression-level", type=int, default=4)
    p.add_argument("--warmup-cache", default="1")
    p.add_argument("--perf-stat", default="0")
    p.add_argument("--pixi-env", default="",
                   help="If set, wrap the tool invocation in `pixi run -e <env>`.")
    p.add_argument("--render-dir", default="",
                   help="Override subdir under tools/ to load render.py from.")
    args = p.parse_args()

    paired = str(args.paired).lower() == "true"

    with open(args.trim_config) as fh:
        trim_cfg = yaml.safe_load(fh)
    with open(args.adapters) as fh:
        adapters = yaml.safe_load(fh)
    # trim_config can override the adapter set; otherwise the sample's set wins.
    adapter_key = trim_cfg.get("adapter_set") or args.adapter_set
    adapter = adapters[adapter_key]

    workdir = Path(args.out_r1).parent
    workdir.mkdir(parents=True, exist_ok=True)

    ctx = {
        "paired": paired,
        "threads": args.threads,
        "trim_cfg": trim_cfg,
        "adapter_r1": adapter["r1"],
        "adapter_r2": adapter["r2"],
        "input_r1": args.in_r1,
        "input_r2": args.in_r2 if paired else None,
        "output_r1": args.out_r1,
        "output_r2": args.out_r2 if paired else None,
        "workdir": str(workdir),
        "chelae_bin": args.chelae_bin,
        "compression_level": args.compression_level,
    }

    plan = load_render(args.tool, args.render_dir or None)(ctx)
    argv = plan["argv"]

    # Dispatch to an alternate pixi env when requested (e.g. fastp-nfcore).
    # We need an absolute path to find pixi.toml since snakemake runs us from
    # the benchmark-pipeline dir already — `pixi run` locates the manifest by
    # walking up from cwd, which is fine here.
    if args.pixi_env:
        argv = ["pixi", "run", "-e", args.pixi_env, "--"] + argv
    # `moves` lets a tool declare "I will write to path A; please move it to path B
    # after I'm done". Needed for tools like Trim Galore / AdapterRemoval that
    # name outputs based on input basename.
    moves = plan.get("moves", {})
    # `cleanup` lets a tool declare files it produces but we don't want to keep.
    # Used for trimmomatic's unpaired_r{1,2}.fastq.gz which aren't snakemake
    # outputs and would otherwise leak disk at scale (~30 GB at 30x WGS).
    cleanup = plan.get("cleanup", [])

    # Warm page cache on inputs. Cheap and removes cold-read variance without
    # requiring root for drop_caches.
    if args.warmup_cache == "1":
        with open(os.devnull, "wb") as devnull:
            inputs = [args.in_r1] + ([args.in_r2] if paired else [])
            for path in inputs:
                subprocess.run(["cat", path], stdout=devnull, check=True)

    gnu_time = shutil.which("time") or "/usr/bin/time"
    cmd = [gnu_time, "-v", "-o", args.time_txt] + argv

    if args.perf_stat == "1" and shutil.which("perf"):
        perf_out = str(Path(args.time_txt).with_name("perf.txt"))
        cmd = ["perf", "stat", "-e", "cycles,instructions,task-clock",
               "-o", perf_out, "--"] + cmd

    with open(args.cmdline, "w") as fh:
        fh.write(shlex.join(argv) + "\n")

    with open(args.log, "wb") as log_fh:
        rc = subprocess.call(cmd, stdout=log_fh, stderr=subprocess.STDOUT)
    if rc != 0:
        sys.exit(f"Tool {args.tool!r} exited with code {rc}; see {args.log}")

    for src, dst in moves.items():
        src_p, dst_p = Path(src), Path(dst)
        if src_p.resolve() == dst_p.resolve():
            continue
        if not src_p.exists():
            sys.exit(f"Tool {args.tool!r} did not produce expected {src}")
        dst_p.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src_p), str(dst_p))

    required = [args.out_r1] + ([args.out_r2] if paired else [])
    for path in required:
        if not Path(path).exists():
            sys.exit(f"Expected output {path} missing after {args.tool!r} run")
    # SE samples still need an empty r2 placeholder to satisfy snakemake outputs.
    if not paired:
        Path(args.out_r2).touch()

    for path in cleanup:
        Path(path).unlink(missing_ok=True)


if __name__ == "__main__":
    main()
