#!/usr/bin/env python3
"""Aggregate per-run outputs into a single wide TSV.

Two modes:
  --kind bench    — walk results/trim/**/time.txt, parse GNU time output,
                    join with host info + sample metadata + tool version.
                    Emits one row per (sample, trim_config, tool, threads, rep).
  --kind accuracy — walk results/eval/**/matrix.tsv, concat with provenance
                    columns. Emits one row per (run × expected × observed).

The merged TSVs are the sole artifact plotting and external analysis consume —
everything needed lives in one file per kind.
"""

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

# parse_gnu_time sits next to this script, not on the default path.
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd
import yaml

from parse_gnu_time import parse as parse_time


# Path schema:
#   results/trim/{sample}/{trim_config}/{tool}/t{threads}/rep{rep}/time.txt
#   results/eval/{sample}/{trim_config}/{tool}/t{threads}/rep{rep}/matrix.tsv
PATH_RE = re.compile(
    r"results/(?:trim|eval)/"
    r"(?P<sample>[^/]+)/"
    r"(?P<trim_config>[^/]+)/"
    r"(?P<tool>[^/]+)/"
    r"t(?P<threads>\d+)/"
    r"rep(?P<rep>\d+)/"
)


def parse_path(path: Path) -> dict:
    m = PATH_RE.search(str(path))
    if not m:
        raise ValueError(f"path doesn't match expected schema: {path}")
    d = m.groupdict()
    d["threads"] = int(d["threads"])
    d["rep"] = int(d["rep"])
    return d


def tool_version(tool: str, tools_cfg: dict, overrides: dict) -> str:
    """Run `version_cmd` for this tool, capture first non-empty line. Best-effort.

    For tools whose binary lives outside the default conda env:
      * `overrides[tool]` — full custom shell command (used for chelae, which
        lives in the workspace at ../target/release/chelae).
      * `tools_cfg[tool]["pixi_env"]` — run `pixi run -e <env> <version_cmd>`
        so we capture the version of the binary actually dispatched to, not
        whatever happens to be on PATH."""
    spec = tools_cfg.get(tool, {})
    if tool in overrides:
        cmd = overrides[tool]
    else:
        cmd = spec.get("version_cmd")
        env = spec.get("pixi_env")
        if cmd and env:
            cmd = f"pixi run -e {env} -- {cmd}"
    if not cmd:
        return spec.get("version", "unknown")
    try:
        out = subprocess.check_output(cmd, shell=True, text=True, timeout=10,
                                      stderr=subprocess.STDOUT)
        for line in out.splitlines():
            line = line.strip()
            if line:
                return line
    except Exception as e:
        return spec.get("version", f"error: {e}")
    return spec.get("version", "unknown")


def sample_stats(sample: str) -> tuple[int, int]:
    """Return (reads_in, bases_in) from the per-sample stats.txt."""
    stats_file = Path(f"results/sim/{sample}/stats.txt")
    reads = bases = 0
    with open(stats_file) as fh:
        for line in fh:
            k, v = line.strip().split("\t")
            if k == "reads":
                reads = int(v)
            elif k == "bases":
                bases = int(v)
    return reads, bases


def merge_bench(host: dict, root: Path, samples_csv: Path, tools_yaml: Path,
                out: Path, version_overrides: dict) -> None:
    samples_df = pd.read_csv(samples_csv, sep="\t").set_index("name")
    with open(tools_yaml) as fh:
        tools_cfg = yaml.safe_load(fh)

    version_cache: dict[str, str] = {}
    sample_cache: dict[str, tuple[int, int]] = {}
    rows: list[dict] = []

    for time_txt in sorted(root.glob("*/*/*/t*/rep*/time.txt")):
        meta = parse_path(time_txt)
        timings = parse_time(time_txt)

        tool = meta["tool"]
        version_cache.setdefault(tool, tool_version(tool, tools_cfg, version_overrides))
        sample_cache.setdefault(meta["sample"], sample_stats(meta["sample"]))
        reads_in, bases_in = sample_cache[meta["sample"]]

        wall = timings.get("wall_s") or 0.0
        row = {
            **meta,
            "tool_version": version_cache[tool],
            "reads_in": reads_in,
            "bases_in": bases_in,
            "reads_per_s": reads_in / wall if wall > 0 else 0.0,
            "bases_per_s": bases_in / wall if wall > 0 else 0.0,
            "wall_s": wall,
            "user_s": timings.get("user_s", 0.0),
            "sys_s": timings.get("sys_s", 0.0),
            "max_rss_kb": timings.get("max_rss_kb", 0),
            "cpu_percent": timings.get("cpu_percent", 0),
            "exit_status": timings.get("exit_status", 0),
            "host_hostname": host.get("hostname"),
            "host_arch": host.get("arch"),
            "host_cpu_model": host.get("cpu_model"),
            "host_cpu_count": host.get("cpu_count_logical"),
            "host_mem_bytes": host.get("total_mem_bytes"),
            "host_os": host.get("os"),
            "host_os_release": host.get("os_release"),
            "host_aws_instance_type": host.get("aws_instance_type"),
        }
        rows.append(row)

    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)


def merge_accuracy(host: dict, root: Path, out: Path) -> None:
    frames: list[pd.DataFrame] = []
    for matrix_tsv in sorted(root.glob("*/*/*/t*/rep*/matrix.tsv")):
        meta = parse_path(matrix_tsv)
        df = pd.read_csv(matrix_tsv, sep="\t", comment="#")
        if df.empty:
            continue
        for k, v in meta.items():
            df[k] = v
        df["host_hostname"] = host.get("hostname")
        df["host_arch"] = host.get("arch")
        df["host_cpu_model"] = host.get("cpu_model")
        frames.append(df)
    if frames:
        pd.concat(frames, ignore_index=True).to_csv(out, sep="\t", index=False)
    else:
        # Emit an empty but valid TSV with the expected columns so downstream
        # doesn't choke on missing files.
        pd.DataFrame(columns=[
            "mate", "expected_trim_len", "observed_trim_len", "dropped", "count",
            "sample", "trim_config", "tool", "threads", "rep",
            "host_hostname", "host_arch", "host_cpu_model",
        ]).to_csv(out, sep="\t", index=False)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--kind", choices=["bench", "accuracy"], required=True)
    p.add_argument("--host", required=True)
    p.add_argument("--root", required=True)
    p.add_argument("--samples", default=None)
    p.add_argument("--tools", default=None)
    p.add_argument("--out", required=True)
    p.add_argument("--version-override", action="append", default=[],
                   help="tool=cmd override for version_cmd (repeatable)")
    args = p.parse_args()

    overrides = {}
    for entry in args.version_override:
        k, _, v = entry.partition("=")
        overrides[k] = v

    with open(args.host) as fh:
        host = json.load(fh)

    root = Path(args.root)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    if args.kind == "bench":
        if not args.samples or not args.tools:
            raise SystemExit("--samples and --tools required for bench merge")
        merge_bench(host, root, Path(args.samples), Path(args.tools), out, overrides)
    else:
        merge_accuracy(host, root, out)


if __name__ == "__main__":
    main()
