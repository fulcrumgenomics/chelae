"""Trim Galore (cutadapt wrapper) render. Trim Galore names outputs after the
input basename + `_val_1.fq.gz` / `_val_2.fq.gz`, so we redirect it to the
workdir and then move outputs to the canonical paths."""

from pathlib import Path


def render(ctx: dict) -> dict:
    cfg = ctx["trim_cfg"]
    workdir = Path(ctx["workdir"])

    # --stringency = trim_galore's name for the cutadapt `-O` minimum-overlap
    # parameter. Default in trim_galore is 1 (treats any 1bp adapter overlap
    # as a hit, flooding accuracy with spurious 1bp trims). We pass the
    # config-unified value so the comparison across tools is fair (see
    # config/trim_configs/adapter_only.yaml).
    stringency = str(cfg.get("min_adapter_overlap", 5))
    # trim_galore's docs warn that speed-up plateaus around 4 cores: each
    # `--cores N` actually spawns ~N+3 worker processes (cutadapt + read/write
    # pigz + fastqc). Passing 8 on an 8-core box oversubscribes and slows the
    # run vs. capping at 4. The benchmark passes ctx["threads"] = vCPU count;
    # we cap here so trim_galore is measured at its sweet-spot config.
    cores = min(int(ctx["threads"]), 4)
    argv = ["trim_galore",
            "--cores", str(cores),
            "--adapter", ctx["adapter_r1"],
            "--stringency", stringency,
            "--gzip",
            "--output_dir", str(workdir)]

    if ctx["paired"]:
        argv += ["--paired", "--adapter2", ctx["adapter_r2"]]

    # Quality trim is on by default (--quality 20). Disable by passing 0.
    argv += ["--quality", str(cfg.get("quality_threshold", 0)) if cfg.get("quality_trim") else "0"]
    # Min length (trim_galore default is 20; pass 1 to effectively disable).
    argv += ["--length", str(cfg["min_length"]) if cfg.get("min_length", 0) > 0 else "1"]

    argv.append(ctx["input_r1"])
    if ctx["paired"]:
        argv.append(ctx["input_r2"])

    # Compute Trim Galore's output names and map to our canonical paths.
    def basename_no_fqgz(path: str) -> str:
        n = Path(path).name
        for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
            if n.endswith(suffix):
                return n[: -len(suffix)]
        return Path(path).stem

    moves: dict[str, str] = {}
    if ctx["paired"]:
        b1 = basename_no_fqgz(ctx["input_r1"])
        b2 = basename_no_fqgz(ctx["input_r2"])
        moves[str(workdir / f"{b1}_val_1.fq.gz")] = ctx["output_r1"]
        moves[str(workdir / f"{b2}_val_2.fq.gz")] = ctx["output_r2"]
    else:
        b1 = basename_no_fqgz(ctx["input_r1"])
        moves[str(workdir / f"{b1}_trimmed.fq.gz")] = ctx["output_r1"]

    return {"argv": argv, "moves": moves}
