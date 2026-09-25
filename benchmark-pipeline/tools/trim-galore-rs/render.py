"""Trim Galore Oxidized Edition (Rust rewrite). CLI is a near-superset of
the Perl trim_galore: same `--paired`, `--adapter`, `--adapter2`,
`--stringency`, `--quality`, `--length`, `--output_dir`, `--cores`. Output
filenames also follow the Perl convention (`<basename>_val_1.fq.gz` for
paired, `<basename>_trimmed.fq.gz` for SE), and output compression mirrors
the (gzipped) input.

Poly-G trimming is auto-detected from the data by default, so it is forced
on or off to match the trim config."""

from pathlib import Path


def render(ctx: dict) -> dict:
    cfg = ctx["trim_cfg"]
    workdir = Path(ctx["workdir"])

    stringency = str(cfg.get("min_adapter_overlap", 5))
    argv = ["trim_galore",
            "--cores", str(ctx["threads"]),
            "--adapter", ctx["adapter_r1"],
            "--stringency", stringency,
            "--compression", str(ctx["compression_level"]),
            "--output_dir", str(workdir)]

    if ctx["paired"]:
        argv += ["--paired", "--adapter2", ctx["adapter_r2"]]

    argv += ["--quality", str(cfg.get("quality_threshold", 0)) if cfg.get("quality_trim") else "0"]
    argv += ["--length", str(cfg["min_length"]) if cfg.get("min_length", 0) > 0 else "1"]
    argv += ["--poly_g" if cfg.get("polyg_trim") else "--no_poly_g"]
    if cfg.get("polyx_trim"):
        argv += ["--poly_a"]
    if cfg.get("filter_n_bases"):
        argv += ["--max_n", str(cfg.get("max_n_bases", 5))]

    argv.append(ctx["input_r1"])
    if ctx["paired"]:
        argv.append(ctx["input_r2"])

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
