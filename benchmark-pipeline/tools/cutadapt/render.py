"""cutadapt render. cutadapt handles adapter trim + quality trim + length
filter natively. It does NOT have polyG/polyX/N-base-filter (those are left
unimplemented and noted in the log if the config asks for them)."""


def render(ctx: dict) -> dict:
    cfg = ctx["trim_cfg"]

    # `-O` (minimum overlap, default 3) is set from the config's unified
    # `min_adapter_overlap` so accuracy comparisons across tools are
    # apples-to-apples; see config/trim_configs/adapter_only.yaml for the
    # rationale. `-e` (mismatch rate, default 0.1) is left at default.
    argv = ["cutadapt",
            "-j", str(ctx["threads"]),
            "--compression-level", str(ctx["compression_level"]),
            "-a", ctx["adapter_r1"],
            "-o", ctx["output_r1"]]
    if ctx["paired"]:
        argv += ["-A", ctx["adapter_r2"], "-p", ctx["output_r2"]]
    if "min_adapter_overlap" in cfg:
        argv += ["-O", str(cfg["min_adapter_overlap"])]

    if cfg.get("quality_trim"):
        argv += ["-q", str(cfg.get("quality_threshold", 20))]
    if cfg.get("min_length", 0) > 0:
        argv += ["--minimum-length", str(cfg["min_length"])]

    # positional inputs last
    argv.append(ctx["input_r1"])
    if ctx["paired"]:
        argv.append(ctx["input_r2"])

    return {"argv": argv, "moves": {}}
