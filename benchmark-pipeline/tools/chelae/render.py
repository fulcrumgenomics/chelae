"""Translate the semantic trim config into a chelae argv list."""


def render(ctx: dict) -> dict:
    cfg = ctx["trim_cfg"]
    argv = [ctx["chelae_bin"], "trim",
            "-t", str(ctx["threads"]),
            "--compression-level", str(ctx["compression_level"]),
            "-i", ctx["input_r1"]]
    if ctx["paired"]:
        argv.append(ctx["input_r2"])
    argv += ["-o", ctx["output_r1"]]
    if ctx["paired"]:
        argv.append(ctx["output_r2"])

    # Adapters
    argv += ["-a", ctx["adapter_r1"]]
    if ctx["paired"]:
        argv.append(ctx["adapter_r2"])

    # Length filter. Default is 15; we pass the config value (0 ⇒ effectively
    # disable — any read with ≥0 bases passes).
    argv += ["--filter-length", str(cfg.get("min_length", 0))]

    # Unified min adapter overlap (see config/trim_configs/adapter_only.yaml).
    if "min_adapter_overlap" in cfg:
        argv += ["--adapter-min-length", str(cfg["min_adapter_overlap"])]

    # polyG: default on at run=10; disable with --trim-polyg 0.
    polyg = 10 if cfg.get("polyg_trim") else 0
    argv += ["--trim-polyg", str(polyg)]

    # polyX: off unless enabled.
    if cfg.get("polyx_trim"):
        argv += ["--trim-polyx", "10"]

    # Sliding-window quality trim on the 3' end.
    if cfg.get("quality_trim"):
        window = cfg.get("quality_window", 4)
        thresh = cfg.get("quality_threshold", 20)
        argv += ["--quality-trim-3p", f"{window}:{thresh}"]

    # N filter.
    if cfg.get("filter_n_bases"):
        argv += ["--filter-max-ns", str(cfg.get("max_n_bases", 5))]

    return {"argv": argv, "moves": {}}
