"""fastp render. Shared across `fastp` (latest bioconda) and `fastp-nfcore`
(pinned 1.1.0); the CLI surface for the flags used here is identical.

Design choices worth calling out:
  - `--detect_adapter_for_pe` is ON in PE mode: this is fastp's "headline"
    mode — overlap-based adapter detection — and matches how nf-core and
    almost all fastp users run it in practice. We want the benchmark to
    reflect the flag set people actually use.
  - `--disable_quality_filtering` is passed whenever neither low-quality
    filtering nor N-base filtering is asked for. When only N-filtering is
    asked (e.g. `wgs`: filter_n_bases=true, filter_low_quality=false), we
    keep fastp's quality filter enabled (fastp couples n_base_limit to it)
    but pin Q thresholds to `--qualified_quality_phred 0
    --unqualified_percent_limit 100` so no read is dropped for quality
    reasons — leaving N-count as the sole per-read filter.
  - fastp's other optional features remain off by default: `--dedup`,
    `--correction`, `--overrepresentation_analysis`,
    `--low_complexity_filter`, `--umi / --umi_loc`.
  - The unified `min_adapter_overlap` config value is intentionally NOT
    threaded into fastp: fastp uses an internal fixed threshold for
    sequence-based adapter matching with no CLI flag exposed. Cross-tool
    accuracy comparisons should account for this asterisk.
"""


def render(ctx: dict) -> dict:
    cfg = ctx["trim_cfg"]
    workdir = ctx["workdir"]

    argv = ["fastp",
            "-i", ctx["input_r1"],
            "-o", ctx["output_r1"],
            "--thread", str(ctx["threads"]),
            "--compression", str(ctx["compression_level"]),
            "--adapter_sequence", ctx["adapter_r1"],
            # fastp writes report files unless we redirect them into the workdir.
            "--json", f"{workdir}/fastp.json",
            "--html", f"{workdir}/fastp.html"]
    if ctx["paired"]:
        argv += ["-I", ctx["input_r2"], "-O", ctx["output_r2"],
                 "--adapter_sequence_r2", ctx["adapter_r2"],
                 "--detect_adapter_for_pe"]

    # Quality + N-base filtering. fastp couples them: `--n_base_limit` only
    # takes effect when quality filtering is enabled, so we keep QF on if
    # N-filtering is requested and neutralize its Q thresholds.
    filter_qual = cfg.get("filter_low_quality", False)
    filter_n = cfg.get("filter_n_bases", False)
    if not filter_qual and not filter_n:
        argv += ["--disable_quality_filtering"]
    elif not filter_qual and filter_n:
        argv += ["--qualified_quality_phred", "0",
                 "--unqualified_percent_limit", "100"]
    # Length filtering
    if cfg.get("min_length", 0) > 0:
        argv += ["--length_required", str(cfg["min_length"])]
    else:
        argv += ["--disable_length_filtering"]
    # polyG
    if not cfg.get("polyg_trim"):
        argv += ["--disable_trim_poly_g"]
    # polyX
    if cfg.get("polyx_trim"):
        argv += ["--trim_poly_x"]
    # Sliding-window quality trim on tail
    if cfg.get("quality_trim"):
        argv += ["--cut_tail",
                 "--cut_tail_window_size", str(cfg.get("quality_window", 4)),
                 "--cut_tail_mean_quality", str(cfg.get("quality_threshold", 20))]
    # N base filter
    if filter_n:
        argv += ["--n_base_limit", str(cfg.get("max_n_bases", 5))]

    return {"argv": argv, "moves": {}}
