"""AdapterRemoval 3 render (`adapterremoval3`).

Several 3.x defaults would do work the semantic trim config didn't ask for,
so each is set explicitly: quality trimming (default `mott`), poly-X
pre-trimming (default `auto`, which trims G tails on 2-colour data), the N
filter (default `--max-ns-fraction 0.05`) and the length filter (default 15).

`--min-overlap` is only threaded through in single-end mode. In paired-end
mode it is the minimum mate-to-mate alignment length (default 11, shared with
`--merge-threshold`), not an adapter overlap, and lowering it admits spurious
short mate alignments.

Singletons are sent to /dev/null uncompressed (an `--out-*` path without a
`.gz` suffix selects plain output) so they cost no compression work, and the
JSON/HTML reports go to the workdir."""

from pathlib import Path


def render(ctx: dict) -> dict:
    cfg = ctx["trim_cfg"]
    workdir = Path(ctx["workdir"])

    argv = ["adapterremoval3",
            "--threads", str(ctx["threads"]),
            "--in-file1", ctx["input_r1"],
            "--out-file1", ctx["output_r1"],
            "--adapter1", ctx["adapter_r1"],
            "--compression-level", str(ctx["compression_level"]),
            # Reference contigs (decoys/HLAs) contain IUPAC codes that propagate
            # into simulated reads; AR aborts on anything outside ACGTN unless
            # told to mask them to N.
            "--mask-degenerate-bases",
            "--out-json", str(workdir / "adapterremoval.json"),
            "--out-html", str(workdir / "adapterremoval.html")]
    if ctx["paired"]:
        argv += ["--in-file2", ctx["input_r2"],
                 "--out-file2", ctx["output_r2"],
                 "--adapter2", ctx["adapter_r2"],
                 "--out-singleton", "/dev/null"]
    elif "min_adapter_overlap" in cfg:
        argv += ["--min-overlap", str(cfg["min_adapter_overlap"])]

    argv += ["--pre-trim-polyx", "G" if cfg.get("polyg_trim") else "off"]
    if cfg.get("polyx_trim"):
        argv += ["--post-trim-polyx"]

    if cfg.get("quality_trim"):
        argv += ["--quality-trimming", "window",
                 "--trim-windows", str(cfg.get("quality_window", 4)),
                 "--trim-min-quality", str(cfg.get("quality_threshold", 20)),
                 "--preserve5p"]
    else:
        argv += ["--quality-trimming", "none"]

    if cfg.get("filter_n_bases"):
        argv += ["--max-ns", str(cfg.get("max_n_bases", 5))]
    else:
        argv += ["--max-ns-fraction", "1"]
    argv += ["--min-length", str(cfg.get("min_length", 0))]

    return {"argv": argv, "moves": {}}
