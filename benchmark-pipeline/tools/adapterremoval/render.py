"""AdapterRemoval v2 render. Uses --file1/--file2 input, --output1/--output2
output. Quality trim is --trimqualities --minquality; length filter is
--minlength.

By default AdapterRemoval also writes `.discarded`, `.singleton.truncated`,
and a `.settings` summary. We redirect the first two to /dev/null (the bytes
would otherwise count as extra gzip work) and point `--settings` at the
workdir so it's captured but not noise in the benchmark surface.
`--minlength` is set explicitly because the default (15) would silently
filter even when our config asks for none."""

from pathlib import Path


def render(ctx: dict) -> dict:
    cfg = ctx["trim_cfg"]
    workdir = Path(ctx["workdir"])

    argv = ["AdapterRemoval",
            "--threads", str(ctx["threads"]),
            "--file1", ctx["input_r1"],
            "--output1", ctx["output_r1"],
            "--adapter1", ctx["adapter_r1"],
            "--gzip",
            "--gzip-level", str(ctx["compression_level"]),
            # Reference contigs (decoys/HLAs) contain IUPAC codes that propagate
            # into simulated reads; AR aborts on anything outside ACGTN unless
            # told to mask them to N.
            "--mask-degenerate-bases",
            "--settings", str(workdir / "adapterremoval.settings"),
            "--discarded", "/dev/null"]
    if ctx["paired"]:
        argv += ["--file2", ctx["input_r2"],
                 "--output2", ctx["output_r2"],
                 "--adapter2", ctx["adapter_r2"],
                 "--singleton", "/dev/null"]

    if cfg.get("quality_trim"):
        argv += ["--trimqualities",
                 "--minquality", str(cfg.get("quality_threshold", 20))]
    argv += ["--minlength", str(cfg.get("min_length", 0))]
    # AdapterRemoval default `--minadapteroverlap` is 0 — every nonzero
    # overlap counts as a hit, which floods short-readthrough accuracy with
    # 1-2bp false trims. Pass the unified value (see config/trim_configs/).
    if "min_adapter_overlap" in cfg:
        argv += ["--minadapteroverlap", str(cfg["min_adapter_overlap"])]

    return {"argv": argv, "moves": {}}
