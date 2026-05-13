"""Trimmomatic render. Trimmomatic needs adapters in a FASTA file, which we
materialize into the workdir at render time. Output naming allows direct
placement at the canonical paths."""

from pathlib import Path


def render(ctx: dict) -> dict:
    cfg = ctx["trim_cfg"]
    workdir = Path(ctx["workdir"])

    adapter_fa = workdir / "adapters.fa"
    with open(adapter_fa, "w") as fh:
        fh.write(f">PrefixPE/1\n{ctx['adapter_r1']}\n")
        if ctx["paired"]:
            fh.write(f">PrefixPE/2\n{ctx['adapter_r2']}\n")

    # Trimmomatic requires unpaired output paths for PE mode. We write them
    # into workdir and ignore them downstream (they're not part of Snakemake
    # outputs — they just satisfy the tool).
    unpaired_r1 = workdir / "unpaired_r1.fastq.gz"
    unpaired_r2 = workdir / "unpaired_r2.fastq.gz"

    mode = "PE" if ctx["paired"] else "SE"
    # The bioconda trimmomatic wrapper hardcodes a 1 GB JVM heap by default,
    # which OOMs on larger inputs (e.g. 2x WGS at 250bp PE — observed crash
    # in GzipParallelCompressor allocating output blocks). Passing -Xmx<size>
    # before the PE/SE subcommand is recognized by the wrapper and routed to
    # the JVM in lieu of the default. 8 GB is generous and fits comfortably
    # alongside the 16 GB box budget since trim_one runs alone (bench=100).
    argv = ["trimmomatic", "-Xmx8g", mode, "-threads", str(ctx["threads"])]

    if ctx["paired"]:
        argv += [ctx["input_r1"], ctx["input_r2"],
                 ctx["output_r1"], str(unpaired_r1),
                 ctx["output_r2"], str(unpaired_r2)]
    else:
        argv += [ctx["input_r1"], ctx["output_r1"]]

    # Pipeline steps (applied in order, as Trimmomatic requires).
    # ILLUMINACLIP:<fa>:<seedMismatches>:<palindromeClipThreshold>:<simpleClipThreshold>:<minAdapterLength>:<keepBothReads>
    #
    # We picked `2:2:7:1:true` after a parameter sweep against pe_cfdna_150
    # (10 configs, results in /tmp/tune_results.log). Findings:
    #
    #   - keepBothReads=true is mandatory: without it, palindrome detection
    #     discards R2 of every readthrough pair (treating it as "redundant
    #     with R1"), producing a 12% drop-rate of paired output on cfDNA.
    #   - minAdapterLength dominates the FN curve. At default 8, trimmomatic
    #     refuses to trim palindrome detections shorter than 8 bp; bringing
    #     it to 1 cut error count from 195k → 133k. Lower is monotonically
    #     better for accuracy on this geometry.
    #   - palindromeClipThreshold is the second strongest lever: 30 → 133k,
    #     15 → 100k, 5 → 69k, 2 → 61k, 1 → 58k. Diminishing returns past 2.
    #   - simpleClipThreshold is irrelevant here — palindrome mode catches
    #     virtually all reads when insert < 2× read_len. Set to 7 anyway
    #     for the rare insert > 300 bp case.
    #   - seedMismatches at 1 vs 2: indistinguishable.
    #
    # `2:2:7:1:true` is near the algorithmic floor (~61k errors). That's
    # still ~100× worse than fastp / bbduk / chelae on the same data —
    # trimmomatic's palindrome-only PE-overlap detection is fundamentally
    # less precise than evidence-combining detectors. nf-core/eager uses
    # `:5:true` (~170k errors here); we go more aggressive because we're
    # benchmarking the tool's best achievable accuracy, not its conventional
    # deployment.
    #
    # simpleClipThreshold is a log-likelihood score — each matched base
    # contributes ≈ 0.6 — so we approximate the unified `min_adapter_overlap`
    # by `round(0.6 × n)`. Not bp-exact but the closest cross-tool unification
    # available; preserves the relative scaling.
    simple_clip = max(2, round(0.6 * cfg.get("min_adapter_overlap", 12)))
    argv.append(f"ILLUMINACLIP:{adapter_fa}:2:2:{simple_clip}:1:true")

    if cfg.get("quality_trim"):
        w = cfg.get("quality_window", 4)
        q = cfg.get("quality_threshold", 20)
        argv.append(f"SLIDINGWINDOW:{w}:{q}")
    if cfg.get("min_length", 0) > 0:
        argv.append(f"MINLEN:{cfg['min_length']}")

    # Unpaired files are an artifact of PE mode; we don't keep them. Listing
    # them under `cleanup` has run_tool.py delete them after the tool runs,
    # matching the temp() lifecycle of the paired outputs.
    cleanup = [str(unpaired_r1), str(unpaired_r2)] if ctx["paired"] else []

    return {"argv": argv, "moves": {}, "cleanup": cleanup}
