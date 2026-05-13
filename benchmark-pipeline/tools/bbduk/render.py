"""BBDuk render. BBDuk uses `in=`/`out=` key-value argv. polyG and polyX are
handled via `trimpolyg`/`trimpolya`; adapter trim via ref= adapter FASTA
with ktrim=r."""

from pathlib import Path


def render(ctx: dict) -> dict:
    cfg = ctx["trim_cfg"]
    workdir = Path(ctx["workdir"])

    adapter_fa = workdir / "adapters.fa"
    with open(adapter_fa, "w") as fh:
        fh.write(f">adapter1\n{ctx['adapter_r1']}\n")
        if ctx["paired"]:
            fh.write(f">adapter2\n{ctx['adapter_r2']}\n")

    # `ordered=t` forces bbduk to emit reads in input order under
    # multi-threading, which the lockstep accuracy eval relies on. `tbo=t`
    # and `tpe=t` enable PE overlap-based adapter detection + symmetric
    # trimming — enabled uniformly across tools that support it, to match
    # how these tools are used in practice. `minlength` is set explicitly
    # because bbduk's default of 10 would silently filter even when our
    # semantic config asks for none. Threads are clamped to a minimum of 2
    # to dodge an internal bbduk assertion ("List size mismatch") that
    # triggers in PE mode under `threads=1` — report bbduk at threads=1 as
    # "min 2" in the final write-up.
    threads = max(int(ctx["threads"]), 2)
    # `k` must be ≤ the shortest adapter sequence we're matching against,
    # or bbduk loads 0 kmers and asserts ("KMER OPERATION WAS CHOSEN BUT
    # NO KMERS WERE LOADED"). Hardcoding k=23 worked for TruSeq (33 bp)
    # but crashed on small-rna (21 bp). Cap dynamically and floor at the
    # adapter length; `mink=11` still provides shorter partial-match kmers.
    k = min(len(ctx["adapter_r1"]), 23)
    if ctx["paired"]:
        k = min(k, len(ctx["adapter_r2"]))
    # `qin=33` forces phred33 quality encoding; bbduk's auto-detection
    # misfires on miRNA-style reads where N-padded tails produce a mix of
    # quality codes that look bimodal between phred33 and phred64. holodeck
    # always emits phred33, so the override is always safe for our data.
    #
    # `mink` controls partial-end (short overlap) k-mer matching — this is
    # the closest analogue to other tools' min-adapter-overlap parameter.
    # bbduk's default is none (only `k` matches count); we set it from the
    # config so the cross-tool comparison is fair. mink is capped at the
    # adapter length minus 1 (bbduk asserts mink < k).
    mink = max(2, min(int(cfg.get("min_adapter_overlap", 11)), k - 1))
    argv = ["bbduk.sh",
            f"threads={threads}",
            "qin=33",
            f"in={ctx['input_r1']}",
            f"out={ctx['output_r1']}",
            f"ref={adapter_fa}",
            "ktrim=r", f"k={k}", f"mink={mink}", "hdist=1",
            "ordered=t",
            f"ziplevel={ctx['compression_level']}",
            "overwrite=t"]
    if ctx["paired"]:
        argv += [f"in2={ctx['input_r2']}", f"out2={ctx['output_r2']}", "tpe=t", "tbo=t"]

    if cfg.get("quality_trim"):
        argv += ["qtrim=r", f"trimq={cfg.get('quality_threshold', 20)}"]
    argv += [f"minlength={cfg.get('min_length', 0)}"]
    if cfg.get("polyg_trim"):
        argv += ["trimpolyg=10"]
    if cfg.get("polyx_trim"):
        argv += ["trimpolya=10"]
    if cfg.get("filter_n_bases"):
        argv += [f"maxns={cfg.get('max_n_bases', 5)}"]

    return {"argv": argv, "moves": {}}
