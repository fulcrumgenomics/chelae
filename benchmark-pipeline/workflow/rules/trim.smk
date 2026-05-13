# Generic trim rule: dispatches to tools/<tool>/render.py via a common runner.
# Every (sample, trim_config, tool, threads, rep) combination instantiates one
# `trim_one` job. The `bench=1` resource serializes them across the DAG.

rule trim_one:
    input:
        r1 = "results/sim/{sample}/r1.fastq.gz",
        r2 = "results/sim/{sample}/r2.fastq.gz",
        sim_stats = "results/sim/{sample}/stats.txt",
    output:
        # Trimmed FASTQs are marked temp() — at 30x WGS the per-tool output
        # is ~30 GB, and we have (samples × trim_configs × tools × threads ×
        # reps) of them. Snakemake deletes each FASTQ as soon as no rule
        # depends on it (i.e. after eval_one runs for adapter_only configs;
        # immediately for wgs configs that don't feed accuracy scoring).
        # Logs and timings are NOT temp — they're the inputs to bench.tsv.
        # Tradeoff: a downstream bug means re-running trim, not just
        # re-aggregating. Worth it for the order-of-magnitude disk savings.
        r1 = temp("results/trim/{sample}/{trim_config}/{tool}/t{nthreads,\\d+}/rep{rep,\\d+}/r1.fastq.gz"),
        r2 = temp("results/trim/{sample}/{trim_config}/{tool}/t{nthreads,\\d+}/rep{rep,\\d+}/r2.fastq.gz"),
        time_txt = "results/trim/{sample}/{trim_config}/{tool}/t{nthreads,\\d+}/rep{rep,\\d+}/time.txt",
        cmdline = "results/trim/{sample}/{trim_config}/{tool}/t{nthreads,\\d+}/rep{rep,\\d+}/cmdline.txt",
        log = "results/trim/{sample}/{trim_config}/{tool}/t{nthreads,\\d+}/rep{rep,\\d+}/tool.log",
    threads: lambda wc: int(wc.nthreads)
    resources:
        # Exclusive-lock pattern: trim_one requests the full `bench` pool
        # (set to 100 on the snakemake CLI by run.sh). Every other CPU-using
        # rule requests bench=1, so when a trim is running nothing else
        # can dispatch — keeping wall-time measurements free of contention
        # from holodeck, eval, aggregation, etc. Outside of trim windows,
        # up to 100 non-trim rules can run in parallel (capped in practice
        # by --cores).
        bench = 100
    params:
        paired = lambda wc: is_paired(wc.sample),
        adapter_set = lambda wc: SAMPLES.loc[wc.sample, "adapter_set"],
        chelae_bin = config["chelae_bin"],
        compression_level = config.get("compression_level", 4),
        warmup = lambda wc: "1" if config.get("warmup_cache", True) else "0",
        perf_stat = lambda wc: "1" if config.get("perf_stat", False) else "0",
        pixi_env = lambda wc: TOOLS.get(wc.tool, {}).get("pixi_env", ""),
        render_dir = lambda wc: TOOLS.get(wc.tool, {}).get("render_dir", ""),
    shell:
        r"""
        python workflow/scripts/run_tool.py \
            --tool {wildcards.tool} \
            --trim-config config/trim_configs/{wildcards.trim_config}.yaml \
            --adapters config/adapters.yaml \
            --adapter-set {params.adapter_set} \
            --threads {wildcards.nthreads} \
            --paired {params.paired} \
            --in-r1 {input.r1} \
            --in-r2 {input.r2} \
            --out-r1 {output.r1} \
            --out-r2 {output.r2} \
            --time-txt {output.time_txt} \
            --cmdline {output.cmdline} \
            --log {output.log} \
            --chelae-bin {params.chelae_bin} \
            --compression-level {params.compression_level} \
            --warmup-cache {params.warmup} \
            --perf-stat {params.perf_stat} \
            --pixi-env "{params.pixi_env}" \
            --render-dir "{params.render_dir}"
        """
