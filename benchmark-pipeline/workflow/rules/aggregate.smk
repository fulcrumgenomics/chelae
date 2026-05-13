# Collects per-run metrics into:
#   bench.tsv             — one row per trim run (raw, includes every replicate)
#   accuracy.tsv          — one row per (run, expected_len, observed_len) bucket
#   bench_summary.tsv     — replicates collapsed; medians + min/max/IQR per cell
#   accuracy_summary.tsv  — headline accuracy stats per (run, mate)
#   plots/                — ggplot bundle from plot.R

# Every CPU rule below requests bench=1; trim_one holds bench=100 (the
# full pool). That makes trim_one exclusive — no other rule can start
# while a timed trim is running — without otherwise serializing the rest
# of the DAG.

rule host_info:
    output:
        "results/host.json"
    resources:
        bench = 1
    shell:
        "python workflow/scripts/host_info.py > {output}"


rule aggregate_bench:
    input:
        trims = all_trim_outputs(),
        host = "results/host.json",
    output:
        "results/bench.tsv"
    params:
        chelae_bin = config["chelae_bin"],
        samples_tsv = str(SAMPLES_TSV),
    resources:
        bench = 1
    shell:
        r"""
        python workflow/scripts/merge_results.py \
            --kind bench \
            --host {input.host} \
            --root results/trim \
            --samples {params.samples_tsv} \
            --tools config/tools.yaml \
            --version-override "chelae={params.chelae_bin} --version" \
            --out {output}
        """


rule aggregate_accuracy:
    input:
        evals = all_eval_outputs(),
        host = "results/host.json",
    output:
        "results/accuracy.tsv"
    resources:
        bench = 1
    shell:
        r"""
        python workflow/scripts/merge_results.py \
            --kind accuracy \
            --host {input.host} \
            --root results/eval \
            --out {output}
        """


rule summarize_bench:
    input:
        "results/bench.tsv"
    output:
        "results/bench_summary.tsv"
    resources:
        bench = 1
    shell:
        r"""
        python workflow/scripts/summarize.py \
            --kind bench --in {input} --out {output}
        """


rule summarize_accuracy:
    input:
        "results/accuracy.tsv"
    output:
        "results/accuracy_summary.tsv"
    resources:
        bench = 1
    shell:
        r"""
        python workflow/scripts/summarize.py \
            --kind accuracy --in {input} --out {output}
        """


# Plotting runs in the `plot` pixi env (R + tidyverse). PDF outputs —
# vector format, scales cleanly for figures + slides without re-rendering.
rule plots:
    input:
        bench    = "results/bench.tsv",
        accuracy = "results/accuracy.tsv",
    output:
        throughput = "results/plots/throughput_vs_threads.pdf",
        wall       = "results/plots/wall_vs_threads.pdf",
        speedup    = "results/plots/speedup_vs_threads.pdf",
        rss        = "results/plots/max_rss.pdf",
        accuracy   = "results/plots/accuracy_heatmap.pdf",
    resources:
        bench = 1
    shell:
        r"""
        mkdir -p results/plots
        pixi run -e plot Rscript workflow/scripts/plot.R \
            {input.bench} results/plots {input.accuracy}
        """
