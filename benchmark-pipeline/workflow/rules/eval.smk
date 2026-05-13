# Per-run accuracy eval. Only instantiated for trim_configs that pass
# is_accuracy_eligible() — see aggregate.smk for the gate on which outputs
# flow into accuracy.tsv.

rule eval_one:
    input:
        sim_r1 = "results/sim/{sample}/r1.fastq.gz",
        sim_r2 = "results/sim/{sample}/r2.fastq.gz",
        trim_r1 = "results/trim/{sample}/{trim_config}/{tool}/t{nthreads}/rep{rep}/r1.fastq.gz",
        trim_r2 = "results/trim/{sample}/{trim_config}/{tool}/t{nthreads}/rep{rep}/r2.fastq.gz",
    output:
        matrix = "results/eval/{sample}/{trim_config}/{tool}/t{nthreads,\\d+}/rep{rep,\\d+}/matrix.tsv",
    # Run eval_one before queuing the next trim_one. eval is the only
    # consumer of the temp() trim FASTQs, so the sooner it runs the sooner
    # snakemake reclaims that ~30 GB. Without this, a wide DAG of pending
    # trim jobs could queue up dozens of trimmed FASTQs on disk before
    # any eval gets scheduled.
    priority: 100
    resources:
        bench = 1  # blocked while a trim_one holds the bench pool
    params:
        paired = lambda wc: is_paired(wc.sample),
        r1_len = lambda wc: int(SAMPLES.loc[wc.sample, "r1_len"]),
        r2_len = lambda wc: int(SAMPLES.loc[wc.sample, "r2_len"]),
    shell:
        r"""
        python workflow/scripts/eval_accuracy.py \
            --sim-r1 {input.sim_r1} \
            --sim-r2 {input.sim_r2} \
            --trim-r1 {input.trim_r1} \
            --trim-r2 {input.trim_r2} \
            --paired {params.paired} \
            --r1-len {params.r1_len} \
            --r2-len {params.r2_len} \
            --out {output.matrix}
        """
