# holodeck simulation — one rule per sample. Always writes r1.fastq.gz; writes
# r2.fastq.gz iff the sample is paired. Also emits stats.txt with read/base
# counts (read from the fastq itself so it reflects what downstream tools see).

rule simulate_sample:
    input:
        reference = REFERENCE_LOCAL,
    output:
        r1 = "results/sim/{sample}/r1.fastq.gz",
        r2 = "results/sim/{sample}/r2.fastq.gz",
        stats = "results/sim/{sample}/stats.txt",
    params:
        holodeck = config["holodeck_bin"],
        reference = REFERENCE_LOCAL,
        prefix = lambda wc: f"results/sim/{wc.sample}/reads",
        row = lambda wc: SAMPLES.loc[wc.sample].to_dict(),
        adapter = lambda wc: sample_adapter(wc.sample),
        single_end_flag = lambda wc: "" if is_paired(wc.sample) else "--single-end",
        paired = lambda wc: "1" if is_paired(wc.sample) else "0",
        # Per-sample holodeck seed (see Snakefile:sample_seed). Order-
        # independent (depends only on sample name) so adding/reordering
        # samples doesn't shift other samples' RNG streams.
        seed = lambda wc: sample_seed(wc.sample),
    threads: 8
    resources:
        bench = 1  # blocked while a trim_one holds the bench pool
    log: "results/sim/{sample}/holodeck.log"
    shell:
        r"""
        {params.holodeck} simulate \
            -r {params.reference} \
            -o {params.prefix} \
            -c {params.row[depth]} \
            -l {params.row[r1_len]} \
            -d {params.row[fragment_mean]} \
            -s {params.row[fragment_stddev]} \
            --adapter-r1 '{params.adapter[r1]}' \
            --adapter-r2 '{params.adapter[r2]}' \
            --min-error-rate {params.row[min_error_rate]} \
            --max-error-rate {params.row[max_error_rate]} \
            --seed {params.seed} \
            -t {threads} \
            {params.single_end_flag} > {log} 2>&1

        mv {params.prefix}.r1.fastq.gz {output.r1}
        if [ "{params.paired}" = "1" ]; then
            mv {params.prefix}.r2.fastq.gz {output.r2}
            python workflow/scripts/fastq_stats.py {output.r1} {output.r2} > {output.stats}
        else
            : > {output.r2}
            python workflow/scripts/fastq_stats.py {output.r1} > {output.stats}
        fi
        """
