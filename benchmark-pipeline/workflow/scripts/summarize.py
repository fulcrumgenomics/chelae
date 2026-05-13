#!/usr/bin/env python3
"""Roll up bench.tsv / accuracy.tsv into compact summary TSVs.

Two modes:
  --kind bench    — collapse replicates: median + min/max wall, median
                    throughput, max RSS per (sample, trim_config, tool, threads).
  --kind accuracy — derive headline accuracy stats per (sample, trim_config,
                    tool, threads, rep, mate) from the melted matrix.

The summaries are what plotting consumes and what humans skim. The full
melted TSVs (bench.tsv, accuracy.tsv) remain the source of truth.
"""

import argparse
import math
import sys
from pathlib import Path

import pandas as pd


def summarize_bench(in_tsv: Path, out_tsv: Path) -> None:
    df = pd.read_csv(in_tsv, sep="\t")
    if df.empty:
        df.to_csv(out_tsv, sep="\t", index=False)
        return

    group_cols = [
        "sample", "trim_config", "tool", "threads", "tool_version",
        "host_hostname", "host_arch", "host_cpu_model",
    ]
    # Some host fields may be missing on aggregated cross-host data; only
    # group by what's present.
    group_cols = [c for c in group_cols if c in df.columns]

    g = df.groupby(group_cols, dropna=False)
    summary = g.agg(
        replicates       = ("rep", "count"),
        wall_s_median    = ("wall_s", "median"),
        wall_s_min       = ("wall_s", "min"),
        wall_s_max       = ("wall_s", "max"),
        # IQR via numeric quantiles. groupby agg lambda is fine here —
        # group counts are tiny (≤ replicates count).
        wall_s_iqr       = ("wall_s", lambda s: s.quantile(0.75) - s.quantile(0.25)),
        user_s_median    = ("user_s", "median"),
        sys_s_median     = ("sys_s", "median"),
        reads_per_s      = ("reads_per_s", "median"),
        bases_per_s      = ("bases_per_s", "median"),
        max_rss_kb       = ("max_rss_kb", "max"),
        cpu_percent      = ("cpu_percent", "median"),
        reads_in         = ("reads_in", "first"),
        bases_in         = ("bases_in", "first"),
    ).reset_index()
    summary.to_csv(out_tsv, sep="\t", index=False)


def summarize_accuracy(in_tsv: Path, out_tsv: Path) -> None:
    df = pd.read_csv(in_tsv, sep="\t")

    out_cols = [
        "sample", "trim_config", "tool", "threads", "rep", "mate",
        "host_hostname", "host_arch", "host_cpu_model",
        "total_reads", "kept_reads", "dropped_reads",
        "exact_match", "exact_match_rate",
        "false_positive", "false_negative",
        "under_trimmed", "over_trimmed",
        "error_events", "avg_error_magnitude",
        "mean_abs_error", "rmse",
    ]

    if df.empty:
        pd.DataFrame(columns=out_cols).to_csv(out_tsv, sep="\t", index=False)
        return

    # Only adapter_only configs feed the accuracy eval; the producer enforces
    # this, but be defensive.
    df = df[df["dropped"].isin([0, 1])].copy()

    group_cols = ["sample", "trim_config", "tool", "threads", "rep", "mate",
                  "host_hostname", "host_arch", "host_cpu_model"]
    group_cols = [c for c in group_cols if c in df.columns]

    rows = []
    for keys, sub in df.groupby(group_cols, dropna=False):
        total = int(sub["count"].sum())
        dropped = int(sub.loc[sub["dropped"] == 1, "count"].sum())
        kept = sub[sub["dropped"] == 0]
        kept_total = int(kept["count"].sum())

        exact = kept[kept["expected_trim_len"] == kept["observed_trim_len"]]
        exact_n = int(exact["count"].sum())

        # False positive: read had no adapter (expected_trim_len == 0) but the
        # tool trimmed something. False negative: adapter present (expected > 0)
        # but the tool trimmed nothing.
        fp = int(kept[
            (kept["expected_trim_len"] == 0) & (kept["observed_trim_len"] > 0)
        ]["count"].sum())
        fn = int(kept[
            (kept["expected_trim_len"] > 0) & (kept["observed_trim_len"] == 0)
        ]["count"].sum())

        # Direction-of-error counts (excluding exact and the FP/FN edge cells).
        under = int(kept[
            (kept["expected_trim_len"] > kept["observed_trim_len"]) &
            ~((kept["expected_trim_len"] > 0) & (kept["observed_trim_len"] == 0))
        ]["count"].sum())
        over = int(kept[
            (kept["observed_trim_len"] > kept["expected_trim_len"]) &
            ~((kept["expected_trim_len"] == 0) & (kept["observed_trim_len"] > 0))
        ]["count"].sum())

        # Weighted MAE and RMSE over kept reads. These are "expected error
        # magnitude per read drawn at random" — counts are normalized by the
        # full kept_total denominator, so a tool with many tiny errors looks
        # similar to one with few large errors. `avg_error_magnitude` below
        # complements that: it's the average size of an error *given that one
        # occurred*.
        if kept_total > 0:
            diff = (kept["observed_trim_len"] - kept["expected_trim_len"]).abs()
            total_abs = float((diff * kept["count"]).sum())
            mae = total_abs / kept_total
            sq = (kept["observed_trim_len"] - kept["expected_trim_len"]).pow(2)
            rmse = math.sqrt((sq * kept["count"]).sum() / kept_total)
        else:
            total_abs = 0.0
            mae = float("nan")
            rmse = float("nan")

        # Count of reads where observed_trim_len != expected_trim_len. Note:
        # FP + FN + under + over equals this exactly (the four categories
        # partition the non-exact kept reads).
        error_events = fp + fn + under + over
        avg_err_mag = (total_abs / error_events) if error_events > 0 else 0.0

        rate = (exact_n / kept_total) if kept_total > 0 else float("nan")

        row = dict(zip(group_cols, keys if isinstance(keys, tuple) else (keys,)))
        row.update(
            total_reads        = total,
            kept_reads         = kept_total,
            dropped_reads      = dropped,
            exact_match        = exact_n,
            exact_match_rate   = rate,
            false_positive     = fp,
            false_negative     = fn,
            under_trimmed      = under,
            over_trimmed       = over,
            error_events       = error_events,
            avg_error_magnitude= avg_err_mag,
            mean_abs_error     = mae,
            rmse               = rmse,
        )
        rows.append(row)

    out_df = pd.DataFrame(rows, columns=out_cols)
    out_df.to_csv(out_tsv, sep="\t", index=False)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--kind", choices=["bench", "accuracy"], required=True)
    p.add_argument("--in", dest="in_tsv", required=True)
    p.add_argument("--out", dest="out_tsv", required=True)
    args = p.parse_args()

    in_tsv = Path(args.in_tsv)
    out_tsv = Path(args.out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    if args.kind == "bench":
        summarize_bench(in_tsv, out_tsv)
    else:
        summarize_accuracy(in_tsv, out_tsv)


if __name__ == "__main__":
    main()
