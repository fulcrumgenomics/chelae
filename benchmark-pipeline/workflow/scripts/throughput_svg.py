#!/usr/bin/env python3
"""Draws the README's throughput chart from a tier-2 bench.tsv: one horizontal bar
per tool, the mean across samples of each sample's read pairs per second at its
median wall time. chelae is drawn in the accent colour and the rest in grey.
Writes a light and a dark variant for a <picture> element to choose between.

usage: throughput_svg.py <bench.tsv> <out_dir>
       (writes <out_dir>/throughput-light.svg and <out_dir>/throughput-dark.svg)"""

import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

THEMES = {
    "light": {"primary": "#0b0b0b", "secondary": "#52514e", "axis": "#c3c2b7",
              "accent": "#2a78d6", "other": "#898781"},
    "dark": {"primary": "#ffffff", "secondary": "#c3c2b7", "axis": "#383835",
             "accent": "#3987e5", "other": "#898781"},
}
FONT = 'system-ui, -apple-system, "Segoe UI", Helvetica, Arial, sans-serif'
WIDTH = 720
NAME_COL = 150    # x of the bars' baseline; tool names are right-aligned just left of it
VALUE_ROOM = 70   # space kept right of the longest bar for its value label
TOP = 72          # y of the first bar row
ROW = 40
BAR = 22          # bar thickness; at most 24 px
RADIUS = 4        # rounded data end; the baseline end stays square
NUMBER_WORDS = {1: "one", 2: "two", 3: "three", 4: "four"}


def throughput(bench_tsv):
    """Returns [(tool, mean M pairs/s)] sorted fastest first, and the samples used."""
    walls = defaultdict(list)
    pairs = {}
    with open(bench_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            walls[(r["tool"], r["sample"])].append(float(r["wall_s"]))
            pairs[r["sample"]] = int(r["reads_in"]) / 2
    samples = sorted(pairs)
    per_tool = defaultdict(list)
    for (tool, sample), w in walls.items():
        per_tool[tool].append(pairs[sample] / statistics.median(w) / 1e6)
    for tool, rates in per_tool.items():
        if len(rates) != len(samples):
            sys.exit(f"{tool} is missing runs for some of {samples}")
    means = [(tool, statistics.mean(rates)) for tool, rates in per_tool.items()]
    return sorted(means, key=lambda tm: -tm[1]), samples


def bar_path(x, y, length, height, r):
    """A bar from the baseline at `x`, square there and rounded at its far end."""
    r = min(r, length, height / 2)
    return (f"M{x},{y} h{length - r:.1f} a{r},{r} 0 0 1 {r},{r} v{height - 2 * r} "
            f"a{r},{r} 0 0 1 {-r},{r} h{-(length - r):.1f} z")


def svg(rows, samples, colours):
    height = TOP + ROW * len(rows) + 8
    scale = (WIDTH - NAME_COL - VALUE_ROOM) / rows[0][1]
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{height}" '
        f'viewBox="0 0 {WIDTH} {height}" font-family=\'{FONT}\' role="img" '
        f'aria-labelledby="title desc">',
        '<title id="title">Throughput at 8 threads</title>',
        '<desc id="desc">' + "; ".join(f"{t}: {v:.2f} million read pairs per second" for t, v in rows)
        + "</desc>",
        f'<text x="0" y="22" font-size="16" font-weight="600" fill="{colours["primary"]}">'
        "Throughput at 8 threads</text>",
        f'<text x="0" y="44" font-size="13" fill="{colours["secondary"]}">Million read pairs per '
        f"second, mean of {NUMBER_WORDS.get(len(samples), len(samples))} paired-end datasets of "
        "50 M pairs each</text>",
    ]
    for i, (tool, value) in enumerate(rows):
        y = TOP + i * ROW
        mid = y + BAR / 2
        length = value * scale
        is_chelae = tool == "chelae"
        weight = ' font-weight="600"' if is_chelae else ""
        ink = colours["primary"] if is_chelae else colours["secondary"]
        fill = colours["accent"] if is_chelae else colours["other"]
        out.append(f'<text x="{NAME_COL - 10}" y="{mid}" dy="0.35em" text-anchor="end" '
                   f'font-size="14"{weight} fill="{ink}">{tool}</text>')
        out.append(f'<path d="{bar_path(NAME_COL, y, length, BAR, RADIUS)}" fill="{fill}"/>')
        out.append(f'<text x="{NAME_COL + length + 8:.1f}" y="{mid}" dy="0.35em" font-size="14"'
                   f'{weight} fill="{ink}">{value:.2f}</text>')
    axis_bottom = TOP + ROW * (len(rows) - 1) + BAR
    out.append(f'<line x1="{NAME_COL}" y1="{TOP - 6}" x2="{NAME_COL}" y2="{axis_bottom + 6}" '
               f'stroke="{colours["axis"]}" stroke-width="1"/>')
    out.append("</svg>")
    return "\n".join(out) + "\n"


def main():
    bench_tsv, out_dir = sys.argv[1], Path(sys.argv[2])
    rows, samples = throughput(bench_tsv)
    for tool, value in rows:
        print(f"{tool}\t{value:.3f}")
    out_dir.mkdir(parents=True, exist_ok=True)
    for theme, colours in THEMES.items():
        (out_dir / f"throughput-{theme}.svg").write_text(svg(rows, samples, colours))


if __name__ == "__main__":
    main()
