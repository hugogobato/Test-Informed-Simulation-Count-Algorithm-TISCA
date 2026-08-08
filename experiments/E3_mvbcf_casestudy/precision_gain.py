#!/usr/bin/env python3
"""Quantify the Monte Carlo precision gain from 500 to 1000 replications.

This is an ancillary diagnostic for P4-T1 / IJDA #12(c).  It reads the complete
1000-row blocks already collected by the E3 confirmatory shards; it does not fit
any model or alter the amended analysis, whose pilot rows (seeds 0..99) remain
discarded from final inference.

For each primary paired contrast, the script compares the estimate and MCSE from
the first ``J_LOW`` rows with those from all ``J_HIGH`` rows.  The expected MCSE
ratio under a stable contrast variance is ``sqrt(J_HIGH / J_LOW)`` and the
corresponding expected reduction is ``1 - sqrt(J_LOW / J_HIGH)``.
"""

from __future__ import annotations

import argparse
import math
import os

import numpy as np
import pandas as pd


_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
DEFAULT_RESULTS = os.path.join(_ROOT, "results", "E3")

CELLS = [(1, 500), (2, 500), (3, 500), (1, 100)]
CONTRASTS = [
    ("C1", "MVBCF vs BCF, PEHE Y1", "mvbcf_pehe1", "bcf_pehe1"),
    ("C2", "MVBCF vs BCF, PEHE Y2", "mvbcf_pehe2", "bcf_pehe2"),
    ("C3", "MVBCF vs BART, PEHE Y1", "mvbcf_pehe1", "bart_pehe1"),
    ("C4", "MVBCF vs BART, PEHE Y2", "mvbcf_pehe2", "bart_pehe2"),
    ("C5", "MVBCF vs MVBART, PEHE Y1", "mvbcf_pehe1", "mvbart_pehe1"),
    ("C6", "MVBCF vs MVBART, PEHE Y2", "mvbcf_pehe2", "mvbart_pehe2"),
]


def load_cell(dgp: int, n: int, results_dir: str, j_high: int) -> pd.DataFrame:
    """Load and validate the raw, complete confirmatory block."""
    path = os.path.join(results_dir, f"DGP{dgp}_n{n}_confirmatory_replications.csv")
    if not os.path.exists(path):
        raise SystemExit(f"missing {path}")
    frame = pd.read_csv(path).sort_values("seed").reset_index(drop=True)
    if len(frame) != j_high:
        raise SystemExit(f"{path}: found {len(frame)} rows, expected {j_high}")
    if frame["seed"].tolist() != list(range(j_high)):
        raise SystemExit(f"{path}: seeds are not a complete 0..{j_high - 1} block")
    if int((frame["converged_flag"] == 0).sum()) != 0:
        raise SystemExit(f"{path}: non-converged rows present")
    return frame


def summarize(values: np.ndarray, j: int, ci_z: float) -> dict[str, float]:
    block = values[:j]
    mcse = float(np.std(block, ddof=1) / math.sqrt(j))
    return {
        "J": j,
        "estimate": float(np.mean(block)),
        "sd": float(np.std(block, ddof=1)),
        "mcse": mcse,
        "ci95_halfwidth": ci_z * mcse,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS)
    parser.add_argument("--j-low", type=int, default=500)
    parser.add_argument("--j-high", type=int, default=1000)
    args = parser.parse_args(argv)
    if not 2 <= args.j_low < args.j_high:
        raise SystemExit("require 2 <= j_low < j_high")

    # This diagnostic reports a conventional 95% precision interval.  It is not
    # the family-adjusted inferential interval in paired_contrasts.csv.
    ci_z = 1.959963984540054
    expected_ratio = math.sqrt(args.j_high / args.j_low)
    expected_reduction = 100.0 * (1.0 - 1.0 / expected_ratio)
    rows = []

    for dgp, n in CELLS:
        frame = load_cell(dgp, n, args.results_dir, args.j_high)
        for contrast, label, mv_col, benchmark_col in CONTRASTS:
            differences = (
                frame[mv_col].to_numpy(float) - frame[benchmark_col].to_numpy(float)
            )
            low = summarize(differences, args.j_low, ci_z)
            high = summarize(differences, args.j_high, ci_z)
            rows.append(
                {
                    "dgp": dgp,
                    "n": n,
                    "contrast": contrast,
                    "label": label,
                    "J_low": args.j_low,
                    "J_high": args.j_high,
                    "estimate_low": low["estimate"],
                    "estimate_high": high["estimate"],
                    "estimate_change_high_minus_low": high["estimate"] - low["estimate"],
                    "sd_low": low["sd"],
                    "sd_high": high["sd"],
                    "mcse_low": low["mcse"],
                    "mcse_high": high["mcse"],
                    "mcse_low_over_high": low["mcse"] / high["mcse"],
                    "mcse_reduction_pct": 100.0 * (1.0 - high["mcse"] / low["mcse"]),
                    "ci95_halfwidth_low": low["ci95_halfwidth"],
                    "ci95_halfwidth_high": high["ci95_halfwidth"],
                    "ci95_halfwidth_reduction_pct": 100.0
                    * (1.0 - high["ci95_halfwidth"] / low["ci95_halfwidth"]),
                    "expected_mcse_low_over_high": expected_ratio,
                    "expected_mcse_reduction_pct": expected_reduction,
                }
            )

    output = pd.DataFrame(rows)
    os.makedirs(args.results_dir, exist_ok=True)
    csv_path = os.path.join(args.results_dir, "precision_gain.csv")
    output.to_csv(csv_path, index=False)

    grouped = output.groupby(["dgp", "n"], sort=False)
    summary_rows = []
    for (dgp, n), group in grouped:
        summary_rows.append(
            {
                "dgp": dgp,
                "n": n,
                "mean_mcse_low_over_high": group["mcse_low_over_high"].mean(),
                "min_mcse_reduction_pct": group["mcse_reduction_pct"].min(),
                "mean_mcse_reduction_pct": group["mcse_reduction_pct"].mean(),
                "max_mcse_reduction_pct": group["mcse_reduction_pct"].max(),
                "max_abs_estimate_change": group[
                    "estimate_change_high_minus_low"
                ].abs().max(),
            }
        )
    summary = pd.DataFrame(summary_rows)

    def f(value: float) -> str:
        return f"{value:.4f}"

    md_path = os.path.join(args.results_dir, "precision_gain.md")
    lines = [
        "# P4-T1 precision gain from 500 to 1000 replications",
        "",
        "This ancillary diagnostic reads the complete 1000-row E3 confirmatory "
        "blocks already collected. It performs no model fits and does not replace "
        "the amended final inference, which discards seeds 0..99 and uses seeds "
        "100..999.",
        "",
        f"For each of the six primary paired PEHE contrasts, `J = {args.j_low}` "
        f"uses the first {args.j_low} seed-sorted rows and `J = {args.j_high}` "
        f"uses all {args.j_high} rows. The expected MCSE ratio "
        f"(`MCSE_{args.j_low} / MCSE_{args.j_high}`) is "
        f"`sqrt({args.j_high}/{args.j_low}) = {f(expected_ratio)}`, "
        f"equivalent to a {f(expected_reduction)}% reduction in MCSE.",
        "",
        "| cell | mean MCSE ratio | MCSE reduction range | mean reduction | "
        "max absolute estimate change |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        lines.append(
            f"| DGP{row['dgp']} n={row['n']} | "
            f"{f(row['mean_mcse_low_over_high'])} | "
            f"{f(row['min_mcse_reduction_pct'])}%--"
            f"{f(row['max_mcse_reduction_pct'])}% | "
            f"{f(row['mean_mcse_reduction_pct'])}% | "
            f"{f(row['max_abs_estimate_change'])} |"
        )

    lines.extend(
        [
            "",
            "Across all 24 cell-by-contrast rows, the observed MCSE ratio ranges "
            f"from `{f(output['mcse_low_over_high'].min())}` to "
            f"`{f(output['mcse_low_over_high'].max())}`, with mean "
            f"`{f(output['mcse_low_over_high'].mean())}`. The observed MCSE "
            f"reduction ranges from `{f(output['mcse_reduction_pct'].min())}%` to "
            f"`{f(output['mcse_reduction_pct'].max())}%`, with mean "
            f"`{f(output['mcse_reduction_pct'].mean())}%`. This is consistent "
            "with the square-root precision law, up to finite-sample variation.",
            "",
            "The row-level details, including estimates, MCSEs, 95% precision "
            "half-widths, and the expected benchmark, are in `precision_gain.csv`.",
        ]
    )
    with open(md_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")

    print(f"wrote {csv_path} ({len(output)} rows)")
    print(f"wrote {md_path}")
    print(
        f"observed MCSE ratio range {output['mcse_low_over_high'].min():.4f}--"
        f"{output['mcse_low_over_high'].max():.4f}; "
        f"mean reduction {output['mcse_reduction_pct'].mean():.4f}%"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
