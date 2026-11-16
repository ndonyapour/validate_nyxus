"""
Combine per-case comparison CSVs into one table: one row per feature, one status column per case.

Usage:
    python summarize_comparisons.py comparisons/*_comparison.csv [--out comparisons/summary.csv]

Each case's status column is named after the CSV (the part before "_comparison.csv").
"worst_status" is the least favourable status across cases, and "max_rel_diff_percent"
is the largest relative difference across cases.
"""

import argparse
import csv
from pathlib import Path

# Least to most favourable.
STATUS_ORDER = ["different", "missing", "close", "identical", "no Nyxus equivalent"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("comparison_csvs", nargs="+")
    parser.add_argument("--out", default="comparisons/summary.csv")
    args = parser.parse_args()

    cases = []
    features = {}  # feature -> {"nyxus_feature": ..., case: (status, rel_diff)}
    for path in args.comparison_csvs:
        case = Path(path).name.removesuffix(".csv").removesuffix("_comparison")
        cases.append(case)
        with open(path) as f:
            for r in csv.DictReader(f):
                entry = features.setdefault(r["pyradiomics_feature"], {"nyxus_feature": r["nyxus_feature"]})
                rel = float(r["rel_diff_percent"]) if r["rel_diff_percent"] else None
                entry[case] = (r["status"], rel)

    rows = []
    for feat, entry in features.items():
        results = [entry.get(c, ("missing", None)) for c in cases]
        worst = min((s for s, _ in results), key=STATUS_ORDER.index)
        rels = [rel for _, rel in results if rel is not None]
        max_rel = max(rels) if rels else ""
        rows.append([feat, entry["nyxus_feature"], worst, max_rel] + [s for s, _ in results])
    rows.sort(key=lambda r: STATUS_ORDER.index(r[2]))

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pyradiomics_feature", "nyxus_feature", "worst_status", "max_rel_diff_percent"] + cases)
        w.writerows(rows)

    print(f"\nSummary over {len(cases)} case(s), worst status per feature:")
    for status in STATUS_ORDER:
        n = sum(1 for r in rows if r[2] == status)
        if n:
            print(f"  {status}: {n}")
    for r in rows:
        if r[2] == "different":
            bad = [c for c, s in zip(cases, r[4:]) if s == "different"]
            print(f"  different: {r[0]:42s} max rel_diff={r[3]:.2f}%  in {', '.join(bad)}")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
