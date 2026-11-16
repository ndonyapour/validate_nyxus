"""
Shared by every dataset script: command-line options, case selection, and comparing
Nyxus with PyRadiomics.

Two values agree when |nyx - pyr| <= atol + tol * max(|nyx|, |pyr|).
Status per feature:
    identical            agrees with tol = 1e-6 (only floating-point noise)
    close                agrees with tol = --rtol (default 0.01 = 1%)
    different            disagrees with tol = --rtol
    missing              one of the tools did not produce the feature
    no Nyxus equivalent  PyRadiomics feature that is not in the mapping
"""

import argparse
import csv
import random
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "results"

IDENTICAL_RTOL = 1e-6
# Least to most favourable; the summary reports the worst status across cases.
STATUS_ORDER = ["different", "missing", "close", "identical", "no Nyxus equivalent"]


def make_parser(description, out):
    """The options every dataset script takes; add dataset-specific ones to the result."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--n", type=int, help="number of random cases (default: all)")
    parser.add_argument("--seed", type=int, default=0, help="random seed for --n (default: 0)")
    parser.add_argument("--case", action="append", help="run only this case (repeatable)")
    parser.add_argument("--rtol", type=float, default=0.01, help="relative tolerance for 'close'")
    parser.add_argument("--atol", type=float, default=1e-9, help="absolute tolerance near zero")
    parser.add_argument("--out", type=Path, default=out,
                        help=f"output folder (default: {out.relative_to(ROOT)})")
    return parser


def select_cases(parser, args, available):
    """--case if given, else --n random cases, else all."""
    available = sorted(available)
    if args.case:
        unknown = sorted(set(args.case) - set(available))
        if unknown:
            parser.error(f"unknown case(s): {', '.join(unknown)}")
        cases = sorted(args.case)
    elif args.n is not None:
        if args.n > len(available):
            parser.error(f"--n {args.n} is larger than the {len(available)} available cases")
        cases = sorted(random.Random(args.seed).sample(available, args.n))
    else:
        cases = available
    print(f"Running {len(cases)} of {len(available)} cases")
    return cases


def to_value(v):
    """PyRadiomics values are 0-d numpy arrays, and some shape features are complex (x+0j)."""
    if isinstance(v, (np.ndarray, np.generic, int, float, complex)) and np.ndim(v) == 0:
        return float(np.real(v))
    return v


def write_pyradiomics(path, result):
    write_csv(path, ["feature", "value"], [[name, to_value(value)] for name, value in result.items()])


def read_pyradiomics(path):
    """Read 'original_<class>_<feature>' rows, returned as '<class>_<feature>'."""
    features = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            if r["feature"].startswith("original_"):
                features[r["feature"].removeprefix("original_")] = float(r["value"])
    return features


def read_nyxus(path, roi):
    """The numeric columns of the row for mask label `roi`, or None if there is no such row."""
    with open(path) as f:
        rows = [r for r in csv.DictReader(f) if str(r["ROI_label"]) == str(roi)]
    if not rows:
        return None
    out = {}
    for k, v in rows[0].items():
        try:
            out[k] = float(v)
        except (ValueError, TypeError):
            pass
    return out


def agrees(a, b, rtol, atol):
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b))


def compare_case(pyr, nyx, mapping, rtol, atol):
    rows = []
    for pyr_name, nyx_name in mapping.items():
        pv = pyr.get(pyr_name)
        nv = nyx.get(nyx_name)
        if pv is None or nv is None:
            rows.append([pyr_name, nyx_name, pv, nv, "", "missing"])
            continue
        if agrees(nv, pv, IDENTICAL_RTOL, atol):
            status = "identical"
        elif agrees(nv, pv, rtol, atol):
            status = "close"
        else:
            status = "different"
        scale = max(abs(nv), abs(pv))
        rel_diff = abs(nv - pv) / scale * 100 if scale else 0.0
        rows.append([pyr_name, nyx_name, pv, nv, rel_diff, status])
    for pyr_name, pv in pyr.items():
        if pyr_name not in mapping:
            rows.append([pyr_name, "", pv, "", "", "no Nyxus equivalent"])
    return rows


def write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def write_case(out_dir, case, roi, rows):
    write_csv(out_dir / f"{case}_label{roi}.csv",
              ["pyradiomics_feature", "nyxus_feature", "pyradiomics_value", "nyxus_value",
               "rel_diff_percent", "status"], rows)
    counts = {s: sum(1 for r in rows if r[5] == s) for s in STATUS_ORDER}
    print(f"  label {roi}: " + ", ".join(f"{s} {n}" for s, n in counts.items() if n))


def write_summary(out_dir, roi, results):
    """One row per feature with the worst status and largest difference across cases."""
    names = {r[0]: r[1] for rows in results.values() for r in rows}
    summary = []
    for feat, nyx_name in names.items():
        per_case = {case: next((r for r in rows if r[0] == feat), None) for case, rows in results.items()}
        statuses = [r[5] if r else "missing" for r in per_case.values()]
        worst = min(statuses, key=STATUS_ORDER.index)
        diffs = [r[4] for r in per_case.values() if r and r[4] != ""]
        summary.append([feat, nyx_name, worst, max(diffs) if diffs else ""] + statuses)
    summary.sort(key=lambda r: STATUS_ORDER.index(r[2]))
    summary_path = out_dir / f"summary_label{roi}.csv"
    write_csv(summary_path,
              ["pyradiomics_feature", "nyxus_feature", "worst_status", "max_rel_diff_percent"]
              + list(results), summary)

    print(f"\nSummary over {len(results)} case(s), label {roi}, worst status per feature:")
    for status in STATUS_ORDER:
        n = sum(1 for r in summary if r[2] == status)
        if n:
            print(f"  {status}: {n}")
    for r in summary:
        if r[2] == "different":
            bad = [case for case, s in zip(results, r[4:]) if s == "different"]
            print(f"  different: {r[0]:42s} max rel_diff {r[3]:6.2f}%  in {', '.join(bad)}")
    print(f"Wrote {summary_path}")
