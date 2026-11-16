"""
Compare Nyxus and PyRadiomics features for each case and summarize across cases.

Usage:
    python src/compare.py                  # all cases, label 1
    python src/compare.py --roi 2          # all cases, label 2
    python src/compare.py --case PE9f4fd4  # one case
    python src/compare.py --rtol 0.001     # stricter tolerance for "close"

Inputs : nyxus_out/<case>_nyxus.csv and pyradiomics_out/<case>_label<roi>.csv
Outputs: comparisons/<case>_label<roi>.csv  one row per feature
         comparisons/summary_label<roi>.csv one row per feature, one status column per case

Two values agree when |nyx - pyr| <= atol + tol * max(|nyx|, |pyr|).
Status per feature:
    identical            agrees with tol = 1e-6 (only floating-point noise)
    close                agrees with tol = --rtol (default 0.01 = 1%)
    different            disagrees with tol = --rtol
    missing              one of the tools did not produce the feature
    no Nyxus equivalent  PyRadiomics feature that is not in MAPPING

No conversions are applied: both tools are configured identically (see
nyxus_features.py and settings_hu1.yaml), so shape features are in mm in both.
"""

import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INTENSITY_DIR = ROOT / "CT_volumes" / "intensity"
NYXUS_DIR = ROOT / "nyxus_out"
PYRADIOMICS_DIR = ROOT / "pyradiomics_out"
OUT_DIR = ROOT / "comparisons"

IDENTICAL_RTOL = 1e-6
# Least to most favourable; the summary reports the worst status across cases.
STATUS_ORDER = ["different", "missing", "close", "identical", "no Nyxus equivalent"]

# PyRadiomics feature -> Nyxus column
# GLCM and GLRLM use the _AVE columns (PyRadiomics averages per-angle values).
MAPPING = {
    # shape
    "shape_VoxelVolume": "3VOXEL_VOLUME",
    "shape_MajorAxisLength": "3MAJOR_AXIS_LEN",
    "shape_MinorAxisLength": "3MINOR_AXIS_LEN",
    "shape_LeastAxisLength": "3LEAST_AXIS_LEN",
    "shape_Elongation": "3ELONGATION",
    "shape_Flatness": "3FLATNESS",
    "shape_MeshVolume": "3MESH_VOLUME",
    "shape_SurfaceArea": "3AREA",
    "shape_Sphericity": "3SPHERICITY",
    "shape_SurfaceVolumeRatio": "3AREA_2_VOLUME",
    # first order
    "firstorder_10Percentile": "3P10",
    "firstorder_90Percentile": "3P90",
    "firstorder_Energy": "3ENERGY",
    "firstorder_Entropy": "3ENTROPY",
    "firstorder_InterquartileRange": "3INTERQUARTILE_RANGE",
    "firstorder_Kurtosis": "3KURTOSIS",
    "firstorder_Maximum": "3MAX",
    "firstorder_MeanAbsoluteDeviation": "3MEAN_ABSOLUTE_DEVIATION",
    "firstorder_Mean": "3MEAN",
    "firstorder_Median": "3MEDIAN",
    "firstorder_Minimum": "3MIN",
    "firstorder_Range": "3RANGE",
    "firstorder_RobustMeanAbsoluteDeviation": "3ROBUST_MEAN_ABSOLUTE_DEVIATION",
    "firstorder_RootMeanSquared": "3ROOT_MEAN_SQUARED",
    "firstorder_Skewness": "3SKEWNESS",
    "firstorder_Uniformity": "3UNIFORMITY",
    "firstorder_Variance": "3VARIANCE",
    # glcm
    "glcm_Autocorrelation": "3GLCM_ACOR_AVE",
    "glcm_ClusterProminence": "3GLCM_CLUPROM_AVE",
    "glcm_ClusterShade": "3GLCM_CLUSHADE_AVE",
    "glcm_ClusterTendency": "3GLCM_CLUTEND_AVE",
    "glcm_Contrast": "3GLCM_CONTRAST_AVE",
    "glcm_Correlation": "3GLCM_CORRELATION_AVE",
    "glcm_DifferenceAverage": "3GLCM_DIFAVE_AVE",
    "glcm_DifferenceEntropy": "3GLCM_DIFENTRO_AVE",
    "glcm_DifferenceVariance": "3GLCM_DIFVAR_AVE",
    "glcm_Id": "3GLCM_ID_AVE",
    "glcm_Idm": "3GLCM_IDM_AVE",
    "glcm_Idmn": "3GLCM_IDMN_AVE",
    "glcm_Idn": "3GLCM_IDN_AVE",
    "glcm_Imc1": "3GLCM_INFOMEAS1_AVE",
    "glcm_Imc2": "3GLCM_INFOMEAS2_AVE",
    "glcm_InverseVariance": "3GLCM_IV_AVE",
    "glcm_JointAverage": "3GLCM_JAVE_AVE",
    "glcm_JointEnergy": "3GLCM_ASM_AVE",
    "glcm_JointEntropy": "3GLCM_JE_AVE",
    "glcm_MaximumProbability": "3GLCM_JMAX_AVE",
    "glcm_SumAverage": "3GLCM_SUMAVERAGE_AVE",
    "glcm_SumEntropy": "3GLCM_SUMENTROPY_AVE",
    "glcm_SumSquares": "3GLCM_VARIANCE_AVE",
    # glrlm
    "glrlm_GrayLevelNonUniformity": "3GLRLM_GLN_AVE",
    "glrlm_GrayLevelNonUniformityNormalized": "3GLRLM_GLNN_AVE",
    "glrlm_GrayLevelVariance": "3GLRLM_GLV_AVE",
    "glrlm_HighGrayLevelRunEmphasis": "3GLRLM_HGLRE_AVE",
    "glrlm_LongRunEmphasis": "3GLRLM_LRE_AVE",
    "glrlm_LongRunHighGrayLevelEmphasis": "3GLRLM_LRHGLE_AVE",
    "glrlm_LongRunLowGrayLevelEmphasis": "3GLRLM_LRLGLE_AVE",
    "glrlm_LowGrayLevelRunEmphasis": "3GLRLM_LGLRE_AVE",
    "glrlm_RunEntropy": "3GLRLM_RE_AVE",
    "glrlm_RunLengthNonUniformity": "3GLRLM_RLN_AVE",
    "glrlm_RunLengthNonUniformityNormalized": "3GLRLM_RLNN_AVE",
    "glrlm_RunPercentage": "3GLRLM_RP_AVE",
    "glrlm_RunVariance": "3GLRLM_RV_AVE",
    "glrlm_ShortRunEmphasis": "3GLRLM_SRE_AVE",
    "glrlm_ShortRunHighGrayLevelEmphasis": "3GLRLM_SRHGLE_AVE",
    "glrlm_ShortRunLowGrayLevelEmphasis": "3GLRLM_SRLGLE_AVE",
    # glszm
    "glszm_GrayLevelNonUniformity": "3GLSZM_GLN",
    "glszm_GrayLevelNonUniformityNormalized": "3GLSZM_GLNN",
    "glszm_GrayLevelVariance": "3GLSZM_GLV",
    "glszm_HighGrayLevelZoneEmphasis": "3GLSZM_HGLZE",
    "glszm_LargeAreaEmphasis": "3GLSZM_LAE",
    "glszm_LargeAreaHighGrayLevelEmphasis": "3GLSZM_LAHGLE",
    "glszm_LargeAreaLowGrayLevelEmphasis": "3GLSZM_LALGLE",
    "glszm_LowGrayLevelZoneEmphasis": "3GLSZM_LGLZE",
    "glszm_SizeZoneNonUniformity": "3GLSZM_SZN",
    "glszm_SizeZoneNonUniformityNormalized": "3GLSZM_SZNN",
    "glszm_SmallAreaEmphasis": "3GLSZM_SAE",
    "glszm_SmallAreaHighGrayLevelEmphasis": "3GLSZM_SAHGLE",
    "glszm_SmallAreaLowGrayLevelEmphasis": "3GLSZM_SALGLE",
    "glszm_ZoneEntropy": "3GLSZM_ZE",
    "glszm_ZonePercentage": "3GLSZM_ZP",
    "glszm_ZoneVariance": "3GLSZM_ZV",
    # gldm
    "gldm_DependenceEntropy": "3GLDM_DE",
    "gldm_DependenceNonUniformity": "3GLDM_DN",
    "gldm_DependenceNonUniformityNormalized": "3GLDM_DNN",
    "gldm_DependenceVariance": "3GLDM_DV",
    "gldm_GrayLevelNonUniformity": "3GLDM_GLN",
    "gldm_GrayLevelVariance": "3GLDM_GLV",
    "gldm_HighGrayLevelEmphasis": "3GLDM_HGLE",
    "gldm_LargeDependenceEmphasis": "3GLDM_LDE",
    "gldm_LargeDependenceHighGrayLevelEmphasis": "3GLDM_LDHGLE",
    "gldm_LargeDependenceLowGrayLevelEmphasis": "3GLDM_LDLGLE",
    "gldm_LowGrayLevelEmphasis": "3GLDM_LGLE",
    "gldm_SmallDependenceEmphasis": "3GLDM_SDE",
    "gldm_SmallDependenceHighGrayLevelEmphasis": "3GLDM_SDHGLE",
    "gldm_SmallDependenceLowGrayLevelEmphasis": "3GLDM_SDLGLE",
    # ngtdm
    "ngtdm_Busyness": "3NGTDM_BUSYNESS",
    "ngtdm_Coarseness": "3NGTDM_COARSENESS",
    "ngtdm_Complexity": "3NGTDM_COMPLEXITY",
    "ngtdm_Contrast": "3NGTDM_CONTRAST",
    "ngtdm_Strength": "3NGTDM_STRENGTH",
}


def read_pyradiomics(path):
    """Read 'original_<class>_<feature>' rows, returned as '<class>_<feature>'."""
    features = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            if r["feature"].startswith("original_"):
                features[r["feature"].removeprefix("original_")] = float(r["value"])
    return features


def read_nyxus(path, roi):
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


def compare_case(pyr, nyx, rtol, atol):
    rows = []
    for pyr_name, nyx_name in MAPPING.items():
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
        if pyr_name not in MAPPING:
            rows.append([pyr_name, "", pv, "", "", "no Nyxus equivalent"])
    return rows


def write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", help="case name, e.g. PE9f4fd4 (default: all cases)")
    parser.add_argument("--roi", type=int, default=1, help="mask label to compare (default: 1)")
    parser.add_argument("--rtol", type=float, default=0.01, help="relative tolerance for 'close'")
    parser.add_argument("--atol", type=float, default=1e-9, help="absolute tolerance near zero")
    args = parser.parse_args()

    if args.case:
        cases = [args.case]
    else:
        cases = sorted(p.name.removesuffix(".nii.gz") for p in INTENSITY_DIR.glob("*.nii.gz"))

    OUT_DIR.mkdir(exist_ok=True)
    results = {}  # case -> rows
    for case in cases:
        pyr_path = PYRADIOMICS_DIR / f"{case}_label{args.roi}.csv"
        nyx_path = NYXUS_DIR / f"{case}_nyxus.csv"
        if not pyr_path.exists() or not nyx_path.exists():
            print(f"Skipping {case}: missing {pyr_path.name if not pyr_path.exists() else nyx_path.name}")
            continue
        nyx = read_nyxus(nyx_path, args.roi)
        if nyx is None:
            print(f"Skipping {case}: no ROI_label {args.roi} in {nyx_path.name}")
            continue

        rows = compare_case(read_pyradiomics(pyr_path), nyx, args.rtol, args.atol)
        write_csv(OUT_DIR / f"{case}_label{args.roi}.csv",
                  ["pyradiomics_feature", "nyxus_feature", "pyradiomics_value", "nyxus_value",
                   "rel_diff_percent", "status"], rows)
        results[case] = rows
        counts = {s: sum(1 for r in rows if r[5] == s) for s in STATUS_ORDER}
        print(f"{case}: " + ", ".join(f"{s} {n}" for s, n in counts.items() if n))

    if not results:
        return

    # Summary: one row per feature with the worst status and largest difference across cases.
    names = {r[0]: r[1] for rows in results.values() for r in rows}
    summary = []
    for feat, nyx_name in names.items():
        per_case = {case: next((r for r in rows if r[0] == feat), None) for case, rows in results.items()}
        statuses = [r[5] if r else "missing" for r in per_case.values()]
        worst = min(statuses, key=STATUS_ORDER.index)
        diffs = [r[4] for r in per_case.values() if r and r[4] != ""]
        summary.append([feat, nyx_name, worst, max(diffs) if diffs else ""] + statuses)
    summary.sort(key=lambda r: STATUS_ORDER.index(r[2]))
    summary_path = OUT_DIR / f"summary_label{args.roi}.csv"
    write_csv(summary_path,
              ["pyradiomics_feature", "nyxus_feature", "worst_status", "max_rel_diff_percent"]
              + list(results), summary)

    print(f"\nSummary over {len(results)} case(s), label {args.roi}, worst status per feature:")
    for status in STATUS_ORDER:
        n = sum(1 for r in summary if r[2] == status)
        if n:
            print(f"  {status}: {n}")
    for r in summary:
        if r[2] == "different":
            bad = [case for case, s in zip(results, r[4:]) if s == "different"]
            print(f"  different: {r[0]:42s} max rel_diff {r[3]:6.2f}%  in {', '.join(bad)}")
    print(f"Wrote {OUT_DIR.relative_to(ROOT)}/<case>_label{args.roi}.csv and {summary_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
