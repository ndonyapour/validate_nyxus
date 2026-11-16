"""
Compare Nyxus 3D features against PyRadiomics features for one case.

Usage:
    python compare_nyxus_pyradiomics.py pyr.txt nyxus.csv [--roi 1] [--rtol 0.01] [--atol 1e-9] [--out comparison.csv]

pyr.txt   : text output of the `pyradiomics` CLI (lines like Case-1_original_glcm_Contrast: 0.45)
nyxus.csv : output of Nyxus3D.featurize_files(...).to_csv(...)
--roi     : which ROI_label row of the Nyxus CSV to use (PyRadiomics only processes one label)
--rtol    : relative tolerance for "close" (default 0.01 = 1%)
--atol    : absolute tolerance, so values near zero are not flagged (default 1e-9)

Two values agree when |nyx - pyr| <= atol + tol * max(|nyx|, |pyr|).
Status per feature:
    identical          agrees with tol = 1e-6 (same implementation, float noise only)
    close              agrees with tol = --rtol
    different          disagrees with tol = --rtol
    missing            one of the tools did not produce the feature
    no Nyxus equivalent  PyRadiomics feature that is not in MAPPING

No conversions are applied: the Nyxus CSV must be computed with the same settings as
PyRadiomics (see validate_nyxus.py), including the voxel spacing, so shape features are in mm.
"""

import argparse
import csv
import re

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

IDENTICAL_RTOL = 1e-6


def parse_number(text):
    """Parse a PyRadiomics value; shape features can be printed as complex numbers like (1.5+0j)."""
    text = text.strip()
    try:
        return float(text)
    except ValueError:
        return complex(text).real


def read_pyradiomics(path):
    features = {}
    spacing = None
    with open(path) as f:
        for line in f:
            m = re.match(r"Case-\d+_original_(\w+?)_(\w+):\s*(\S+)\s*$", line)
            if m:
                try:
                    features[f"{m.group(1)}_{m.group(2)}"] = parse_number(m.group(3))
                except ValueError:
                    pass
                continue
            m = re.match(r"Case-\d+_diagnostics_Image-original_Spacing:\s*\(([^)]*)\)", line)
            if m:
                spacing = [float(s) for s in m.group(1).split(",")]
    return features, spacing


def read_nyxus(path, roi):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    if roi is not None:
        rows = [r for r in rows if str(r["ROI_label"]) == str(roi)]
    if not rows:
        raise SystemExit(f"No row with ROI_label={roi} in {path}")
    row = rows[0]
    out = {}
    for k, v in row.items():
        try:
            out[k] = float(v)
        except (ValueError, TypeError):
            pass
    return out


def agrees(a, b, rtol, atol):
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b))


def rel_diff(a, b):
    scale = max(abs(a), abs(b))
    return abs(a - b) / scale if scale else 0.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pyradiomics_txt")
    parser.add_argument("nyxus_csv")
    parser.add_argument("--roi", default="1")
    parser.add_argument("--rtol", type=float, default=0.01)
    parser.add_argument("--atol", type=float, default=1e-9)
    parser.add_argument("--out", default="comparison.csv")
    args = parser.parse_args()

    pyr, spacing = read_pyradiomics(args.pyradiomics_txt)
    nyx = read_nyxus(args.nyxus_csv, args.roi)

    rows = []
    for pyr_name, nyx_name in MAPPING.items():
        pv = pyr.get(pyr_name)
        nv = nyx.get(nyx_name)
        if pv is None or nv is None:
            rows.append([pyr_name, nyx_name, pv, nv, "", "", "missing"])
            continue
        if agrees(nv, pv, IDENTICAL_RTOL, args.atol):
            status = "identical"
        elif agrees(nv, pv, args.rtol, args.atol):
            status = "close"
        else:
            status = "different"
        rows.append([pyr_name, nyx_name, pv, nv, abs(nv - pv), rel_diff(nv, pv) * 100, status])

    for pyr_name in pyr:
        if pyr_name not in MAPPING:
            rows.append([pyr_name, "", pyr[pyr_name], "", "", "", "no Nyxus equivalent"])

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pyradiomics_feature", "nyxus_feature", "pyradiomics_value", "nyxus_value",
                    "abs_diff", "rel_diff_percent", "status"])
        w.writerows(rows)

    counts = {}
    for r in rows:
        counts[r[6]] = counts.get(r[6], 0) + 1
    print(f"ROI_label: {args.roi}  rtol: {args.rtol}  atol: {args.atol}  spacing: {spacing}")
    for status, n in sorted(counts.items()):
        print(f"  {status}: {n}")
    print(f"Wrote {args.out}")

    hits = [r for r in rows if r[6] == "different"]
    if hits:
        print("\nDifferent:")
        for r in hits:
            print(f"  {r[0]:45s} pyr={r[2]:<14.6g} nyx={r[3]:<14.6g} rel_diff={r[5]:.2f}%")


if __name__ == "__main__":
    main()
