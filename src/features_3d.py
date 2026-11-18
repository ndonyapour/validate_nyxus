"""
Shared 3D pipeline: extract Nyxus and PyRadiomics features from NIfTI files and compare them.

Used by ct_3d.py and adni_3d.py. For each case, paths(case) returns the intensity and mask
files; every label in `labels` is compared separately.

Outputs (under out_dir, or its bins<N>_d<D> subfolder for non-default settings):
    settings.txt                                bin count and GLCM distance used
    nyxus_out/<case>_nyxus.csv                  one row per mask label
    pyradiomics_out/<case>_label<roi>.csv
    comparisons/<case>_label<roi>.csv and comparisons/summary_label<roi>.csv

Both tools are configured to compute the same quantities, so no conversions are applied:
    coarse_gray_depth=-N   PyRadiomics-style binning into N bins (binCount: N) for all texture
                           families and 3ENTROPY / 3UNIFORMITY; N = 0 bins nothing in Nyxus and
                           uses binWidth: 1 in PyRadiomics
    3glcm/offset=D         GLCM neighbour distance (distances: [D], GLCM only)
    anisotropy_x/y/z       the voxel spacing from the NIfTI header, so shape features are in mm
                           like PyRadiomics instead of voxel units
For N = 0 and D > 1 the tools define some features differently; see make_pyradiomics in common.py.

Nyxus 0.12.1 truncates floating-point NIfTI intensities to integers; PyRadiomics does not.
"""

import logging
import time

import SimpleITK as sitk
from nyxus import Nyxus3D

from common import (ROOT, compare_case, make_pyradiomics, read_nyxus, read_pyradiomics,
                    settings_dir, write_case, write_pyradiomics, write_settings, write_summary)

SETTINGS = ROOT / "settings_3d.yaml"

# PyRadiomics feature -> Nyxus column
# GLCM and GLRLM use the _AVE columns (PyRadiomics averages per-angle values).
MAPPING_3D = {
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


def run(cases, paths, labels, out_dir, rtol, atol, bin_count, glcm_distance):
    out_dir = settings_dir(out_dir, bin_count, glcm_distance)
    nyxus_dir = out_dir / "nyxus_out"
    pyradiomics_dir = out_dir / "pyradiomics_out"
    comparisons_dir = out_dir / "comparisons"
    for d in (nyxus_dir, pyradiomics_dir, comparisons_dir):
        d.mkdir(parents=True, exist_ok=True)
    write_settings(out_dir, SETTINGS, bin_count, glcm_distance)

    logging.getLogger("radiomics").setLevel(logging.ERROR)
    pyradiomics = make_pyradiomics(SETTINGS, bin_count, glcm_distance)

    results = {roi: {} for roi in labels}  # roi -> case -> rows
    for case in cases:
        intensity, mask = paths(case)
        if not intensity.exists() or not mask.exists():
            print(f"Skipping {case}: missing {intensity if not intensity.exists() else mask}")
            continue

        # Spacing is fixed when Nyxus3D is constructed, so build one instance per case.
        reader = sitk.ImageFileReader()
        reader.SetFileName(str(intensity))
        reader.ReadImageInformation()  # header only
        sx, sy, sz = reader.GetSpacing()
        nyx = Nyxus3D(["*3D_ALL*"], coarse_gray_depth=-bin_count,
                      anisotropy_x=sx, anisotropy_y=sy, anisotropy_z=sz)
        nyx.set_metaparam(f"3glcm/offset={glcm_distance}")

        # Nyxus resolves relative paths against /, so pass absolute ones.
        start_time = time.perf_counter()
        df = nyx.featurize_files(intensity_files=[str(intensity.absolute())],
                                 mask_files=[str(mask.absolute())], single_roi=False)
        nyx_time = time.perf_counter() - start_time
        nyx_path = nyxus_dir / f"{case}_nyxus.csv"
        df.to_csv(nyx_path, index=False)
        print(f"{case}: spacing ({sx:.4g}, {sy:.4g}, {sz:.4g}) mm, Nyxus {nyx_time:.1f} s (all labels)")

        for roi in labels:
            nyx_features = read_nyxus(nyx_path, roi)
            if nyx_features is None:
                print(f"  label {roi}: not in the mask, skipped")
                continue
            start_time = time.perf_counter()
            result = pyradiomics(str(intensity), str(mask), roi)
            pyr_time = time.perf_counter() - start_time
            pyr_path = pyradiomics_dir / f"{case}_label{roi}.csv"
            write_pyradiomics(pyr_path, result)

            print(f"  label {roi}: PyRadiomics {pyr_time:.1f} s")
            rows = compare_case(read_pyradiomics(pyr_path), nyx_features, MAPPING_3D, rtol, atol)
            write_case(comparisons_dir, case, roi, rows)
            results[roi][case] = rows

    for roi in labels:
        if results[roi]:
            write_summary(comparisons_dir, roi, results[roi])
