"""
Shared 2D pipeline: extract Nyxus and PyRadiomics features from in-memory images and compare them.

Used by xray_2d.py and grape_2d.py. For each case, load(case) returns a 2D integer image and a
2D label mask of the same shape; every label in `labels` is compared separately.

Outputs (under out_dir, or its bins<N>_d<D> subfolder for non-default settings):
    settings.txt                                bin count and GLCM distance used
    nyxus_out/<case>_nyxus.csv                  one row per mask label
    pyradiomics_out/<case>_label<roi>.csv
    comparisons/<case>_label<roi>.csv and comparisons/summary_label<roi>.csv

Both tools are configured to compute the same quantities:
    coarse_gray_depth=-N   PyRadiomics-style binning into N bins (binCount: N); N = 0 bins
                           nothing in Nyxus and uses binWidth: 1 in PyRadiomics
    glcm/offset=D          GLCM neighbour distance (distances: [D], GLCM only)
    force2D                PyRadiomics gets a single-slice volume and computes 2D features
For N = 0 and D > 1 the tools define some features differently; see make_pyradiomics in common.py.
"""

import logging
import time

import numpy as np
import SimpleITK as sitk
from nyxus import Nyxus

from common import (ROOT, compare_case, make_pyradiomics, read_nyxus, read_pyradiomics,
                    settings_dir, write_case, write_pyradiomics, write_settings, write_summary)
from features_3d import MAPPING_3D

SETTINGS = ROOT / "settings_2d.yaml"

# PyRadiomics feature -> Nyxus 2D column.
# First order and texture: the 3D names without their leading "3".
MAPPING_2D = {
    "shape2D_PixelSurface": "AREA_PIXELS_COUNT",
    "shape2D_MajorAxisLength": "MAJOR_AXIS_LENGTH",
    "shape2D_MinorAxisLength": "MINOR_AXIS_LENGTH",
    "shape2D_Elongation": "ELONGATION",
    "shape2D_MaximumDiameter": "STAT_FERET_DIAM_MAX",  # Nyxus samples the Feret diameter over angles
    "shape2D_Perimeter": "PERIMETER",                   # PyRadiomics: marching-squares mesh perimeter
    "shape2D_Sphericity": "CIRCULARITY",                # both 2*sqrt(pi*A)/P, with their own perimeter
    **{pyr: nyx.removeprefix("3") for pyr, nyx in MAPPING_3D.items() if not pyr.startswith("shape_")},
}
SHAPE_FEATURES = [nyx for pyr, nyx in MAPPING_2D.items() if pyr.startswith("shape2D_")]
NYXUS_FEATURES = ["*ALL_INTENSITY*", "*ALL_GLCM*", "*ALL_GLRLM*", "*ALL_GLSZM*", "*ALL_GLDM*",
                  "*ALL_NGTDM*"] + SHAPE_FEATURES


def run(cases, load, labels, out_dir, rtol, atol, bin_count, glcm_distance):
    out_dir = settings_dir(out_dir, bin_count, glcm_distance)
    nyxus_dir = out_dir / "nyxus_out"
    pyradiomics_dir = out_dir / "pyradiomics_out"
    comparisons_dir = out_dir / "comparisons"
    for d in (nyxus_dir, pyradiomics_dir, comparisons_dir):
        d.mkdir(parents=True, exist_ok=True)
    write_settings(out_dir, SETTINGS, bin_count, glcm_distance)

    logging.getLogger("radiomics").setLevel(logging.ERROR)
    pyradiomics = make_pyradiomics(SETTINGS, bin_count, glcm_distance)
    nyx = Nyxus(NYXUS_FEATURES, coarse_gray_depth=-bin_count)
    nyx.set_metaparam(f"glcm/offset={glcm_distance}")

    results = {roi: {} for roi in labels}  # roi -> case -> rows
    for case in cases:
        image, mask = load(case)
        if image.shape != mask.shape:
            print(f"Skipping {case}: image {image.shape} and mask {mask.shape} differ")
            continue

        start_time = time.perf_counter()
        df = nyx.featurize(image[None].astype(np.uint32), mask[None].astype(np.uint32))
        nyx_time = time.perf_counter() - start_time
        nyx_path = nyxus_dir / f"{case}_nyxus.csv"
        df.to_csv(nyx_path, index=False)
        print(f"{case}: {image.shape[1]} x {image.shape[0]} px, Nyxus {nyx_time:.1f} s (all labels)")

        # PyRadiomics wants a volume; force2D treats the single slice as a 2D image.
        sitk_image = sitk.GetImageFromArray(image[None].astype(np.float64))
        sitk_mask = sitk.GetImageFromArray(mask[None].astype(np.uint8))
        for roi in labels:
            nyx_features = read_nyxus(nyx_path, roi)
            if nyx_features is None:
                print(f"  label {roi}: not in the mask, skipped")
                continue
            start_time = time.perf_counter()
            result = pyradiomics(sitk_image, sitk_mask, roi)
            pyr_time = time.perf_counter() - start_time
            pyr_path = pyradiomics_dir / f"{case}_label{roi}.csv"
            write_pyradiomics(pyr_path, result)

            print(f"  label {roi}: PyRadiomics {pyr_time:.1f} s")
            rows = compare_case(read_pyradiomics(pyr_path), nyx_features, MAPPING_2D, rtol, atol)
            write_case(comparisons_dir, case, roi, rows)
            results[roi][case] = rows

    for roi in labels:
        if results[roi]:
            write_summary(comparisons_dir, roi, results[roi])
