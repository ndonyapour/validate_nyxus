"""
Extract Nyxus 3D features for the CT volumes in CT_volumes/.

Usage:
    python src/nyxus_features.py                  # all cases
    python src/nyxus_features.py --case PE9f4fd4  # one case

Inputs : CT_volumes/intensity/<case>.nii.gz and CT_volumes/segmentations/<case>.nii
Output : nyxus_out/<case>_nyxus.csv (one row per mask label)

Settings are chosen to match settings_hu1.yaml (PyRadiomics):
    coarse_gray_depth=-20  PyRadiomics-style binning into 20 bins (binCount: 20) for all
                           texture families and for 3ENTROPY / 3UNIFORMITY
    anisotropy_x/y/z       the voxel spacing from the NIfTI header, so shape features are
                           in mm like PyRadiomics instead of voxel units
"""

import argparse
import time
from pathlib import Path

import SimpleITK as sitk
from nyxus import Nyxus3D

BIN_COUNT = 20  # must match binCount in settings_hu1.yaml

ROOT = Path(__file__).resolve().parents[1]
INTENSITY_DIR = ROOT / "CT_volumes" / "intensity"
MASK_DIR = ROOT / "CT_volumes" / "segmentations"
OUT_DIR = ROOT / "nyxus_out"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", help="case name, e.g. PE9f4fd4 (default: all cases)")
    args = parser.parse_args()

    if args.case:
        cases = [args.case]
    else:
        cases = sorted(p.name.removesuffix(".nii.gz") for p in INTENSITY_DIR.glob("*.nii.gz"))

    OUT_DIR.mkdir(exist_ok=True)
    for case in cases:
        intensity = INTENSITY_DIR / f"{case}.nii.gz"
        mask = MASK_DIR / f"{case}.nii"
        if not intensity.exists() or not mask.exists():
            print(f"Skipping {case}: missing {intensity if not intensity.exists() else mask}")
            continue

        # Spacing is fixed when Nyxus3D is constructed, so build one instance per case.
        reader = sitk.ImageFileReader()
        reader.SetFileName(str(intensity))
        reader.ReadImageInformation()  # header only
        sx, sy, sz = reader.GetSpacing()
        nyx = Nyxus3D(
            ["*3D_ALL*"],
            coarse_gray_depth=-BIN_COUNT,
            anisotropy_x=sx,
            anisotropy_y=sy,
            anisotropy_z=sz,
        )

        start_time = time.perf_counter()
        df = nyx.featurize_files(
            intensity_files=[str(intensity)],
            mask_files=[str(mask)],
            single_roi=False,
        )
        elapsed = time.perf_counter() - start_time

        output = OUT_DIR / f"{case}_nyxus.csv"
        df.to_csv(output, index=False)
        print(f"{case}: spacing ({sx:.4g}, {sy:.4g}, {sz:.4g}) mm, extraction took {elapsed:.3f} seconds")
        print(f"  Saved to: {output}")


if __name__ == "__main__":
    main()
