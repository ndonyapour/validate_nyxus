"""
Extract PyRadiomics features for the CT volumes in CT_volumes/.

Usage:
    python src/pyradiomics_features.py                  # all cases, label 1
    python src/pyradiomics_features.py --roi 2          # all cases, label 2
    python src/pyradiomics_features.py --case PE9f4fd4  # one case

Inputs : CT_volumes/intensity/<case>.nii.gz and CT_volumes/segmentations/<case>.nii
Output : pyradiomics_out/<case>_label<roi>.csv (columns: feature, value)

Settings come from settings_hu1.yaml; --roi overrides its label.
"""

import argparse
import csv
import logging
import time
from pathlib import Path

import numpy as np
from radiomics import featureextractor

ROOT = Path(__file__).resolve().parents[1]
INTENSITY_DIR = ROOT / "CT_volumes" / "intensity"
MASK_DIR = ROOT / "CT_volumes" / "segmentations"
OUT_DIR = ROOT / "pyradiomics_out"
SETTINGS = ROOT / "settings_hu1.yaml"


def to_value(v):
    """Feature values are 0-d numpy arrays, and some shape features are complex (x+0j)."""
    if isinstance(v, (np.ndarray, np.generic, int, float, complex)) and np.ndim(v) == 0:
        return float(np.real(v))
    return v


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", help="case name, e.g. PE9f4fd4 (default: all cases)")
    parser.add_argument("--roi", type=int, default=1, help="mask label to process (default: 1)")
    args = parser.parse_args()

    if args.case:
        cases = [args.case]
    else:
        cases = sorted(p.name.removesuffix(".nii.gz") for p in INTENSITY_DIR.glob("*.nii.gz"))

    logging.getLogger("radiomics").setLevel(logging.ERROR)
    extractor = featureextractor.RadiomicsFeatureExtractor(str(SETTINGS))

    OUT_DIR.mkdir(exist_ok=True)
    for case in cases:
        intensity = INTENSITY_DIR / f"{case}.nii.gz"
        mask = MASK_DIR / f"{case}.nii"
        output = OUT_DIR / f"{case}_label{args.roi}.csv"
        if not intensity.exists() or not mask.exists():
            print(f"Skipping {case}: missing {intensity if not intensity.exists() else mask}")
            continue

        start_time = time.perf_counter()
        result = extractor.execute(str(intensity), str(mask), label=args.roi)
        elapsed = time.perf_counter() - start_time

        with open(output, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["feature", "value"])
            for name, value in result.items():
                w.writerow([name, to_value(value)])
        print(f"{case}: label {args.roi}, extraction took {elapsed:.1f} seconds")
        print(f"  Saved to: {output}")


if __name__ == "__main__":
    main()
