"""
Segment the five lung lobes in the INSPECT CTPA scans with TotalSegmentator.

Usage:
    python src/segment_inspect.py                  # all scans
    python src/segment_inspect.py --n 10           # 10 random scans (seed 0)
    python src/segment_inspect.py --case PEc2b4b8  # one scan
    python src/segment_inspect.py --device mps     # Apple GPU; default: CUDA GPU, else CPU

Needs the segmentation extra: uv sync --extra segmentation
The first run downloads the TotalSegmentator model weights to ~/.totalsegmentator.

Inputs : <root>/<case>.nii.gz
Outputs: <out>/<case>.nii.gz, one multilabel mask per scan with the labels in LOBES.
Scans that already have a mask are skipped unless --force is given.

Same settings as TotalSegmentator -ta total_v3 -rs <lobes> --ml.
"""

import argparse
import time
from pathlib import Path

import numpy as np
import SimpleITK as sitk

from common import RESULTS_DIR, select_cases

INSPECT_ROOT = Path("/Users/donyapourn2/Desktop/projects/Inspect_data/INSPECT_Sample_Datatset/CTPA")
SEGMENTATION_DIR = RESULTS_DIR / "inspect_3d" / "segmentations"
TASK = "total_v3"
# TotalSegmentator class -> label in the multilabel mask (class_map["total_v3"])
LOBES = {
    "lung_upper_lobe_left": 10,
    "lung_lower_lobe_left": 11,
    "lung_upper_lobe_right": 12,
    "lung_middle_lobe_right": 13,
    "lung_lower_lobe_right": 14,
}


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--n", type=int, help="number of random scans (default: all)")
    parser.add_argument("--seed", type=int, default=0, help="random seed for --n (default: 0)")
    parser.add_argument("--case", action="append", help="segment only this scan (repeatable)")
    parser.add_argument("--root", type=Path, default=INSPECT_ROOT, help=f"CTPA folder (default: {INSPECT_ROOT})")
    parser.add_argument("--out", type=Path, default=SEGMENTATION_DIR,
                        help=f"mask folder (default: {SEGMENTATION_DIR})")
    parser.add_argument("--device", default="gpu", help="gpu, gpu:N, mps or cpu (default: gpu, falls back to cpu)")
    parser.add_argument("--force", action="store_true", help="segment again even if the mask exists")
    args = parser.parse_args()

    # Imported here so the comparison script can use LOBES without TotalSegmentator installed.
    from totalsegmentator.map_to_binary import class_map
    from totalsegmentator.python_api import totalsegmentator

    expected = {name: class_map[TASK].get(label) for name, label in LOBES.items()}
    assert all(name == found for name, found in expected.items()), f"label ids changed: {expected}"

    images = {p.name.removesuffix(".nii.gz"): p for p in args.root.glob("*.nii.gz")}
    cases = select_cases(parser, args, images)
    args.out.mkdir(parents=True, exist_ok=True)

    failed = []
    for i, case in enumerate(cases, 1):
        output = args.out / f"{case}.nii.gz"
        if output.exists() and not args.force:
            print(f"[{i}/{len(cases)}] {case}: mask exists, skipped")
            continue
        # Write to a temporary name so an interrupted run never leaves a partial mask behind.
        partial = args.out / f"{case}.partial.nii.gz"
        start_time = time.perf_counter()
        try:
            totalsegmentator(images[case], partial, task=TASK, roi_subset=list(LOBES), ml=True,
                             device=args.device, quiet=True)
        except Exception as e:
            print(f"[{i}/{len(cases)}] {case}: FAILED: {e}")
            failed.append(case)
            partial.unlink(missing_ok=True)
            continue
        elapsed = time.perf_counter() - start_time

        labels = set(np.unique(sitk.GetArrayViewFromImage(sitk.ReadImage(str(partial)))).tolist())
        missing = [name for name, label in LOBES.items() if label not in labels]
        if len(missing) == len(LOBES):
            print(f"[{i}/{len(cases)}] {case}: FAILED: the mask is empty")
            failed.append(case)
            partial.unlink()
            continue
        partial.rename(output)
        note = f" (no {', '.join(missing)})" if missing else ""
        print(f"[{i}/{len(cases)}] {case}: {elapsed:.0f} s -> {output}{note}")

    if failed:
        raise SystemExit(f"Failed: {', '.join(failed)}")


if __name__ == "__main__":
    main()
