"""
Compare Nyxus and PyRadiomics 3D features on INSPECT CTPA scans, for each lung lobe.

Usage:
    python src/segment_inspect.py           # first: segment the lobes (see that script)
    python src/inspect_3d.py                # all segmented scans, all five lobes
    python src/inspect_3d.py --n 10         # 10 random segmented scans
    python src/inspect_3d.py --lobes lung_upper_lobe_left lung_lower_lobe_left

Inputs : <root>/<case>.nii.gz and <seg-dir>/<case>.nii.gz (from segment_inspect.py)
Outputs: results/inspect_3d/{nyxus_out,pyradiomics_out,comparisons}/, see features_3d.py

Each lobe is compared as its own label (TotalSegmentator total_v3 ids):
    10 lung_upper_lobe_left     11 lung_lower_lobe_left
    12 lung_upper_lobe_right    13 lung_middle_lobe_right    14 lung_lower_lobe_right
"""

from pathlib import Path

import features_3d
from common import RESULTS_DIR, make_parser, select_cases
from segment_inspect import INSPECT_ROOT, LOBES, SEGMENTATION_DIR


def main():
    parser = make_parser(__doc__.strip().splitlines()[0], RESULTS_DIR / "inspect_3d")
    parser.add_argument("--lobes", nargs="+", choices=LOBES, default=list(LOBES),
                        help="lobes to compare (default: all five)")
    parser.add_argument("--root", type=Path, default=INSPECT_ROOT, help=f"CTPA folder (default: {INSPECT_ROOT})")
    parser.add_argument("--seg-dir", type=Path, default=SEGMENTATION_DIR,
                        help=f"masks from segment_inspect.py (default: {SEGMENTATION_DIR})")
    args = parser.parse_args()

    images = {p.name.removesuffix(".nii.gz") for p in args.root.glob("*.nii.gz")}
    masks = {p.name.removesuffix(".nii.gz") for p in args.seg_dir.glob("*.nii.gz")
             if not p.name.endswith(".partial.nii.gz")}
    if images - masks:
        print(f"{len(images - masks)} of {len(images)} scans have no mask yet; run segment_inspect.py first")
    cases = select_cases(parser, args, images & masks)

    def paths(case):
        return args.root / f"{case}.nii.gz", args.seg_dir / f"{case}.nii.gz"

    labels = [LOBES[lobe] for lobe in args.lobes]
    features_3d.run(cases, paths, labels, args.out, args.rtol, args.atol,
                    args.bin_count, args.glcm_distance)


if __name__ == "__main__":
    main()
