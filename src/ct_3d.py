"""
Compare Nyxus and PyRadiomics 3D features on the CT volumes.

Usage:
    python src/ct_3d.py                  # all cases, label 1
    python src/ct_3d.py --labels 1 2     # labels 1 and 2
    python src/ct_3d.py --case PE9f4fd4  # one case
    python src/ct_3d.py --n 5            # 5 random cases

Inputs : <root>/intensity/<case>.nii.gz and <root>/segmentations/<case>.nii (or .nii.gz)
Outputs: results/ct_3d/{nyxus_out,pyradiomics_out,comparisons}/, see features_3d.py
"""

from pathlib import Path

import features_3d
from common import RESULTS_DIR, make_parser, select_cases

CT_ROOT = Path("/polus2/donyapourn2/data/CT_volumes")


def main():
    parser = make_parser(__doc__.strip().splitlines()[0], RESULTS_DIR / "ct_3d")
    parser.add_argument("--labels", type=int, nargs="+", default=[1], help="mask labels (default: 1)")
    parser.add_argument("--root", type=Path, default=CT_ROOT,
                        help=f"folder with intensity/ and segmentations/ (default: {CT_ROOT})")
    args = parser.parse_args()

    intensity_dir = args.root / "intensity"
    mask_dir = args.root / "segmentations"
    available = [p.name.removesuffix(".nii.gz") for p in intensity_dir.glob("*.nii.gz")]
    cases = select_cases(parser, args, available)

    def paths(case):
        mask = mask_dir / f"{case}.nii"
        if not mask.exists():
            mask = mask_dir / f"{case}.nii.gz"
        return intensity_dir / f"{case}.nii.gz", mask

    features_3d.run(cases, paths, args.labels, args.out, args.rtol, args.atol,
                    args.bin_count, args.glcm_distance)


if __name__ == "__main__":
    main()
