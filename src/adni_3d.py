"""
Compare Nyxus and PyRadiomics 3D features on ADNI brain MRIs.

Usage:
    python src/adni_3d.py --n 10                     # 10 random cases, whole-brain mask
    python src/adni_3d.py --n 10 --mask hippocampus  # hippocampus mask instead
    python src/adni_3d.py --n 10 --seed 3            # a different 10

Inputs : <image-root>/<cohort>/images/<case>_0000.nii.gz  (Nyul-normalized MRI, float32)
         and a binary mask (label 1) per case:
           --mask brain        BrainIAC_preprocessed/<cohort>/masks/<case>_0000_mask.nii.gz
           --mask hippocampus  SynthSeg/<cohort>/hippocampus_masks/<case>_0000_hippocampus_mask.nii.gz
Outputs: results/adni_3d/<mask>/{nyxus_out,pyradiomics_out,comparisons}/, see features_3d.py

The images are floating point, which Nyxus 0.12.1 truncates to integers. That causes 1-4%
differences in some texture features on small ROIs such as the hippocampus.
"""

from pathlib import Path

import features_3d
from common import RESULTS_DIR, make_parser, select_cases

DATA = Path("/polus2/donyapourn2/data/ADNI_dataset_processed")
IMAGE_ROOT = DATA / "Nyul_normalized_clipped"
# mask -> (mask root, mask path relative to it; {cohort} and {stem} are filled in per image)
MASKS = {
    "brain": (DATA / "BrainIAC_preprocessed", "{cohort}/masks/{stem}_mask.nii.gz"),
    "hippocampus": (DATA / "SynthSeg", "{cohort}/hippocampus_masks/{stem}_hippocampus_mask.nii.gz"),
}


def main():
    parser = make_parser(__doc__.strip().splitlines()[0], RESULTS_DIR / "adni_3d")
    parser.add_argument("--mask", choices=MASKS, default="brain", help="ROI mask (default: brain)")
    parser.add_argument("--image-root", type=Path, default=IMAGE_ROOT)
    parser.add_argument("--mask-root", type=Path, help="override the root of --mask")
    args = parser.parse_args()
    mask_root, mask_pattern = MASKS[args.mask]
    mask_root = args.mask_root or mask_root

    pairs = {}  # case -> (image, mask)
    for image in args.image_root.glob("*/images/*_0000.nii.gz"):
        cohort = image.parts[-3]
        stem = image.name.removesuffix(".nii.gz")
        mask = mask_root / mask_pattern.format(cohort=cohort, stem=stem)
        if mask.exists():
            pairs[stem.removesuffix("_0000")] = (image, mask)
    cases = select_cases(parser, args, pairs)

    features_3d.run(cases, pairs.get, [1], args.out / args.mask, args.rtol, args.atol,
                    args.bin_count, args.glcm_distance)


if __name__ == "__main__":
    main()
