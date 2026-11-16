"""
Compare Nyxus and PyRadiomics 2D features on ChestX-ray14 images with CheXmask lung masks.

Usage:
    python src/xray_2d.py --n 10              # 10 random images (seed 0)
    python src/xray_2d.py --n 10 --seed 3     # a different 10
    python src/xray_2d.py --case 00026101_000 # one image

Inputs : <image-root>/images_*/<case>.png               (8-bit grayscale, 1024 x 1024)
         <mask-root>/<case>_lung_mask.npz                (array "label_mask", lungs = 1)
Outputs: results/xray_2d/{nyxus_out,pyradiomics_out,comparisons}/, see features_2d.py

The images are integers, so Nyxus's truncation of float intensities does not apply.
"""

from pathlib import Path

import numpy as np
import SimpleITK as sitk

import features_2d
from common import RESULTS_DIR, make_parser, select_cases

DATA = Path("/polus2/donyapourn2/data/Chest_Xray_Images")
IMAGE_ROOT = DATA / "ChestX-ray14"
MASK_ROOT = DATA / "Chest_Xray_Segmentations" / "ChestX-ray14_Chexmask_Segmentations"


def read_image(path):
    """A few ChestX-ray14 PNGs are RGB(A); their channels are equal, so keep the first."""
    a = sitk.GetArrayFromImage(sitk.ReadImage(str(path)))
    return a[..., 0] if a.ndim == 3 else a


def main():
    parser = make_parser(__doc__.strip().splitlines()[0], RESULTS_DIR / "xray_2d")
    parser.add_argument("--image-root", type=Path, default=IMAGE_ROOT)
    parser.add_argument("--mask-root", type=Path, default=MASK_ROOT)
    args = parser.parse_args()

    images = {p.stem: p for p in args.image_root.glob("images_*/*.png")}
    masks = {p.name.removesuffix("_lung_mask.npz"): p for p in args.mask_root.glob("*_lung_mask.npz")}
    cases = select_cases(parser, args, images.keys() & masks.keys())

    def load(case):
        return read_image(images[case]), np.load(masks[case])["label_mask"]

    features_2d.run(cases, load, [1], args.out, args.rtol, args.atol)


if __name__ == "__main__":
    main()
