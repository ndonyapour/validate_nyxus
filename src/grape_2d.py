"""
Compare Nyxus and PyRadiomics 2D features on GRAPE fundus images, for the optic disc rim and cup.

Usage:
    python src/grape_2d.py --n 10                 # 10 random images (seed 0), green channel
    python src/grape_2d.py --n 10 --seed 3        # a different 10
    python src/grape_2d.py --n 10 --channel gray  # luminance instead of the green channel
    python src/grape_2d.py                        # all 631 images

Inputs : <root>/ROI_images/<case>.jpg   (RGB crop around the optic disc)
         <root>/json/<case>.json        (LabelMe polygons "OD" and "OC" in ROI_images pixels)
Outputs: results/grape_2d/{nyxus_out,pyradiomics_out,comparisons}/, see features_2d.py

The cup lies inside the disc and a pixel can carry one label, so the mask is:
    label 1  neuroretinal rim  (OD polygon minus OC polygon)
    label 2  optic cup         (OC polygon)
"""

import json
from pathlib import Path

import numpy as np
import SimpleITK as sitk

import features_2d
from common import RESULTS_DIR, make_parser, select_cases

GRAPE_ROOT = Path("/polus2/donyapourn2/data/GRAPE")
LABELS = {"OD": 1, "OC": 2}  # drawn in this order, so OC overwrites the centre of OD


def polygon_mask(points, shape):
    """Pixels whose centre lies inside the polygon (even-odd rule)."""
    h, w = shape
    y, x = np.mgrid[0:h, 0:w]
    inside = np.zeros(shape, bool)
    px, py = np.asarray(points, float).T
    for x0, y0, x1, y1 in zip(px, py, np.roll(px, -1), np.roll(py, -1)):
        if y0 == y1:
            continue
        crosses = (y0 > y) != (y1 > y)
        x_at_y = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
        inside ^= crosses & (x < x_at_y)
    return inside


def read_image(path, channel):
    rgb = sitk.GetArrayFromImage(sitk.ReadImage(str(path)))  # (h, w, 3), uint8
    if channel == "gray":
        return np.round(rgb @ [0.299, 0.587, 0.114]).astype(np.uint8)
    return rgb[..., "rgb".index(channel[0])]


def read_mask(path, shape):
    mask = np.zeros(shape, np.uint8)
    shapes = {s["label"]: s["points"] for s in json.load(open(path))["shapes"]}
    for name, label in LABELS.items():
        mask[polygon_mask(shapes[name], shape)] = label
    return mask


def main():
    parser = make_parser(__doc__.strip().splitlines()[0], RESULTS_DIR / "grape_2d")
    parser.add_argument("--channel", choices=["green", "red", "blue", "gray"], default="green",
                        help="image channel to analyse (default: green)")
    parser.add_argument("--root", type=Path, default=GRAPE_ROOT)
    args = parser.parse_args()

    images = {p.stem: p for p in (args.root / "ROI_images").glob("*.jpg")}
    annotations = {p.stem: p for p in (args.root / "json").glob("*.json")}
    cases = select_cases(parser, args, images.keys() & annotations.keys())

    def load(case):
        image = read_image(images[case], args.channel)
        return image, read_mask(annotations[case], image.shape)

    features_2d.run(cases, load, sorted(LABELS.values()), args.out, args.rtol, args.atol)


if __name__ == "__main__":
    main()
