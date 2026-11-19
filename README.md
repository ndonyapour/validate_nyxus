# validate_nyxus

Compares the radiomics features from [Nyxus](https://github.com/PolusAI/nyxus) with those from [PyRadiomics](https://github.com/AIM-Harvard/pyradiomics), in 3D and 2D, and reports where the two tools disagree.

## Setup

Requirements:
- [uv](https://docs.astral.sh/uv/)
- `git`
- a C compiler, which PyRadiomics needs to build. On macOS, install it with `xcode-select --install`.

```bash
uv sync
```

This creates `.venv` with Python 3.12 and installs Nyxus and PyRadiomics.

The INSPECT lung lobes are segmented with TotalSegmentator, which needs PyTorch. Install it only if you need to segment:

```bash
uv sync --extra segmentation
```

## Run

There is one script per dataset. Each one extracts features with both tools and compares them:

| Script | Dataset | Labels |
|---|---|---|
| `src/ct_3d.py` | CT volumes | `--labels` (default 1) |
| `src/adni_3d.py` | ADNI brain MRI, Nyul-normalized | 1: whole brain or hippocampus (`--mask`) |
| `src/inspect_3d.py` | INSPECT CTPA with TotalSegmentator lung lobes | 10–14: the five lobes (`--lobes`) |
| `src/xray_2d.py` | ChestX-ray14 with CheXmask lung masks | 1: lungs |
| `src/grape_2d.py` | GRAPE fundus photos | 1: optic disc rim, 2: optic cup |

All scripts take the same options:

| Option | Meaning |
|---|---|
| `--n N` | N random cases (default: all) |
| `--seed S` | random seed for `--n` (default: 0) |
| `--case NAME` | only this case; repeat for more |
| `--rtol`, `--atol` | tolerances, see [Results](#results) |
| `--bin-count N` | grey-level bins for texture features; 0 = no binning, integer images only (default: 20) |
| `--glcm-distance D` | GLCM neighbour distance (default: 1) |
| `--out DIR` | output folder (default: `results/<dataset>/`) |

```bash
uv run src/ct_3d.py                          # all CT cases, label 1
uv run src/ct_3d.py --labels 1 2             # labels 1 and 2
uv run src/adni_3d.py --n 10                 # 10 random ADNI cases, whole brain
uv run src/adni_3d.py --n 10 --mask hippocampus
uv run src/inspect_3d.py --n 10              # 10 random INSPECT scans, all five lobes
uv run src/xray_2d.py --n 10
uv run src/grape_2d.py --n 10                # green channel; --channel red/blue/gray
uv run src/xray_2d.py --n 10 --bin-count 64  # writes to results/xray_2d/bins64_d1/
```

Run any script with `--help` to see all its options.

Each run writes, for every case and label:

```
results/<dataset>/
  settings.txt                                 bin count and GLCM distance used
  nyxus_out/<case>_nyxus.csv                   one row per label
  pyradiomics_out/<case>_label<roi>.csv
  comparisons/<case>_label<roi>.csv            one row per feature
  comparisons/summary_label<roi>.csv           one row per feature: worst status, number of
                                               cases that agree, largest difference, and one
                                               status column per case
```

ADNI writes to `results/adni_3d/<mask>/`. A non-default `--bin-count` or `--glcm-distance` writes to a `bins<N>_d<D>/` subfolder, so runs with different settings don't overwrite each other.

## Code

| File | Contents |
|---|---|
| `src/common.py` | command-line options, case selection, comparison and summary |
| `src/features_3d.py` | 3D extraction from NIfTI files, and the PyRadiomics to Nyxus feature mapping |
| `src/features_2d.py` | 2D extraction from in-memory images, and the 2D feature mapping |
| `src/<dataset>_3d.py`, `src/<dataset>_2d.py` | where each dataset's images and masks are, and how to read them |
| `src/segment_inspect.py` | TotalSegmentator lung-lobe segmentation for INSPECT |

## Data

The datasets are read from `/polus2/donyapourn2/data/`; see the top of each script, or pass the dataset's root option (e.g. `--root`) to read them from elsewhere. The CT volumes use this layout:

```
CT_volumes/
  intensity/<case>.nii.gz
  segmentations/<case>.nii       (or .nii.gz)
```

- **INSPECT**: the CTPA scans are read from `<root>/<case>.nii.gz`. Segment them first with `uv run --extra segmentation src/segment_inspect.py`, which runs TotalSegmentator (`total_v3`, lobes only, one multilabel mask per scan) and writes the masks to `results/inspect_3d/segmentations/`. It skips scans that already have a mask, so you can stop and restart it. It is slow on a CPU; use `--device gpu` (the default) on a CUDA machine or `--device mps` on a Mac. The labels are 10 left upper, 11 left lower, 12 right upper, 13 right middle and 14 right lower lobe.
- **ADNI**: the BrainIAC masks cover the whole brain. The hippocampus masks merge the left and right SynthSeg hippocampus labels into label 1.
- **ChestX-ray14**: the 8-bit PNGs are paired with the `label_mask` array in each CheXmask `.npz`.
- **GRAPE**: the LabelMe polygons in `json/` are drawn on `ROI_images/`. The cup lies inside the disc, so label 1 is the rim (disc minus cup) and label 2 is the cup.

## Settings

Both tools are configured to compute the same quantities:

| | PyRadiomics | Nyxus |
|---|---|---|
| Settings | `settings_3d.yaml`, `settings_2d.yaml` | `src/features_3d.py`, `src/features_2d.py` |
| Grey-level binning (`--bin-count N`) | `binCount: N` | `coarse_gray_depth=-N` |
| No binning (`--bin-count 0`) | `binWidth: 1` | `coarse_gray_depth=0` |
| GLCM distance (`--glcm-distance D`) | `distances: [D]`, GLCM only | `glcm/offset=D`, `3glcm/offset=D` |
| Voxel spacing (3D) | from the image | `anisotropy_x/y/z` from the image |
| Resampling | none | none |
| 2D | `force2D` on a single-slice volume | `Nyxus` (2D) on the array |

The command-line options override the yaml files, which keep the defaults (20 bins, distance 1).

## Results

Each feature gets one of these statuses:

| Status | Meaning |
|---|---|
| `identical` | Agrees to a relative tolerance of 1e-6 |
| `close` | Agrees within 1%. Use `--rtol` to change this tolerance. |
| `different` | Differs by more than 1% |
| `missing` | One of the tools did not produce the feature |
| `no Nyxus equivalent` | No matching Nyxus feature |

The percentage is `|nyxus - pyradiomics| / max(|nyxus|, |pyradiomics|) * 100`.

## Known differences

- **Float images (ADNI)**: Nyxus 0.12.1 truncates floating-point NIfTI intensities to integers ([raw_nifti.h#L373](https://github.com/PolusAI/nyxus/blob/v0.12.1/src/nyx/raw_nifti.h#L373)), while PyRadiomics keeps them. This causes 1-4% differences in some texture features on small ROIs such as the hippocampus.
- **GLCM with empty grey-level bins**: Nyxus builds the GLCM sum and difference distributions from positions in the list of grey levels present, not from the grey levels themselves ([3d_glcm.cpp#L376-L403](https://github.com/PolusAI/nyxus/blob/v0.12.1/src/nyx/features/3d_glcm.cpp#L376-L403), same code in 2D in [glcm.cpp#L485-L512](https://github.com/PolusAI/nyxus/blob/v0.12.1/src/nyx/features/glcm.cpp#L485-L512)). When a ROI leaves a bin empty, DifferenceAverage, DifferenceVariance, InverseVariance and SumAverage are wrong (up to 98% off on CT); Contrast and JointAverage are not affected.
- **2D GLSZM**: the Nyxus 2D zone search only steps right and down ([glszm.cpp#L124-L165](https://github.com/PolusAI/nyxus/blob/v0.12.1/src/nyx/features/glszm.cpp#L124-L165)), so zones that bend back up or left are split. All 2D GLSZM features differ.
- **No binning (`--bin-count 0`)**: PyRadiomics renumbers the grey levels from 1 at the ROI minimum, while Nyxus keeps the image's values. Features of the level values themselves (JointAverage, Autocorrelation, the `*GrayLevel*Emphasis` family) therefore differ, while features of level differences or counts (Contrast, the run, zone and dependence counts) agree. Nyxus' `ENTROPY` and `UNIFORMITY` also switch to 24 bins.
- **GLCM distance above 1**: PyRadiomics uses every offset at Chebyshev distance D (8 directions in 2D at D = 2, 49 in 3D), while Nyxus scales its 4 (2D) or 13 (3D) unit directions by D. The two definitions agree only at D = 1; at D = 2 GLCM features differ by 1-4% on an X-ray and up to 50% on CT.
- **2D perimeter**: `shape2D_Perimeter` is mapped to `PERIMETER` and `shape2D_Sphericity` to `CIRCULARITY`. Each pair uses the same formula, but the two tools measure the perimeter differently, so these pairs disagree by design.
