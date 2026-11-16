# validate_nyxus

Compares the 3D radiomics features from [Nyxus](https://github.com/PolusAI/nyxus) with those from [PyRadiomics](https://github.com/AIM-Harvard/pyradiomics) on CT volumes, and reports where the two tools disagree.

## Setup

Requirements:
- [uv](https://docs.astral.sh/uv/)
- `git`
- a C compiler, which PyRadiomics needs to build. On macOS, install it with `xcode-select --install`.

```bash
uv sync
```

This creates `.venv` with Python 3.12 and installs Nyxus and PyRadiomics.

## Data

The CT volumes are not included in the repository. Put them in this layout:

```
CT_volumes/
  intensity/<case>.nii.gz
  segmentations/<case>.nii
```

## Run

```bash
./run_all.sh            # label 1
./run_all.sh --roi 2    # label 2
```

This runs the three scripts in `src/` for every case:

| Script | Output |
|---|---|
| `src/nyxus_features.py` | `nyxus_out/<case>_nyxus.csv` |
| `src/pyradiomics_features.py` | `pyradiomics_out/<case>_label<roi>.csv` |
| `src/compare.py` | `comparisons/<case>_label<roi>.csv` and `comparisons/summary_label<roi>.csv` |

You can also run each script on its own, e.g. `uv run src/compare.py --case PE9f4fd4`. Run any script with `--help` to see its options.

## Settings

Both tools are configured to compute the same quantities:

| | PyRadiomics (`settings_hu1.yaml`) | Nyxus (`src/nyxus_features.py`) |
|---|---|---|
| Grey-level binning | `binCount: 20` | `coarse_gray_depth=-20` |
| Voxel spacing | from the image | `anisotropy_x/y/z` from the image |
| Resampling | none | none |

If you change `binCount`, change `BIN_COUNT` in `src/nyxus_features.py` to match.

## Results

Each feature gets one of these statuses:

| Status | Meaning |
|---|---|
| `identical` | Agrees to a relative tolerance of 1e-6 |
| `close` | Agrees within 1%. Use `--rtol` to change this tolerance. |
| `different` | Differs by more than 1% |
| `missing` | One of the tools did not produce the feature |
| `no Nyxus equivalent` | No matching Nyxus feature |
