# validate_nyxus

Compares the 3D radiomics features from [Nyxus](https://github.com/PolusAI/nyxus) with those from [PyRadiomics](https://github.com/AIM-Harvard/pyradiomics) for each case in `CT_volumes/`. It classifies each feature as identical, close or different.

## Files

| File | Purpose |
|---|---|
| `run_all.sh` | Runs the whole pipeline (Nyxus, PyRadiomics, comparison, summary) for every case |
| `compare_nyxus_pyradiomics.py` | Matches each PyRadiomics feature to its Nyxus column, compares the two values, and writes a comparison CSV for one case |
| `summarize_comparisons.py` | Combines the per-case comparison CSVs into one summary table |
| `validate_nyxus.py` | Runs Nyxus 3D on the cases in `CT_volumes/` and saves one feature CSV per case |
| `CT_volumes/` | CT intensity volumes (`intensity/`) and their masks (`segmentations/`) |
| `settings_hu1.yaml` | PyRadiomics parameter file (`binCount: 20`, label 1, original image only, all feature classes) |
| `command.txt` | The three pipeline steps for one case, as separate commands |

## Installation

You need **Python 3.10–3.13**, because Nyxus only publishes wheels for those versions. Python 3.14 will not work.

PyRadiomics is installed from GitHub instead of PyPI. The PyPI release (3.1.0) has no wheels for Python 3.10 or later, and its source build fails on modern Python. To build PyRadiomics, you need `git` and a C compiler. On macOS, install the compiler with `xcode-select --install`.

### Option A: uv (recommended)

```bash
cd validate_nyxus
uv sync
```

`uv sync` creates `.venv` with the Python version in `.python-version` (3.12), downloading that Python if it isn't already installed. It then installs the exact versions recorded in `uv.lock`. To run a command in the environment, use `uv run`, for example `uv run validate_nyxus.py`. Or activate the environment with `source .venv/bin/activate` and use the commands as written below.

### Option B: venv + pip

```bash
cd validate_nyxus
python3.12 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -e .
```

### Option C: conda

```bash
conda create -n validate_nyxus python=3.12
conda activate validate_nyxus
pip install -e .
```

To check the installation, run:

```bash
python -c "import nyxus, radiomics; print('nyxus', nyxus.__version__, '| pyradiomics', radiomics.__version__)"
```

Installing the project adds two commands: `pyradiomics` and `compare-nyxus-pyradiomics`.

## Running everything

```bash
./run_all.sh                      # all cases, label 1
./run_all.sh --roi 2              # all cases, label 2
./run_all.sh PE9f4fd4 PE9f577d    # only these cases
./run_all.sh --force              # recompute PyRadiomics even if its output exists
./run_all.sh --rtol 0.001         # stricter "close" tolerance (default 0.01 = 1%)
```

For each case, the script runs Nyxus, then PyRadiomics, then the comparison. At the end it writes a summary across all cases. It uses `uv run`, so you don't need to activate the environment. PyRadiomics can take several minutes per case, so if a case's PyRadiomics output already exists, the script reuses it. Use `--force` to recompute it. The label comes from `--roi` and is passed to PyRadiomics with `--setting label:<roi>`, so you don't need to edit the YAML file.

The outputs are:

| Path | Content |
|---|---|
| `nyxus_out/<case>_nyxus.csv` | Nyxus features, one row per mask label |
| `pyradiomics_out/<case>_label<roi>_pyr.txt` | PyRadiomics output for one label |
| `comparisons/<case>_label<roi>_comparison.csv` | Feature-by-feature comparison for one case |
| `comparisons/summary_label<roi>.csv` | One row per feature: the worst status, the largest relative difference, and the status in each case |

## Running the steps by hand

`run_all.sh` runs these three steps for each case. Steps 1 and 2 must use the **same** intensity image and mask.

### 1. Extract features with PyRadiomics

```bash
mkdir -p pyradiomics_out
pyradiomics CT_volumes/intensity/PE9f4fd4.nii.gz CT_volumes/segmentations/PE9f4fd4.nii \
    --param settings_hu1.yaml --setting label:1 -o pyradiomics_out/PE9f4fd4_label1_pyr.txt
```

The output is a text file with lines like `Case-1_original_glcm_Contrast: 0.45`. PyRadiomics only processes one label: the one given by `--setting label:N`, or `label: 1` in the YAML file.

### 2. Extract features with Nyxus

The script reads `CT_volumes/intensity/<case>.nii.gz` and `CT_volumes/segmentations/<case>.nii`. It writes `nyxus_out/<case>_nyxus.csv`, with one row per mask label.

```bash
python validate_nyxus.py                   # all cases in CT_volumes/intensity
python validate_nyxus.py --case PE9f4fd4   # one case
```

The script computes all 3D features (`*3D_ALL*`). Its settings match `settings_hu1.yaml`, so both tools compute the same quantities:

| Setting | PyRadiomics (`settings_hu1.yaml`) | Nyxus (`validate_nyxus.py`) |
|---|---|---|
| Grey-level binning | `binCount: 20` over the ROI's range | `coarse_gray_depth=-20` (PyRadiomics-style, 20 bins). This applies to every texture family and to `3ENTROPY`/`3UNIFORMITY`. |
| Voxel spacing | Read from the image, shape in mm | `anisotropy_x/y/z` set to the spacing in the NIfTI header, shape in mm |
| Resampling | none | none |
| Label | `label: 1` (or `--setting label:N`) | every label. The comparison uses the row given by `--roi`. |

If you change `binCount` in the YAML file, change `BIN_COUNT` in `validate_nyxus.py` to match.

### 3. Compare the two outputs

```bash
python compare_nyxus_pyradiomics.py pyradiomics_out/PE9f4fd4_label1_pyr.txt nyxus_out/PE9f4fd4_nyxus.csv \
    --roi 1 --out comparisons/PE9f4fd4_label1_comparison.csv
```

| Argument | Default | Meaning |
|---|---|---|
| `pyradiomics_txt` | required | PyRadiomics text output from step 1 |
| `nyxus_csv` | required | Nyxus CSV from step 2 |
| `--roi` | `1` | Which `ROI_label` row of the Nyxus CSV to compare. Use the same label as PyRadiomics. |
| `--rtol` | `0.01` | Relative tolerance for `close` (0.01 = 1%) |
| `--atol` | `1e-9` | Absolute tolerance, so values near zero aren't flagged because of floating-point noise |
| `--out` | `comparison.csv` | Path of the output CSV |

To combine several comparisons, run:

```bash
python summarize_comparisons.py comparisons/*_label1_comparison.csv --out comparisons/summary_label1.csv
```

## How features are compared

Two values agree when

```
|nyxus − pyradiomics| ≤ atol + tol · max(|nyxus|, |pyradiomics|)
```

This test works the same whichever tool you treat as the reference, and it behaves sensibly when the values are zero or near zero. A plain percent difference doesn't: it divides by the PyRadiomics value and becomes infinite when that value is 0. Each feature gets one of these statuses:

| Status | Meaning |
|---|---|
| `identical` | Agrees with `tol = 1e-6`, so the only difference is floating-point noise |
| `close` | Agrees with `tol = --rtol` but not `1e-6` |
| `different` | Disagrees by more than `--rtol`. **Investigate these.** |
| `missing` | One of the tools did not produce the feature |
| `no Nyxus equivalent` | A PyRadiomics feature with no entry in `MAPPING` |

The comparison CSV also has `abs_diff` and `rel_diff_percent` columns. `rel_diff_percent` is `|nyxus − pyradiomics| / max(|nyxus|, |pyradiomics|) · 100`.

The script applies no unit conversions or exceptions. The two tools are configured identically, so any difference larger than the tolerance is reported as `different`.

The GLCM and GLRLM features are compared against the Nyxus `_AVE` columns, because PyRadiomics averages over angles. To add a feature pair or change one, edit `MAPPING` in `compare_nyxus_pyradiomics.py`.
