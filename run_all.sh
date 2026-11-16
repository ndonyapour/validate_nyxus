#!/usr/bin/env bash
# Run the full Nyxus vs PyRadiomics comparison for every case in CT_volumes/.
#
# Usage:
#   ./run_all.sh                  # all cases, label 1
#   ./run_all.sh --roi 2          # all cases, label 2
#   ./run_all.sh PE9f4fd4 PE9f577d  # selected cases
#   ./run_all.sh --force          # recompute PyRadiomics even if its output exists
#
# PyRadiomics is slow, so its output is reused when it already exists.
# Outputs: nyxus_out/, pyradiomics_out/, comparisons/ (one CSV per case + summary.csv)

set -euo pipefail
cd "$(dirname "$0")"

ROI=1
FORCE=0
RTOL=0.01
CASES=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --roi) ROI="$2"; shift 2 ;;
    --rtol) RTOL="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    -h|--help) sed -n '2,11p' "$0"; exit 0 ;;
    *) CASES+=("$1"); shift ;;
  esac
done

if [[ ${#CASES[@]} -eq 0 ]]; then
  for img in CT_volumes/intensity/*.nii.gz; do
    CASES+=("$(basename "$img" .nii.gz)")
  done
fi

mkdir -p nyxus_out pyradiomics_out comparisons
OUTPUTS=()

for case in "${CASES[@]}"; do
  echo "=== $case (label $ROI) ==="
  img="CT_volumes/intensity/$case.nii.gz"
  mask="CT_volumes/segmentations/$case.nii"
  pyr="pyradiomics_out/${case}_label${ROI}_pyr.txt"
  nyx="nyxus_out/${case}_nyxus.csv"
  out="comparisons/${case}_label${ROI}_comparison.csv"

  uv run validate_nyxus.py --case "$case"

  if [[ $FORCE -eq 1 || ! -s "$pyr" ]]; then
    echo "Running PyRadiomics (this can take several minutes)..."
    uv run pyradiomics "$img" "$mask" --param settings_hu1.yaml --setting "label:$ROI" -o "$pyr"
  else
    echo "Reusing $pyr (pass --force to recompute)"
  fi

  uv run compare_nyxus_pyradiomics.py "$pyr" "$nyx" --roi "$ROI" --rtol "$RTOL" --out "$out"
  OUTPUTS+=("$out")
  echo
done

uv run summarize_comparisons.py "${OUTPUTS[@]}" --out "comparisons/summary_label${ROI}.csv"
