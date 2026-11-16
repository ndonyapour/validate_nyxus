#!/usr/bin/env bash
# Run the full Nyxus vs PyRadiomics comparison for every case in CT_volumes/.
#
# Usage:
#   ./run_all.sh            # label 1
#   ./run_all.sh --roi 2    # label 2

set -euo pipefail
cd "$(dirname "$0")"

ROI=1
if [[ "${1:-}" == "--roi" ]]; then
  ROI="$2"
fi

uv run src/nyxus_features.py
uv run src/pyradiomics_features.py --roi "$ROI"
uv run src/compare.py --roi "$ROI"
