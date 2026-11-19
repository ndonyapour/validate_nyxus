#!/usr/bin/env bash
# Run the validation on a random sample of each dataset (not all cases):
#   1. 20 bins (the default) on N2D cases per 2D dataset and N3D per ADNI mask (CT: all 5)
#   2. 8, 32 and 64 bins on a smaller sample of NSWEEP cases
# Each run's output is logged to results/logs/.
#
#   ./run_all.sh                          # N2D=100, N3D=50, NSWEEP=20
#   N2D=10 N3D=5 NSWEEP=3 ./run_all.sh    # quick check
#   nohup ./run_all.sh > run_all.log 2>&1 &
#
# All runs use the same --seed. A failed run is reported and the script moves on.

set -u
cd "$(dirname "$0")"

N2D=${N2D:-100}
N3D=${N3D:-50}
NSWEEP=${NSWEEP:-20}
SWEEP_BINS=${SWEEP_BINS:-"8 32 64"}
SEED=${SEED:-0}
LOGS=results/logs
mkdir -p "$LOGS"

failed=()

run() {
    local name=$1; shift
    echo "[$(date '+%F %T')] $name"
    if ! uv run "$@" --seed "$SEED" > "$LOGS/$name.log" 2>&1; then
        echo "  FAILED, see $LOGS/$name.log"
        failed+=("$name")
    fi
}

min() { echo $(( $1 < $2 ? $1 : $2 )); }
# number of INSPECT scans with a lobe mask
n_inspect() { ls results/inspect_3d/segmentations/ 2>/dev/null | grep -v partial | grep -c '\.nii\.gz$'; }

# run_datasets <n 2D> <n 3D> <bins>
run_datasets() {
    local n2d=$1 n3d=$2 bins=$3
    run "xray_2d_bins$bins"             src/xray_2d.py  --n "$n2d" --bin-count "$bins"
    run "grape_2d_bins$bins"            src/grape_2d.py --n "$n2d" --bin-count "$bins"
    run "ct_3d_bins$bins"               src/ct_3d.py    --bin-count "$bins"  # only 5 cases: all of them
    run "adni_3d_hippocampus_bins$bins" src/adni_3d.py  --n "$n3d" --bin-count "$bins" --mask hippocampus
    run "adni_3d_brain_bins$bins"       src/adni_3d.py  --n "$n3d" --bin-count "$bins" --mask brain
    # needs the lobe masks from src/segment_inspect.py; uses every segmented scan if fewer than n3d
    run "inspect_3d_bins$bins"          src/inspect_3d.py --n "$(min "$n3d" "$(n_inspect)")" --bin-count "$bins"
}

run_datasets "$N2D" "$N3D" 20
for bins in $SWEEP_BINS; do
    run_datasets "$NSWEEP" "$NSWEEP" "$bins"
done

echo "[$(date '+%F %T')] done"
if [ ${#failed[@]} -gt 0 ]; then
    echo "Failed: ${failed[*]}"
    exit 1
fi
