#!/usr/bin/env bash
# Full pipeline. Usage: bash run_all.sh   (from ~/landlside; activates the 'landslide' env)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-$HOME/.conda/envs/landslide/bin/python}
mkdir -p logs
$PY -m pytest -q tests
for s in 01_to_parquet 02_eda_filter 03_labels_glc 04_calibrate 05_run_physics 06_compare_baselines 07_maps; do
  echo "=== $s"
  $PY scripts/$s.py 2>&1 | tee logs/$s.log
done
