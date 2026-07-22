#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."
EXECUTABLE_NAME="simulation_fitted_mu_$$"
cleanup() {
    if [ -f "./$EXECUTABLE_NAME" ]; then
        rm -f "./$EXECUTABLE_NAME"
    fi
}
trap cleanup EXIT

gcc -fopenmp simulation.c func.c helpers.c -o "$EXECUTABLE_NAME" -lm

for CONFIG_FILE in ./experiments/configs/ltdb_fitted_mu/*.conf; do
    NUM_THREADS="$(awk -F= '/^num_threads=/ {print $2; exit}' "$CONFIG_FILE")"
    export OMP_NUM_THREADS="${NUM_THREADS:-1}"
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    echo "Running $CONFIG_FILE"
    "./$EXECUTABLE_NAME" "$CONFIG_FILE"
done
