#!/bin/bash
set -euo pipefail

# === OAR directives ===
#OAR -l /nodes=1/core=50,walltime=3:00:00
#OAR -O OAR_%jobid%.out
#OAR -E OAR_%jobid%.err
#OAR -n fixed-surface

cd "$(dirname "$0")/.."

SOURCE_FILE="simulation.c"
LIB_FILE="func.c helpers.c"
EXECUTABLE_NAME="simulation_executable_$$"  # Executable with unique name
CONFIG_FILE="./experiments/configs/detection_time_fixed_surface.conf"

cleanup() {
    if [ -f "./$EXECUTABLE_NAME" ]; then
        echo "Deleting the executable $EXECUTABLE_NAME..."
        rm -f "./$EXECUTABLE_NAME"
    fi
}
trap cleanup EXIT

NUM_THREADS="$(awk -F= '/^num_threads=/ {print $2; exit}' "$CONFIG_FILE")"
export OMP_NUM_THREADS="${NUM_THREADS:-1}"
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "Compiling $SOURCE_FILE..."
gcc -fopenmp "$SOURCE_FILE" $LIB_FILE -o "$EXECUTABLE_NAME" -lm
echo "Compilation successful. Executable: $EXECUTABLE_NAME"

# === Running the simulation ===
echo "Running the simulation with OMP_NUM_THREADS=$OMP_NUM_THREADS..."
"./$EXECUTABLE_NAME" "$CONFIG_FILE"
echo "Simulation completed successfully."
