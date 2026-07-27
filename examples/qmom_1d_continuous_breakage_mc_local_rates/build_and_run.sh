#!/usr/bin/env bash
set -e

script_dir="$(cd "$(dirname "$0")" && pwd)"
build_dir="${script_dir}/build"
target="qmom_1d_continuous_breakage_mc_local_rates_example"
prefix_path="${1:-${CMAKE_PREFIX_PATH:-}}"

echo "Building QMOM 1D local moment-rate example..."
if [[ -n "${prefix_path}" ]]; then
    cmake -S "${script_dir}" -B "${build_dir}" -DCMAKE_BUILD_TYPE=Debug -DCMAKE_PREFIX_PATH="${prefix_path}"
else
    cmake -S "${script_dir}" -B "${build_dir}" -DCMAKE_BUILD_TYPE=Debug
fi
cmake --build "${build_dir}"

echo "Running example..."
(
    cd "${build_dir}"
    "./${target}"
)

echo "Outputs written to:"
echo "${build_dir}/output"
