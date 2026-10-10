#!/usr/bin/env python3
"""Self-test for the Monte Carlo moment-history comparison utility."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path


script_path = Path(__file__).with_name("compare_moments.py")
spec = importlib.util.spec_from_file_location("compare_moments", script_path)
if spec is None or spec.loader is None:
    raise RuntimeError("Could not import compare_moments.py")
compare_moments = importlib.util.module_from_spec(spec)
sys.modules["compare_moments"] = compare_moments
spec.loader.exec_module(compare_moments)


def check(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def write_text(path: Path, header: list[str], rows: list[list[float | int]]) -> None:
    with path.open("w") as handle:
        handle.write("# " + " ".join(header) + "\n")
        for row in rows:
            handle.write(" ".join(str(value) for value in row) + "\n")


def read_rows(path: Path) -> list[dict[str, str]]:
    fieldnames: list[str] | None = None
    rows: list[dict[str, str]] = []
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                candidate = stripped[1:].strip().split()
                if candidate and candidate[0] == "time":
                    fieldnames = candidate
                continue
            check(fieldnames is not None, f"{path} did not provide a data header before rows")
            rows.append(dict(zip(fieldnames, stripped.split())))
    return rows


def run_compare(reference: Path, comparison: Path, output_dir: Path, statistics: Path | None = None, moments: list[str] | None = None) -> list[dict[str, str]]:
    argv = [str(reference), str(comparison), "--output-dir", str(output_dir)]
    if statistics is not None:
        argv += ["--mc-statistics", str(statistics)]
    if moments is not None:
        argv += ["--moments", *moments]
    status = compare_moments.main(argv)
    check(status == 0, f"compare_moments.py returned {status}")
    return read_rows(output_dir / "moment_comparison.txt")


def test_identical_files(workdir: Path) -> None:
    path = workdir / "identical.txt"
    write_text(path, ["time", "M0", "M1"], [[0.0, 2.0, 3.0], [0.5, 2.0, 4.0]])
    rows = run_compare(path, path, workdir / "identical_out")
    check(all(float(row["signed_error"]) == 0.0 for row in rows), "identical files produced nonzero signed error")
    check(all(float(row["absolute_error"]) == 0.0 for row in rows), "identical files produced nonzero absolute error")
    check(all(float(row["relative_error"]) == 0.0 for row in rows), "identical files produced nonzero relative error")


def test_perturbed_moments(workdir: Path) -> None:
    reference = workdir / "reference.txt"
    comparison = workdir / "comparison.txt"
    write_text(reference, ["time", "M0", "M1"], [[0.0, 2.0, 3.0], [0.5, 2.0, 4.0]])
    write_text(comparison, ["time", "M0", "M1"], [[0.0, 2.0, 3.5], [0.5, 2.1, 3.0]])
    rows = run_compare(reference, comparison, workdir / "perturbed_out")
    check(any(float(row["signed_error"]) > 0.0 for row in rows), "perturbed files did not produce a positive signed error")
    check(any(float(row["signed_error"]) < 0.0 for row in rows), "perturbed files did not produce a negative signed error")
    check(all(float(row["absolute_error"]) >= 0.0 for row in rows), "absolute error column contained a negative value")
    check(
        any(row["moment"] == "M1" and abs(float(row["signed_error"]) - 0.5) < 1.0e-14 for row in rows),
        "M1 signed error at t=0.0 was not 0.5",
    )
    check(
        any(row["moment"] == "M1" and abs(float(row["absolute_error"]) - 1.0) < 1.0e-14 for row in rows),
        "M1 absolute error at t=0.5 was not 1.0",
    )


def test_mc_statistics(workdir: Path) -> None:
    reference = workdir / "stats_reference.txt"
    comparison = workdir / "stats_comparison.txt"
    statistics = workdir / "statistics.txt"
    write_text(reference, ["time", "M0", "M1", "M2"], [[0.0, 2.0, 3.0, 0.01]])
    write_text(comparison, ["time", "M0", "M1", "M2"], [[0.0, 2.0, 3.2, 0.03]])
    write_text(statistics, ["time", "sample_count", "seed", "number_density", "SE_M0", "SE_M1", "SE_M2"], [[0.0, 10, 7, 2.0, 0.0, 0.1, 0.05]])
    rows = run_compare(reference, comparison, workdir / "stats_out", statistics=statistics)
    m1 = next(row for row in rows if row["moment"] == "M1")
    check(abs(float(m1["z_score"]) - 2.0) < 1.0e-14, "M1 z-score was not 2.0")
    check(abs(float(m1["relative_error"]) - 0.2 / 3.0) < 1.0e-14, "M1 relative error was not measured against |reference|")
    # A near-zero reference below its standard error is measured against the
    # standard error, as in compare_moments(...) of the C++ package.
    m2 = next(row for row in rows if row["moment"] == "M2")
    check(abs(float(m2["relative_error"]) - 0.02 / 0.05) < 1.0e-14, "M2 relative error was not measured against the standard error")


def test_multidimensional_labels(workdir: Path) -> None:
    reference = workdir / "tensor_reference.txt"
    comparison = workdir / "tensor_comparison.txt"
    write_text(reference, ["time", "M00", "M01", "M10"], [[0.0, 2.0, 0.5, 1.5]])
    write_text(comparison, ["time", "M00", "M01", "M10"], [[0.0, 2.0, 0.25, 1.75]])
    rows = run_compare(reference, comparison, workdir / "tensor_out", moments=["M00", "M01", "M10"])
    check([row["moment"] for row in rows] == ["M00", "M01", "M10"], "multidimensional moment labels were not preserved")
    check(any(row["moment"] == "M01" and float(row["signed_error"]) == -0.25 for row in rows), "M01 signed error was not -0.25")
    check(any(row["moment"] == "M01" and float(row["absolute_error"]) == 0.25 for row in rows), "M01 absolute error was not 0.25")


def main() -> int:
    with tempfile.TemporaryDirectory() as directory:
        workdir = Path(directory)
        test_identical_files(workdir)
        test_perturbed_moments(workdir)
        test_mc_statistics(workdir)
        test_multidimensional_labels(workdir)
    print("compare_moments.py self-test passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
