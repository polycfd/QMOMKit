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
            assert fieldnames is not None
            rows.append(dict(zip(fieldnames, stripped.split())))
    return rows


def run_compare(reference: Path, comparison: Path, output_dir: Path, statistics: Path | None = None, moments: list[str] | None = None) -> list[dict[str, str]]:
    argv = [str(reference), str(comparison), "--output-dir", str(output_dir)]
    if statistics is not None:
        argv += ["--mc-statistics", str(statistics)]
    if moments is not None:
        argv += ["--moments", *moments]
    status = compare_moments.main(argv)
    assert status == 0
    return read_rows(output_dir / "moment_comparison.txt")


def test_identical_files(workdir: Path) -> None:
    path = workdir / "identical.txt"
    write_text(path, ["time", "M0", "M1"], [[0.0, 2.0, 3.0], [0.5, 2.0, 4.0]])
    rows = run_compare(path, path, workdir / "identical_out")
    assert all(float(row["signed_error"]) == 0.0 for row in rows)
    assert all(float(row["absolute_error"]) == 0.0 for row in rows)
    assert all(float(row["relative_error"]) == 0.0 for row in rows)


def test_perturbed_moments(workdir: Path) -> None:
    reference = workdir / "reference.txt"
    comparison = workdir / "comparison.txt"
    write_text(reference, ["time", "M0", "M1"], [[0.0, 2.0, 3.0], [0.5, 2.0, 4.0]])
    write_text(comparison, ["time", "M0", "M1"], [[0.0, 2.0, 3.5], [0.5, 2.1, 3.0]])
    rows = run_compare(reference, comparison, workdir / "perturbed_out")
    assert any(float(row["signed_error"]) > 0.0 for row in rows)
    assert any(float(row["signed_error"]) < 0.0 for row in rows)
    assert all(float(row["absolute_error"]) >= 0.0 for row in rows)
    assert any(row["moment"] == "M1" and abs(float(row["signed_error"]) - 0.5) < 1.0e-14 for row in rows)
    assert any(row["moment"] == "M1" and abs(float(row["absolute_error"]) - 1.0) < 1.0e-14 for row in rows)


def test_mc_statistics(workdir: Path) -> None:
    reference = workdir / "stats_reference.txt"
    comparison = workdir / "stats_comparison.txt"
    statistics = workdir / "statistics.txt"
    write_text(reference, ["time", "M0", "M1"], [[0.0, 2.0, 3.0]])
    write_text(comparison, ["time", "M0", "M1"], [[0.0, 2.0, 3.2]])
    write_text(statistics, ["time", "sample_count", "seed", "number_density", "SE_M0", "SE_M1"], [[0.0, 10, 7, 2.0, 0.0, 0.1]])
    rows = run_compare(reference, comparison, workdir / "stats_out", statistics=statistics)
    m1 = next(row for row in rows if row["moment"] == "M1")
    assert abs(float(m1["z_score"]) - 2.0) < 1.0e-14


def test_multidimensional_labels(workdir: Path) -> None:
    reference = workdir / "tensor_reference.txt"
    comparison = workdir / "tensor_comparison.txt"
    write_text(reference, ["time", "M00", "M01", "M10"], [[0.0, 2.0, 0.5, 1.5]])
    write_text(comparison, ["time", "M00", "M01", "M10"], [[0.0, 2.0, 0.25, 1.75]])
    rows = run_compare(reference, comparison, workdir / "tensor_out", moments=["M00", "M01", "M10"])
    assert [row["moment"] for row in rows] == ["M00", "M01", "M10"]
    assert any(row["moment"] == "M01" and float(row["signed_error"]) == -0.25 for row in rows)
    assert any(row["moment"] == "M01" and float(row["absolute_error"]) == 0.25 for row in rows)


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
