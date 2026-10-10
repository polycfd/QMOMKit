#!/usr/bin/env python3
"""Compare two QMOMKit/Monte Carlo moment-history text tables.

The first table is the reference, such as a QMOMKit history, and the second the
Monte Carlo estimate. The errors follow `compare_moments(...)` of the C++
package: signed_error = second - first, and the relative error divides by the
largest of |reference|, the Monte Carlo standard error and a positive floor.
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path


default_epsilon = 1.0e-300
metadata_columns = {"time", "sample_count", "seed", "number_density"}


@dataclass(frozen=True)
class MomentHistory:
    path: Path
    fieldnames: list[str]
    rows: list[dict[str, str]]


@dataclass(frozen=True)
class StatisticsHistory:
    path: Path
    fieldnames: list[str]
    rows: list[dict[str, str]]


@dataclass(frozen=True)
class ComparisonPoint:
    time: float
    moment: str
    reference: float
    comparison: float
    signed_error: float
    absolute_error: float
    relative_error: float
    z_score: float | None


def is_moment_column(name: str) -> bool:
    return name.startswith("M") and name not in metadata_columns and not name.startswith("SE_")


def read_text_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    fieldnames: list[str] | None = None
    rows: list[dict[str, str]] = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                candidate = stripped[1:].strip().split()
                if candidate and candidate[0] == "time":
                    fieldnames = candidate
                continue
            if fieldnames is None:
                raise ValueError(f"{path} data row appears before a '# time ...' header at line {line_number}")
            values = stripped.split()
            if len(values) != len(fieldnames):
                raise ValueError(
                    f"{path} line {line_number} has {len(values)} columns, but the header has {len(fieldnames)}"
                )
            rows.append(dict(zip(fieldnames, values)))

    if fieldnames is None:
        raise ValueError(f"{path} does not contain a commented '# time ...' header")
    if "time" not in fieldnames:
        raise ValueError(f"{path} does not contain a 'time' column")
    if not rows:
        raise ValueError(f"{path} does not contain any data rows")
    return fieldnames, rows


def read_moment_history(path: Path) -> MomentHistory:
    fieldnames, rows = read_text_rows(path)
    moments = moment_columns(fieldnames)
    if not moments:
        raise ValueError(f"{path} does not contain any moment columns")
    return MomentHistory(path=path, fieldnames=fieldnames, rows=rows)


def read_statistics_history(path: Path) -> StatisticsHistory:
    fieldnames, rows = read_text_rows(path)
    if not any(name.startswith("SE_") for name in fieldnames):
        raise ValueError(f"{path} does not contain any standard-error columns named SE_<moment>")
    return StatisticsHistory(path=path, fieldnames=fieldnames, rows=rows)


def moment_columns(fieldnames: list[str]) -> list[str]:
    return [name for name in fieldnames if is_moment_column(name)]


def parse_float(row: dict[str, str], column: str, path: Path) -> float:
    try:
        return float(row[column])
    except KeyError as error:
        raise ValueError(f"{path} does not contain column {column!r}") from error
    except ValueError as error:
        raise ValueError(f"{path} column {column!r} contains a non-numeric value: {row[column]!r}") from error


def validate_time_grids(reference: MomentHistory, comparison: MomentHistory) -> list[float]:
    if len(reference.rows) != len(comparison.rows):
        raise ValueError(
            f"text tables have different row counts: {reference.path} has {len(reference.rows)}, "
            f"{comparison.path} has {len(comparison.rows)}"
        )

    times: list[float] = []
    for row_index, (reference_row, comparison_row) in enumerate(zip(reference.rows, comparison.rows)):
        reference_time = parse_float(reference_row, "time", reference.path)
        comparison_time = parse_float(comparison_row, "time", comparison.path)
        if not math.isclose(reference_time, comparison_time, rel_tol=1.0e-12, abs_tol=1.0e-12):
            raise ValueError(
                f"Time grids differ at row {row_index}: reference time {reference_time:.17g}, "
                f"comparison time {comparison_time:.17g}"
            )
        times.append(reference_time)
    return times


def validate_statistics_grid(statistics: StatisticsHistory, times: list[float]) -> None:
    if len(statistics.rows) != len(times):
        raise ValueError(
            f"Statistics text table has {len(statistics.rows)} rows, but the moment histories have {len(times)} rows"
        )

    for row_index, (row, expected_time) in enumerate(zip(statistics.rows, times)):
        statistics_time = parse_float(row, "time", statistics.path)
        if not math.isclose(statistics_time, expected_time, rel_tol=1.0e-12, abs_tol=1.0e-12):
            raise ValueError(
                f"Statistics time grid differs at row {row_index}: expected {expected_time:.17g}, "
                f"got {statistics_time:.17g}"
            )


def select_moments(reference: MomentHistory, comparison: MomentHistory, requested: list[str] | None) -> list[str]:
    reference_moments = moment_columns(reference.fieldnames)
    comparison_moments = moment_columns(comparison.fieldnames)

    if requested is None:
        if reference_moments != comparison_moments:
            reference_set = set(reference_moments)
            comparison_set = set(comparison_moments)
            missing_from_comparison = [moment for moment in reference_moments if moment not in comparison_set]
            missing_from_reference = [moment for moment in comparison_moments if moment not in reference_set]
            messages = []
            if missing_from_comparison:
                messages.append("comparison missing " + ", ".join(missing_from_comparison))
            if missing_from_reference:
                messages.append("reference missing " + ", ".join(missing_from_reference))
            if not messages:
                messages.append("moment columns appear in different orders")
            raise ValueError("Moment columns do not match: " + "; ".join(messages))
        return reference_moments

    missing_reference = [moment for moment in requested if moment not in reference_moments]
    missing_comparison = [moment for moment in requested if moment not in comparison_moments]
    messages = []
    if missing_reference:
        messages.append("reference missing " + ", ".join(missing_reference))
    if missing_comparison:
        messages.append("comparison missing " + ", ".join(missing_comparison))
    if messages:
        raise ValueError("Requested moment columns are unavailable: " + "; ".join(messages))
    return requested


def standard_error(statistics_row: dict[str, str], moment: str, statistics_path: Path) -> float | None:
    column = "SE_" + moment
    value = statistics_row.get(column)
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError as error:
        raise ValueError(f"{statistics_path} column {column!r} contains a non-numeric value: {value!r}") from error


def compute_z_score(error: float, error_scale: float | None) -> float | None:
    if error_scale is None:
        return None
    if error_scale == 0.0:
        return 0.0 if error == 0.0 else math.copysign(math.inf, error)
    return error / error_scale


def compare_histories(
    reference: MomentHistory,
    comparison: MomentHistory,
    moments: list[str],
    epsilon: float,
    statistics: StatisticsHistory | None = None,
) -> list[ComparisonPoint]:
    if epsilon <= 0.0 or not math.isfinite(epsilon):
        raise ValueError("epsilon must be finite and positive")

    times = validate_time_grids(reference, comparison)
    if statistics is not None:
        validate_statistics_grid(statistics, times)

    points: list[ComparisonPoint] = []
    for row_index, (time, reference_row, comparison_row) in enumerate(zip(times, reference.rows, comparison.rows)):
        statistics_row = None if statistics is None else statistics.rows[row_index]
        for moment in moments:
            reference_value = parse_float(reference_row, moment, reference.path)
            comparison_value = parse_float(comparison_row, moment, comparison.path)
            signed_error = comparison_value - reference_value
            absolute_error = abs(signed_error)
            error_scale = None if statistics_row is None else standard_error(statistics_row, moment, statistics.path)
            # Near-zero moments with visible sampling noise are measured against
            # that noise rather than against their own magnitude.
            relative_error = absolute_error / max(abs(reference_value), 0.0 if error_scale is None else error_scale, epsilon)
            points.append(
                ComparisonPoint(
                    time=time,
                    moment=moment,
                    reference=reference_value,
                    comparison=comparison_value,
                    signed_error=signed_error,
                    absolute_error=absolute_error,
                    relative_error=relative_error,
                    z_score=compute_z_score(signed_error, error_scale),
                )
            )
    return points


def write_comparison_text(path: Path, points: list[ComparisonPoint], include_z_score: bool) -> None:
    with path.open("w") as handle:
        header = ["time", "moment", "reference", "computed", "signed_error", "absolute_error", "relative_error"]
        if include_z_score:
            header.append("z_score")
        handle.write("# " + " ".join(header) + "\n")
        for point in points:
            row = [
                f"{point.time:.17g}",
                point.moment,
                f"{point.reference:.17g}",
                f"{point.comparison:.17g}",
                f"{point.signed_error:.17g}",
                f"{point.absolute_error:.17g}",
                f"{point.relative_error:.17g}",
            ]
            if include_z_score:
                row.append("" if point.z_score is None else f"{point.z_score:.17g}")
            handle.write(" ".join(row) + "\n")


def print_summary(points: list[ComparisonPoint], moments: list[str], include_z_score: bool) -> None:
    header = "# moment max_absolute_error max_relative_error final_signed_error final_absolute_error final_relative_error"
    if include_z_score:
        header += " max_abs_z_score final_z_score"
    print(header)

    for moment in moments:
        moment_points = [point for point in points if point.moment == moment]
        max_absolute = max(abs(point.absolute_error) for point in moment_points)
        max_relative = max(point.relative_error for point in moment_points)
        final = moment_points[-1]
        row = (
            f"{moment} {max_absolute:.6g} {max_relative:.6g} "
            f"{final.signed_error:.6g} {final.absolute_error:.6g} {final.relative_error:.6g}"
        )
        if include_z_score:
            z_values = [abs(point.z_score) for point in moment_points if point.z_score is not None]
            max_abs_z = "" if not z_values else f"{max(z_values):.6g}"
            final_z = "" if final.z_score is None else f"{final.z_score:.6g}"
            row += f" {max_abs_z} {final_z}"
        print(row)


def plot_metric(output_dir: Path, points: list[ComparisonPoint], moments: list[str], metric: str, ylabel: str, filename: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    for moment in moments:
        moment_points = [point for point in points if point.moment == moment]
        values = [getattr(point, metric) for point in moment_points]
        ax.plot([point.time for point in moment_points], values, marker="o", markersize=2.5, linewidth=1.4, label=moment)
    ax.set_xlabel("time")
    ax.set_ylabel(ylabel)
    if metric != "z_score" and any(getattr(point, metric) > 0.0 for point in points):
        ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", ncols=1 if len(moments) <= 8 else 2)
    fig.savefig(output_dir / filename, bbox_inches="tight")
    plt.close(fig)


def write_plots(output_dir: Path, points: list[ComparisonPoint], moments: list[str], include_z_score: bool) -> None:
    try:
        plot_metric(output_dir, points, moments, "absolute_error", "absolute error", "absolute_error.pdf")
        plot_metric(output_dir, points, moments, "relative_error", "relative error", "relative_error.pdf")
        if include_z_score:
            z_points = [point for point in points if point.z_score is not None and math.isfinite(point.z_score)]
            if z_points:
                plot_metric(output_dir, z_points, moments, "z_score", "z-score", "z_score.pdf")
    except ImportError as error:
        raise RuntimeError("matplotlib is required for --plots") from error


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare two raw moment-history text tables.")
    parser.add_argument("reference_table", type=Path, help="Reference moment-history text table, such as a QMOMKit history.")
    parser.add_argument("comparison_table", type=Path, help="Monte Carlo moment-history text table.")
    parser.add_argument("--mc-statistics", type=Path, help="Optional Monte Carlo statistics text table with SE_<moment> columns.")
    parser.add_argument("--output-dir", type=Path, default=Path("."), help="Directory for moment_comparison.txt and optional plots.")
    parser.add_argument("--moments", nargs="+", help="Moment columns to compare, for example M0 M2 M5 or M00 M01 M10.")
    parser.add_argument("--epsilon", type=float, default=default_epsilon, help="Positive denominator floor for relative error.")
    parser.add_argument("--plots", action="store_true", help="Write matplotlib PDF plots of the error histories.")
    return parser.parse_args(argv)


def run(args: argparse.Namespace) -> int:
    reference = read_moment_history(args.reference_table)
    comparison = read_moment_history(args.comparison_table)
    statistics = None if args.mc_statistics is None else read_statistics_history(args.mc_statistics)
    moments = select_moments(reference, comparison, args.moments)
    points = compare_histories(reference, comparison, moments, args.epsilon, statistics)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    comparison_path = args.output_dir / "moment_comparison.txt"
    include_z_score = statistics is not None
    write_comparison_text(comparison_path, points, include_z_score)
    if args.plots:
        write_plots(args.output_dir, points, moments, include_z_score)

    print_summary(points, moments, include_z_score)
    print(f"Wrote comparison rows to {comparison_path}", file=sys.stderr)
    return 0


def main(argv: list[str] | None = None) -> int:
    try:
        return run(parse_args(argv))
    except Exception as error:
        print(f"compare_moments.py: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
