#include "qmomkit/realizability.hpp"

#include "internal/linear_algebra.hpp"
#include "internal/numeric_utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace qmomkit
{

bool RealizabilityResult::realizable() const noexcept
{
    return status == RealizabilityStatus::strictly_realizable || status == RealizabilityStatus::boundary;
}

bool RealizabilityResult::on_boundary() const noexcept
{
    return status == RealizabilityStatus::boundary;
}

namespace
{

using internal::DenseMatrix;
using internal::EigenDecomposition;

enum class MatrixStatus
{
    positive_definite,
    positive_semidefinite,
    not_positive_semidefinite,
    numerically_indeterminate
};

struct MatrixAnalysis
{
    std::string name;
    std::size_t order = 0;
    MatrixStatus status = MatrixStatus::not_positive_semidefinite;
    double minimum_eigenvalue = 0.0;
    double tolerance = 0.0;
    EigenDecomposition eigensystem{std::vector<double>{}, DenseMatrix(0, 0)};
};

struct ClassificationState
{
    RealizabilityResult result;
    bool boundary = false;
};

[[nodiscard]] const char* support_name(const MomentSupport support)
{
    switch (support)
    {
    case MomentSupport::real_line:
        return "real_line";
    case MomentSupport::nonnegative:
        return "nonnegative";
    case MomentSupport::bounded_interval:
        return "bounded_interval";
    }
    return "unknown";
}

[[nodiscard]] RealizabilityResult make_result(const RealizabilityStatus status, const std::size_t checked_order, const std::string& matrix_name,
                                              const std::size_t matrix_order, const double minimum_pivot, const double tolerance, std::string message)
{
    return RealizabilityResult{status, checked_order, matrix_order, matrix_name, minimum_pivot, tolerance, std::move(message)};
}

[[nodiscard]] double vector_scale(const std::vector<double>& values)
{
    double scale = 0.0;
    for (const double value : values)
    {
        scale = std::max(scale, std::abs(value));
    }
    return scale;
}

[[nodiscard]] EigenDecomposition symmetric_jacobi_eigen(DenseMatrix matrix, const double matrix_scale)
{
    const std::size_t n = matrix.rows();
    constexpr std::size_t max_sweeps_per_node = 80;
    constexpr double diagonalization_tolerance = 1.0e-14;
    const std::size_t max_sweeps = max_sweeps_per_node * std::max<std::size_t>(n, 1);

    if (matrix_scale == 0.0 || !std::isfinite(matrix_scale))
    {
        return internal::symmetric_jacobi_eigen(std::move(matrix), {diagonalization_tolerance, max_sweeps});
    }

    // The Jacobi convergence tolerance is absolute. Scale this temporary
    // matrix so its convergence test is independent of the physical moment
    // magnitude, then restore the original eigenvalue scale for all
    // realizability comparisons. Positive scalar scaling leaves eigenvectors
    // unchanged.
    for (std::size_t row = 0; row < matrix.rows(); ++row)
    {
        for (std::size_t col = 0; col < matrix.cols(); ++col)
        {
            matrix(row, col) /= matrix_scale;
        }
    }

    EigenDecomposition eigensystem = internal::symmetric_jacobi_eigen(std::move(matrix), {diagonalization_tolerance, max_sweeps});
    for (double& eigenvalue : eigensystem.eigenvalues)
    {
        eigenvalue *= matrix_scale;
    }
    return eigensystem;
}

[[nodiscard]] MatrixAnalysis analyze_matrix(std::string name, const std::size_t order, const DenseMatrix& matrix, const RealizabilityOptions& options)
{
    MatrixAnalysis analysis{std::move(name), order};
    double normalization_scale = 0.0;
    for (std::size_t row = 0; row < matrix.rows(); ++row)
    {
        for (std::size_t col = 0; col < matrix.cols(); ++col)
        {
            const double value = std::abs(matrix(row, col));
            if (!std::isfinite(value))
            {
                analysis.status = MatrixStatus::numerically_indeterminate;
                return analysis;
            }
            normalization_scale = std::max(normalization_scale, value);
        }
    }

    double normalized_infinity_norm = 0.0;
    if (normalization_scale > 0.0)
    {
        for (std::size_t row = 0; row < matrix.rows(); ++row)
        {
            double row_sum = 0.0;
            for (std::size_t col = 0; col < matrix.cols(); ++col)
            {
                row_sum += std::abs(matrix(row, col)) / normalization_scale;
            }
            normalized_infinity_norm = std::max(normalized_infinity_norm, row_sum);
        }
    }

    analysis.tolerance = options.relative_tolerance ? options.tolerance * normalization_scale * normalized_infinity_norm : options.tolerance;
    analysis.eigensystem = symmetric_jacobi_eigen(matrix, normalization_scale);
    analysis.minimum_eigenvalue = analysis.eigensystem.eigenvalues.empty() ? 0.0 : analysis.eigensystem.eigenvalues.front();

    // Positive semidefinite matrices indicate realizability on the boundary.
    // Eigenvalues within the user tolerance but outside the tighter roundoff
    // band are reported as indeterminate.
    const double diagonalization_roundoff = 1000.0 * std::numeric_limits<double>::epsilon() * normalization_scale * normalized_infinity_norm;

    if (!internal::all_finite(analysis.eigensystem.eigenvalues))
    {
        analysis.status = MatrixStatus::numerically_indeterminate;
    }
    else if (analysis.minimum_eigenvalue < -analysis.tolerance)
    {
        analysis.status = MatrixStatus::not_positive_semidefinite;
    }
    else if (std::abs(analysis.minimum_eigenvalue) <= diagonalization_roundoff)
    {
        analysis.status = MatrixStatus::positive_semidefinite;
    }
    else if (analysis.minimum_eigenvalue < analysis.tolerance)
    {
        analysis.status = MatrixStatus::numerically_indeterminate;
    }
    else
    {
        analysis.status = MatrixStatus::positive_definite;
    }

    return analysis;
}

void update_limiting_matrix(ClassificationState& state, const MatrixAnalysis& analysis)
{
    const double current_normalized =
        std::abs(state.result.minimum_pivot) / std::max(state.result.tolerance_used, std::numeric_limits<double>::min());
    const double candidate_normalized = std::abs(analysis.minimum_eigenvalue) / std::max(analysis.tolerance, std::numeric_limits<double>::min());
    if (state.result.failed_matrix.empty() || candidate_normalized < current_normalized)
    {
        state.result.failed_matrix = analysis.name;
        state.result.failed_matrix_order = analysis.order;
        state.result.minimum_pivot = analysis.minimum_eigenvalue;
        state.result.tolerance_used = analysis.tolerance;
    }
}

[[nodiscard]] bool incorporate_matrix(ClassificationState& state, const MatrixAnalysis& analysis, const std::size_t checked_order)
{
    update_limiting_matrix(state, analysis);

    if (analysis.status == MatrixStatus::not_positive_semidefinite)
    {
        state.result =
            make_result(RealizabilityStatus::not_realizable, checked_order, analysis.name, analysis.order, analysis.minimum_eigenvalue,
                        analysis.tolerance, "Moment sequence is not realizable: " + analysis.name + " has a negative eigenvalue beyond tolerance.");
        return false;
    }

    if (analysis.status == MatrixStatus::numerically_indeterminate)
    {
        state.result =
            make_result(RealizabilityStatus::numerically_indeterminate, checked_order, analysis.name, analysis.order, analysis.minimum_eigenvalue,
                        analysis.tolerance,
                        "Realizability is numerically indeterminate: " + analysis.name + " has an eigenvalue within the numerical tolerance.");
        return false;
    }

    if (analysis.status == MatrixStatus::positive_semidefinite)
    {
        state.boundary = true;
    }
    return true;
}

[[nodiscard]] DenseMatrix hankel_matrix(const MomentSet& moments, const std::size_t order, const std::size_t shift)
{
    // Hankel moment matrices encode nonnegativity of squared polynomials on
    // the support; shifted variants enforce support-specific constraints.
    DenseMatrix matrix(order + 1, order + 1);
    for (std::size_t row = 0; row <= order; ++row)
    {
        for (std::size_t col = 0; col <= order; ++col)
        {
            matrix(row, col) = moments[row + col + shift];
        }
    }
    return matrix;
}

[[nodiscard]] DenseMatrix interval_quadratic_localizing_matrix(const MomentSet& moments, const BoundedInterval interval, const std::size_t order)
{
    // Localizes by (x - lower)(upper - x), enforcing support inside the bounded
    // interval for even-order Hausdorff checks.
    DenseMatrix matrix(order + 1, order + 1);
    for (std::size_t row = 0; row <= order; ++row)
    {
        for (std::size_t col = 0; col <= order; ++col)
        {
            const std::size_t index = row + col;
            matrix(row, col) =
                (interval.lower + interval.upper) * moments[index + 1] - interval.lower * interval.upper * moments[index] - moments[index + 2];
        }
    }
    return matrix;
}

[[nodiscard]] DenseMatrix interval_linear_localizing_matrix(const MomentSet& moments, const BoundedInterval interval, const std::size_t order,
                                                            const bool lower_bound)
{
    // Linear localizing matrices encode nonnegativity of x - lower and
    // upper - x for odd-order bounded-interval checks.
    DenseMatrix matrix(order + 1, order + 1);
    for (std::size_t row = 0; row <= order; ++row)
    {
        for (std::size_t col = 0; col <= order; ++col)
        {
            const std::size_t index = row + col;
            matrix(row, col) =
                lower_bound ? moments[index + 1] - interval.lower * moments[index] : interval.upper * moments[index] - moments[index + 1];
        }
    }
    return matrix;
}

[[nodiscard]] std::vector<double> extension_column(const MomentSet& moments, const std::size_t matrix_order, const std::size_t shift)
{
    std::vector<double> values(matrix_order + 1, 0.0);
    for (std::size_t row = 0; row <= matrix_order; ++row)
    {
        values[row] = moments[row + matrix_order + 1 + shift];
    }
    return values;
}

[[nodiscard]] double projection_on_eigenvector(const DenseMatrix& eigenvectors, const std::size_t column, const std::vector<double>& vector)
{
    double projection = 0.0;
    for (std::size_t row = 0; row < vector.size(); ++row)
    {
        projection += eigenvectors(row, column) * vector[row];
    }
    return projection;
}

[[nodiscard]] RealizabilityStatus check_range_consistency(const MatrixAnalysis& analysis, const std::vector<double>& vector,
                                                          const RealizabilityOptions& options, double& maximum_kernel_projection)
{
    // For a singular moment/localizing matrix, the next extension column must
    // lie in its range. Projection onto nullspace eigenvectors detects failure.
    maximum_kernel_projection = 0.0;
    const double tolerance = options.relative_tolerance ? std::max(analysis.tolerance, options.tolerance * vector_scale(vector)) : options.tolerance;

    for (std::size_t col = 0; col < analysis.eigensystem.eigenvalues.size(); ++col)
    {
        if (analysis.eigensystem.eigenvalues[col] <= analysis.tolerance)
        {
            const double projection = std::abs(projection_on_eigenvector(analysis.eigensystem.eigenvectors, col, vector));
            maximum_kernel_projection = std::max(maximum_kernel_projection, projection);
        }
    }

    if (maximum_kernel_projection > tolerance)
    {
        return RealizabilityStatus::not_realizable;
    }
    if (maximum_kernel_projection > 0.1 * tolerance)
    {
        return RealizabilityStatus::numerically_indeterminate;
    }
    return RealizabilityStatus::boundary;
}

[[nodiscard]] bool incorporate_range_consistency(ClassificationState& state, const MatrixAnalysis& analysis, const std::vector<double>& vector,
                                                 const std::string& name, const std::size_t checked_order, const RealizabilityOptions& options)
{
    if (analysis.status == MatrixStatus::positive_definite)
    {
        return true;
    }

    double maximum_kernel_projection = 0.0;
    const RealizabilityStatus status = check_range_consistency(analysis, vector, options, maximum_kernel_projection);
    state.boundary = true;
    if (status == RealizabilityStatus::boundary)
    {
        return true;
    }

    const double tolerance = options.relative_tolerance ? std::max(analysis.tolerance, options.tolerance * vector_scale(vector)) : options.tolerance;
    state.result =
        make_result(status, checked_order, name, analysis.order, maximum_kernel_projection, tolerance,
                    status == RealizabilityStatus::not_realizable
                        ? "Moment sequence is not realizable: " + name + " violates singular-matrix kernel consistency."
                        : "Realizability is numerically indeterminate: " + name + " is too close to a singular kernel-consistency boundary.");
    return false;
}

[[nodiscard]] std::size_t numerical_rank(const MatrixAnalysis& analysis)
{
    return static_cast<std::size_t>(std::count_if(analysis.eigensystem.eigenvalues.begin(), analysis.eigensystem.eigenvalues.end(),
                                                  [&analysis](const double eigenvalue) { return eigenvalue > analysis.tolerance; }));
}

[[nodiscard]] RealizabilityResult finish_classification(ClassificationState state, const std::size_t checked_order, const MomentSupport support)
{
    state.result.checked_order = checked_order;
    if (state.boundary)
    {
        state.result.status = RealizabilityStatus::boundary;
        state.result.message = std::string("Moment sequence is boundary-realizable on support ") + support_name(support) + ".";
    }
    else
    {
        state.result.status = RealizabilityStatus::strictly_realizable;
        state.result.message = std::string("Moment sequence is strictly realizable on support ") + support_name(support) + ".";
    }
    return state.result;
}

[[nodiscard]] RealizabilityResult validate_common_inputs(const MomentSet& moments, const RealizabilityOptions& options)
{
    if (moments.empty())
    {
        return make_result(RealizabilityStatus::invalid_input, 0, "", 0, 0.0, 0.0, "Moment sequence is empty.");
    }
    if (!internal::all_finite(moments.values()))
    {
        return make_result(RealizabilityStatus::invalid_input, moments.size() - 1, "", 0, 0.0, 0.0, "Moment sequence contains a non-finite value.");
    }
    if (!std::isfinite(options.tolerance) || options.tolerance <= 0.0)
    {
        return make_result(RealizabilityStatus::invalid_input, moments.size() - 1, "", 0, 0.0, 0.0,
                           "Realizability tolerance must be finite and positive.");
    }
    if (moments[0] < 0.0)
    {
        return make_result(RealizabilityStatus::not_realizable, moments.size() - 1, "m0", 0, moments[0], options.tolerance,
                           "Moment sequence is not realizable: m0 is negative.");
    }
    if (moments[0] == 0.0)
    {
        const bool all_zero = std::all_of(moments.values().begin(), moments.values().end(), [](const double value) { return value == 0.0; });
        if (all_zero)
        {
            return make_result(RealizabilityStatus::boundary, moments.size() - 1, "zero_measure", 0, 0.0, options.tolerance,
                               "All moments are zero; the zero measure is boundary-realizable.");
        }
        return make_result(RealizabilityStatus::not_realizable, moments.size() - 1, "m0", 0, moments[0], options.tolerance,
                           "Moment sequence is not realizable: m0 is zero but at least one higher moment is nonzero.");
    }

    return RealizabilityResult{RealizabilityStatus::strictly_realizable, moments.size() - 1, 0, "", 0.0, options.tolerance, ""};
}

[[nodiscard]] RealizabilityResult check_hamburger(const MomentSet& moments, const RealizabilityOptions& options)
{
    // Hamburger realizability on the real line is governed by PSD Hankel moment
    // matrices. A singular even truncation must be a flat extension of its
    // preceding moment matrix; odd truncations require next-column range validity.
    const std::size_t checked_order = moments.size() - 1;
    const std::size_t hankel_order = checked_order / 2;
    ClassificationState state;
    state.result.checked_order = checked_order;

    MatrixAnalysis moment_matrix =
        analyze_matrix("H_" + std::to_string(hankel_order), hankel_order, hankel_matrix(moments, hankel_order, 0), options);
    if (!incorporate_matrix(state, moment_matrix, checked_order))
    {
        return state.result;
    }

    if (checked_order % 2 == 1)
    {
        const std::vector<double> next_column = extension_column(moments, hankel_order, 0);
        if (!incorporate_range_consistency(state, moment_matrix, next_column, "H_" + std::to_string(hankel_order) + "_next_column", checked_order,
                                           options))
        {
            return state.result;
        }
    }
    else if (hankel_order > 0 && moment_matrix.status == MatrixStatus::positive_semidefinite)
    {
        const MatrixAnalysis preceding_matrix =
            analyze_matrix("H_" + std::to_string(hankel_order - 1), hankel_order - 1, hankel_matrix(moments, hankel_order - 1, 0), options);
        if (!incorporate_matrix(state, preceding_matrix, checked_order))
        {
            return state.result;
        }
        if (numerical_rank(moment_matrix) != numerical_rank(preceding_matrix))
        {
            return make_result(RealizabilityStatus::not_realizable, checked_order, "H_" + std::to_string(hankel_order) + "_flat_extension",
                               hankel_order, moment_matrix.minimum_eigenvalue, moment_matrix.tolerance,
                               "Moment sequence is not realizable: singular Hamburger moment matrix is not a flat extension.");
        }
    }

    return finish_classification(std::move(state), checked_order, MomentSupport::real_line);
}

[[nodiscard]] RealizabilityResult check_stieltjes(const MomentSet& moments, const RealizabilityOptions& options)
{
    // Stieltjes realizability adds nonnegative-support constraints through
    // shifted Hankel matrices; the tested matrix family depends on parity.
    const std::size_t checked_order = moments.size() - 1;
    ClassificationState state;
    state.result.checked_order = checked_order;

    if (checked_order == 0)
    {
        const DenseMatrix h0 = hankel_matrix(moments, 0, 0);
        const MatrixAnalysis h_analysis = analyze_matrix("H_0", 0, h0, options);
        if (!incorporate_matrix(state, h_analysis, checked_order))
        {
            return state.result;
        }
        return finish_classification(std::move(state), checked_order, MomentSupport::nonnegative);
    }

    const std::size_t order = checked_order / 2;
    MatrixAnalysis moment_matrix = analyze_matrix("H_" + std::to_string(order), order, hankel_matrix(moments, order, 0), options);
    if (!incorporate_matrix(state, moment_matrix, checked_order))
    {
        return state.result;
    }

    if (checked_order % 2 == 1)
    {
        MatrixAnalysis shifted_matrix = analyze_matrix("S_" + std::to_string(order), order, hankel_matrix(moments, order, 1), options);
        if (!incorporate_matrix(state, shifted_matrix, checked_order))
        {
            return state.result;
        }
        const std::vector<double> next_column = extension_column(moments, order, 0);
        if (!incorporate_range_consistency(state, moment_matrix, next_column, "H_" + std::to_string(order) + "_next_column", checked_order, options))
        {
            return state.result;
        }
    }
    else if (order > 0)
    {
        MatrixAnalysis shifted_matrix = analyze_matrix("S_" + std::to_string(order - 1), order - 1, hankel_matrix(moments, order - 1, 1), options);
        if (!incorporate_matrix(state, shifted_matrix, checked_order))
        {
            return state.result;
        }
        const std::vector<double> next_shifted_column = extension_column(moments, order - 1, 1);
        if (!incorporate_range_consistency(state, shifted_matrix, next_shifted_column, "S_" + std::to_string(order - 1) + "_next_column",
                                           checked_order, options))
        {
            return state.result;
        }
    }

    return finish_classification(std::move(state), checked_order, MomentSupport::nonnegative);
}

[[nodiscard]] RealizabilityResult check_hausdorff(const MomentSet& moments, const BoundedInterval interval, const RealizabilityOptions& options)
{
    // Hausdorff realizability on [lower, upper] alternates between endpoint
    // localizing matrices for odd order and moment/interval matrices for even.
    const std::size_t checked_order = moments.size() - 1;
    ClassificationState state;
    state.result.checked_order = checked_order;

    if (checked_order == 0)
    {
        const MatrixAnalysis h_analysis = analyze_matrix("H_0", 0, hankel_matrix(moments, 0, 0), options);
        if (!incorporate_matrix(state, h_analysis, checked_order))
        {
            return state.result;
        }
        return finish_classification(std::move(state), checked_order, MomentSupport::bounded_interval);
    }

    const std::size_t order = checked_order / 2;
    if (checked_order % 2 == 1)
    {
        MatrixAnalysis lower_matrix =
            analyze_matrix("L_lower_" + std::to_string(order), order, interval_linear_localizing_matrix(moments, interval, order, true), options);
        if (!incorporate_matrix(state, lower_matrix, checked_order))
        {
            return state.result;
        }
        MatrixAnalysis upper_matrix =
            analyze_matrix("L_upper_" + std::to_string(order), order, interval_linear_localizing_matrix(moments, interval, order, false), options);
        if (!incorporate_matrix(state, upper_matrix, checked_order))
        {
            return state.result;
        }
    }
    else
    {
        MatrixAnalysis moment_matrix = analyze_matrix("H_" + std::to_string(order), order, hankel_matrix(moments, order, 0), options);
        if (!incorporate_matrix(state, moment_matrix, checked_order))
        {
            return state.result;
        }
        if (order > 0)
        {
            MatrixAnalysis interval_matrix = analyze_matrix("L_interval_" + std::to_string(order - 1), order - 1,
                                                            interval_quadratic_localizing_matrix(moments, interval, order - 1), options);
            if (!incorporate_matrix(state, interval_matrix, checked_order))
            {
                return state.result;
            }
        }
    }

    return finish_classification(std::move(state), checked_order, MomentSupport::bounded_interval);
}

[[nodiscard]] RealizabilityResult with_support_context(RealizabilityResult result, const MomentSupport support)
{
    if (!result.message.empty() && result.message.find("support") == std::string::npos && result.message.find("Support") == std::string::npos)
    {
        result.message += std::string(" Support: ") + support_name(support) + ".";
    }
    return result;
}

} // namespace

RealizabilityResult check_realizability(const MomentSet& moments, const MomentSupport support, const RealizabilityOptions& options)
{
    RealizabilityResult validation = validate_common_inputs(moments, options);
    if (validation.status != RealizabilityStatus::strictly_realizable)
    {
        return validation;
    }

    switch (support)
    {
    case MomentSupport::real_line:
        return with_support_context(check_hamburger(moments, options), support);
    case MomentSupport::nonnegative:
        return with_support_context(check_stieltjes(moments, options), support);
    case MomentSupport::bounded_interval:
        return make_result(RealizabilityStatus::invalid_input, moments.empty() ? 0 : moments.size() - 1, "bounded_interval", 0, 0.0,
                           options.tolerance, "Bounded-interval realizability requires explicit interval bounds; use the BoundedInterval overload.");
    }

    return make_result(RealizabilityStatus::invalid_input, moments.empty() ? 0 : moments.size() - 1, "", 0, 0.0, options.tolerance,
                       "Unknown moment support.");
}

RealizabilityResult check_realizability(const MomentSet& moments, const BoundedInterval interval, const RealizabilityOptions& options)
{
    RealizabilityResult validation = validate_common_inputs(moments, options);
    if (!std::isfinite(interval.lower) || !std::isfinite(interval.upper) || !(interval.lower < interval.upper))
    {
        return make_result(RealizabilityStatus::invalid_input, moments.empty() ? 0 : moments.size() - 1, "bounded_interval", 0, 0.0,
                           options.tolerance, "Bounded interval must have finite bounds with lower < upper.");
    }
    if (validation.status != RealizabilityStatus::strictly_realizable)
    {
        return with_support_context(std::move(validation), MomentSupport::bounded_interval);
    }

    return with_support_context(check_hausdorff(moments, interval, options), MomentSupport::bounded_interval);
}

} // namespace qmomkit
