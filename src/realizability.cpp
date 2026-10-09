#include "qmomkit/realizability.hpp"

#include "internal/linear_algebra.hpp"
#include "internal/numeric_utils.hpp"
#include "internal/realizability.hpp"
#include "qmomkit/realizability_policy.hpp"

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
    // The threshold an eigenvalue must exceed in magnitude to count as
    // positive or negative: the user tolerance plus the data uncertainty.
    double tolerance = 0.0;
    // Eigenvalues within this band of zero count as zero: the roundoff of the
    // eigen solve plus the data uncertainty.
    double zero_band = 0.0;
    // Bound on how far the data uncertainty can move an eigenvalue.
    double uncertainty = 0.0;
    // Eigenvalues within this band are zero for the data as given.
    double roundoff = 0.0;
    EigenDecomposition eigensystem{std::vector<double>{}, DenseMatrix(0, 0)};
};

// A moment-like sequence with a bound on the absolute uncertainty of each
// entry: the standardized moments carry the propagated rounding of the raw
// input moments. Uncertainty that reaches the information scale leaves a
// matrix of the sequence without a single reliable digit: for the standardized
// moments, it is their mass; for a localizing sequence, that mass times the size of the
// coefficients of the localizing polynomial.
struct MomentSequence
{
    std::vector<double> values;
    std::vector<double> uncertainties;
    double information_scale = std::numeric_limits<double>::infinity();
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

// The roundoff allowance of the eigen solve, relative to the scale of the
// matrix. Nothing within it is resolved, whatever tolerance is requested.
constexpr double roundoff_allowance = 1000.0 * std::numeric_limits<double>::epsilon();

enum class EigenvalueClass
{
    negative,
    zero,
    ambiguous,
    positive
};

// One partition serves the definiteness decision, the rank and the kernel:
// an eigenvalue within the zero band is zero, one beyond minus the tolerance
// is negative, one at or above the tolerance is positive, and any other is
// ambiguous. The zero band is tested first because a tolerance below its
// roundoff allowance makes the tolerance smaller than it, and an eigenvalue
// within the roundoff of the eigen solve has no sign.
[[nodiscard]] EigenvalueClass classify_eigenvalue(const MatrixAnalysis& analysis, const double eigenvalue)
{
    if (std::abs(eigenvalue) <= analysis.zero_band)
    {
        return EigenvalueClass::zero;
    }
    if (eigenvalue < -analysis.tolerance)
    {
        return EigenvalueClass::negative;
    }
    if (eigenvalue < analysis.tolerance)
    {
        return EigenvalueClass::ambiguous;
    }
    return EigenvalueClass::positive;
}

[[nodiscard]] DenseMatrix hankel_matrix(const MomentSequence& sequence, const std::size_t order)
{
    // Hankel moment matrices encode nonnegativity of squared polynomials on
    // the support; Hankel matrices of localizing sequences enforce the
    // support bounds.
    DenseMatrix matrix(order + 1, order + 1);
    for (std::size_t row = 0; row <= order; ++row)
    {
        for (std::size_t col = 0; col <= order; ++col)
        {
            matrix(row, col) = sequence.values[row + col];
        }
    }
    return matrix;
}

// The Frobenius norm of the entry uncertainties, which bounds how far they
// can move any eigenvalue of the Hankel matrix.
[[nodiscard]] double hankel_uncertainty(const MomentSequence& sequence, const std::size_t order)
{
    double sum = 0.0;
    for (std::size_t row = 0; row <= order; ++row)
    {
        for (std::size_t col = 0; col <= order; ++col)
        {
            sum += sequence.uncertainties[row + col] * sequence.uncertainties[row + col];
        }
    }
    return std::sqrt(sum);
}

[[nodiscard]] MatrixAnalysis analyze_matrix(std::string name, const std::size_t order, const MomentSequence& sequence,
                                            const RealizabilityOptions& options)
{
    MatrixAnalysis analysis{std::move(name), order};
    const DenseMatrix matrix = hankel_matrix(sequence, order);
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

    const double user_tolerance = options.tolerance * normalization_scale * normalized_infinity_norm;
    const double diagonalization_roundoff = roundoff_allowance * normalization_scale * normalized_infinity_norm;
    analysis.uncertainty = hankel_uncertainty(sequence, order);
    analysis.roundoff = diagonalization_roundoff;
    analysis.tolerance = user_tolerance + analysis.uncertainty;
    analysis.zero_band = diagonalization_roundoff + analysis.uncertainty;
    analysis.eigensystem = internal::symmetric_jacobi_eigen(matrix, {});
    analysis.minimum_eigenvalue = analysis.eigensystem.eigenvalues.empty() ? 0.0 : analysis.eigensystem.eigenvalues.front();

    if (!analysis.eigensystem.converged || !internal::all_finite(analysis.eigensystem.eigenvalues) || !std::isfinite(analysis.uncertainty) ||
        analysis.uncertainty >= sequence.information_scale)
    {
        analysis.status = MatrixStatus::numerically_indeterminate;
        return analysis;
    }

    // Any negative eigenvalue decides; otherwise any ambiguous one leaves the
    // matrix indeterminate; otherwise any zero one makes it singular.
    bool has_zero = false;
    bool has_ambiguous = false;
    for (const double eigenvalue : analysis.eigensystem.eigenvalues)
    {
        switch (classify_eigenvalue(analysis, eigenvalue))
        {
        case EigenvalueClass::negative:
            analysis.status = MatrixStatus::not_positive_semidefinite;
            return analysis;
        case EigenvalueClass::zero:
            has_zero = true;
            break;
        case EigenvalueClass::ambiguous:
            has_ambiguous = true;
            break;
        case EigenvalueClass::positive:
            break;
        }
    }
    if (has_ambiguous)
    {
        analysis.status = MatrixStatus::numerically_indeterminate;
    }
    else if (has_zero)
    {
        analysis.status = MatrixStatus::positive_semidefinite;
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
        state.result = make_result(RealizabilityStatus::numerically_indeterminate, checked_order, analysis.name, analysis.order,
                                   analysis.minimum_eigenvalue, analysis.tolerance,
                                   "Realizability is numerically indeterminate: " + analysis.name +
                                       " has an eigenvalue within the numerical tolerance, or the input moments carry too few digits to decide.");
        return false;
    }

    if (analysis.status == MatrixStatus::positive_semidefinite)
    {
        state.boundary = true;
    }
    return true;
}

// Localizing sequence of (x - lower): entry k is m_(k+1) - lower m_k.
[[nodiscard]] MomentSequence lower_localizing_sequence(const MomentSequence& moments, const double lower)
{
    MomentSequence sequence;
    sequence.information_scale = (1.0 + std::abs(lower)) * moments.information_scale;
    for (std::size_t index = 0; index + 1 < moments.values.size(); ++index)
    {
        sequence.values.push_back(moments.values[index + 1] - lower * moments.values[index]);
        sequence.uncertainties.push_back(moments.uncertainties[index + 1] + std::abs(lower) * moments.uncertainties[index]);
    }
    return sequence;
}

// Localizing sequence of (upper - x): entry k is upper m_k - m_(k+1).
[[nodiscard]] MomentSequence upper_localizing_sequence(const MomentSequence& moments, const double upper)
{
    MomentSequence sequence;
    sequence.information_scale = (1.0 + std::abs(upper)) * moments.information_scale;
    for (std::size_t index = 0; index + 1 < moments.values.size(); ++index)
    {
        sequence.values.push_back(upper * moments.values[index] - moments.values[index + 1]);
        sequence.uncertainties.push_back(std::abs(upper) * moments.uncertainties[index] + moments.uncertainties[index + 1]);
    }
    return sequence;
}

// Localizing sequence of (x - lower)(upper - x), enforcing support inside the
// bounded interval for even-order Hausdorff checks: entry k is
// (lower + upper) m_(k+1) - lower upper m_k - m_(k+2).
[[nodiscard]] MomentSequence interval_localizing_sequence(const MomentSequence& moments, const double lower, const double upper)
{
    MomentSequence sequence;
    sequence.information_scale = (1.0 + std::abs(lower + upper) + std::abs(lower * upper)) * moments.information_scale;
    for (std::size_t index = 0; index + 2 < moments.values.size(); ++index)
    {
        sequence.values.push_back((lower + upper) * moments.values[index + 1] - lower * upper * moments.values[index] - moments.values[index + 2]);
        sequence.uncertainties.push_back(std::abs(lower + upper) * moments.uncertainties[index + 1] +
                                         std::abs(lower * upper) * moments.uncertainties[index] + moments.uncertainties[index + 2]);
    }
    return sequence;
}

// The column that extends the Hankel matrix of the given order by one.
[[nodiscard]] MomentSequence extension_column(const MomentSequence& sequence, const std::size_t matrix_order)
{
    MomentSequence column;
    for (std::size_t row = 0; row <= matrix_order; ++row)
    {
        column.values.push_back(sequence.values[row + matrix_order + 1]);
        column.uncertainties.push_back(sequence.uncertainties[row + matrix_order + 1]);
    }
    return column;
}

[[nodiscard]] double euclidean_norm(const std::vector<double>& values)
{
    double sum = 0.0;
    for (const double value : values)
    {
        sum += value * value;
    }
    return std::sqrt(sum);
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

// The threshold for a kernel projection of the extension column. Beyond the
// tolerance, it covers the data uncertainty: that of the column itself, and
// the rotation of the kernel vectors, which the data uncertainty can turn by
// up to its size over the gap between the zero and the positive eigenvalues
// (Davis-Kahan); zero eigenvalues lie within the zero band, so the smallest
// positive eigenvalue minus the band bounds that gap from below.
[[nodiscard]] double kernel_projection_tolerance(const MatrixAnalysis& analysis, const MomentSequence& column, const RealizabilityOptions& options)
{
    const double tolerance = std::max(analysis.tolerance, options.tolerance * vector_scale(column.values));
    double smallest_positive = std::numeric_limits<double>::infinity();
    for (const double eigenvalue : analysis.eigensystem.eigenvalues)
    {
        if (classify_eigenvalue(analysis, eigenvalue) == EigenvalueClass::positive)
        {
            smallest_positive = std::min(smallest_positive, eigenvalue);
        }
    }
    const double kernel_rotation = std::isfinite(smallest_positive) ? analysis.uncertainty / (smallest_positive - analysis.zero_band) : 0.0;
    const double uncertainty = euclidean_norm(column.uncertainties) + euclidean_norm(column.values) * kernel_rotation;
    // The factor ten keeps the uncertainty below the tenth of the tolerance
    // under which a projection counts as consistent.
    return std::max(tolerance, 10.0 * uncertainty);
}

[[nodiscard]] RealizabilityStatus check_range_consistency(const MatrixAnalysis& analysis, const MomentSequence& column, const double tolerance,
                                                          double& maximum_kernel_projection)
{
    // For a singular moment/localizing matrix, the next extension column must
    // lie in its range. Projection onto nullspace eigenvectors detects failure.
    maximum_kernel_projection = 0.0;
    for (std::size_t col = 0; col < analysis.eigensystem.eigenvalues.size(); ++col)
    {
        if (classify_eigenvalue(analysis, analysis.eigensystem.eigenvalues[col]) == EigenvalueClass::zero)
        {
            const double projection = std::abs(projection_on_eigenvector(analysis.eigensystem.eigenvectors, col, column.values));
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

[[nodiscard]] bool incorporate_range_consistency(ClassificationState& state, const MatrixAnalysis& analysis, const MomentSequence& column,
                                                 const std::string& name, const std::size_t checked_order, const RealizabilityOptions& options)
{
    if (analysis.status == MatrixStatus::positive_definite)
    {
        return true;
    }

    const double tolerance = kernel_projection_tolerance(analysis, column, options);
    double maximum_kernel_projection = 0.0;
    const RealizabilityStatus status = check_range_consistency(analysis, column, tolerance, maximum_kernel_projection);
    state.boundary = true;
    if (status == RealizabilityStatus::boundary)
    {
        return true;
    }

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
                                                  [&analysis](const double eigenvalue)
                                                  { return classify_eigenvalue(analysis, eigenvalue) == EigenvalueClass::positive; }));
}

// Whether every zero eigenvalue of the matrix is zero for the data as given,
// within the roundoff of the eigen solve, rather than only within the data
// uncertainty.
[[nodiscard]] bool singular_as_given(const MatrixAnalysis& analysis)
{
    return std::all_of(analysis.eigensystem.eigenvalues.begin(), analysis.eigensystem.eigenvalues.end(), [&analysis](const double eigenvalue)
                       { return classify_eigenvalue(analysis, eigenvalue) != EigenvalueClass::zero || std::abs(eigenvalue) <= analysis.roundoff; });
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

[[nodiscard]] const char* tolerance_error(const RealizabilityOptions& options)
{
    return std::isfinite(options.tolerance) && options.tolerance > 0.0 ? nullptr : "Realizability tolerance must be finite and positive.";
}

[[nodiscard]] const char* interval_error(const BoundedInterval interval)
{
    return std::isfinite(interval.lower) && std::isfinite(interval.upper) && interval.lower < interval.upper
               ? nullptr
               : "Bounded interval must have finite bounds with lower < upper.";
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
    if (const char* error = tolerance_error(options); error != nullptr)
    {
        return make_result(RealizabilityStatus::invalid_input, moments.size() - 1, "", 0, 0.0, 0.0, error);
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

[[nodiscard]] RealizabilityResult check_hamburger(const MomentSequence& moments, const RealizabilityOptions& options)
{
    // Hamburger realizability on the real line is governed by PSD Hankel moment
    // matrices. A singular even truncation must be a flat extension of its
    // preceding moment matrix; odd truncations require next-column range validity.
    const std::size_t checked_order = moments.values.size() - 1;
    const std::size_t hankel_order = checked_order / 2;
    ClassificationState state;
    state.result.checked_order = checked_order;

    MatrixAnalysis moment_matrix = analyze_matrix("H_" + std::to_string(hankel_order), hankel_order, moments, options);
    if (!incorporate_matrix(state, moment_matrix, checked_order))
    {
        return state.result;
    }

    if (checked_order % 2 == 1)
    {
        if (!incorporate_range_consistency(state, moment_matrix, extension_column(moments, hankel_order),
                                           "H_" + std::to_string(hankel_order) + "_next_column", checked_order, options))
        {
            return state.result;
        }
    }
    else if (hankel_order > 0 && moment_matrix.status == MatrixStatus::positive_semidefinite)
    {
        const MatrixAnalysis preceding_matrix = analyze_matrix("H_" + std::to_string(hankel_order - 1), hankel_order - 1, moments, options);
        if (!incorporate_matrix(state, preceding_matrix, checked_order))
        {
            return state.result;
        }
        // H_n is a flat extension if it has no more nonzero eigenvalues than
        // H_(n-1). The zero band of H_n is the wider one, since H_(n-1) is its
        // principal submatrix, so H_n can count fewer nonzero eigenvalues than
        // H_(n-1); it then only resolves them less finely, and the extension
        // is flat within its precision. More nonzero eigenvalues decide that
        // the extension is not flat only if the zeros of H_(n-1) hold for the
        // data as given. If they rest on the data uncertainty, H_(n-1) may as
        // well be positive definite, and the sequence realizable whatever H_n.
        if (numerical_rank(moment_matrix) > numerical_rank(preceding_matrix))
        {
            const bool decided = singular_as_given(preceding_matrix);
            return make_result(decided ? RealizabilityStatus::not_realizable : RealizabilityStatus::numerically_indeterminate, checked_order,
                               "H_" + std::to_string(hankel_order) + "_flat_extension", hankel_order, moment_matrix.minimum_eigenvalue,
                               moment_matrix.tolerance,
                               decided ? "Moment sequence is not realizable: singular Hamburger moment matrix is not a flat extension."
                                       : "Realizability is numerically indeterminate: the flat extension of the Hamburger moment matrix is "
                                         "decided only within the data uncertainty.");
        }
    }

    return finish_classification(std::move(state), checked_order, MomentSupport::real_line);
}

[[nodiscard]] RealizabilityResult check_stieltjes(const MomentSequence& moments, const double lower, const RealizabilityOptions& options)
{
    // Stieltjes realizability on [lower, infinity) adds the localizing matrix
    // of (x - lower); the tested matrix family depends on parity.
    const std::size_t checked_order = moments.values.size() - 1;
    ClassificationState state;
    state.result.checked_order = checked_order;

    if (checked_order == 0)
    {
        const MatrixAnalysis h_analysis = analyze_matrix("H_0", 0, moments, options);
        if (!incorporate_matrix(state, h_analysis, checked_order))
        {
            return state.result;
        }
        return finish_classification(std::move(state), checked_order, MomentSupport::nonnegative);
    }

    const std::size_t order = checked_order / 2;
    MatrixAnalysis moment_matrix = analyze_matrix("H_" + std::to_string(order), order, moments, options);
    if (!incorporate_matrix(state, moment_matrix, checked_order))
    {
        return state.result;
    }

    const MomentSequence localized = lower_localizing_sequence(moments, lower);
    if (checked_order % 2 == 1)
    {
        MatrixAnalysis shifted_matrix = analyze_matrix("S_" + std::to_string(order), order, localized, options);
        if (!incorporate_matrix(state, shifted_matrix, checked_order))
        {
            return state.result;
        }
        if (!incorporate_range_consistency(state, moment_matrix, extension_column(moments, order), "H_" + std::to_string(order) + "_next_column",
                                           checked_order, options))
        {
            return state.result;
        }
    }
    else if (order > 0)
    {
        MatrixAnalysis shifted_matrix = analyze_matrix("S_" + std::to_string(order - 1), order - 1, localized, options);
        if (!incorporate_matrix(state, shifted_matrix, checked_order))
        {
            return state.result;
        }
        if (!incorporate_range_consistency(state, shifted_matrix, extension_column(localized, order - 1),
                                           "S_" + std::to_string(order - 1) + "_next_column", checked_order, options))
        {
            return state.result;
        }
    }

    return finish_classification(std::move(state), checked_order, MomentSupport::nonnegative);
}

[[nodiscard]] RealizabilityResult check_hausdorff(const MomentSequence& moments, const double lower, const double upper,
                                                  const RealizabilityOptions& options)
{
    // Hausdorff realizability on [lower, upper] alternates between endpoint
    // localizing matrices for odd order and moment/interval matrices for even.
    const std::size_t checked_order = moments.values.size() - 1;
    ClassificationState state;
    state.result.checked_order = checked_order;

    if (checked_order == 0)
    {
        const MatrixAnalysis h_analysis = analyze_matrix("H_0", 0, moments, options);
        if (!incorporate_matrix(state, h_analysis, checked_order))
        {
            return state.result;
        }
        return finish_classification(std::move(state), checked_order, MomentSupport::bounded_interval);
    }

    const std::size_t order = checked_order / 2;
    if (checked_order % 2 == 1)
    {
        MatrixAnalysis lower_matrix = analyze_matrix("L_lower_" + std::to_string(order), order, lower_localizing_sequence(moments, lower), options);
        if (!incorporate_matrix(state, lower_matrix, checked_order))
        {
            return state.result;
        }
        MatrixAnalysis upper_matrix = analyze_matrix("L_upper_" + std::to_string(order), order, upper_localizing_sequence(moments, upper), options);
        if (!incorporate_matrix(state, upper_matrix, checked_order))
        {
            return state.result;
        }
    }
    else
    {
        MatrixAnalysis moment_matrix = analyze_matrix("H_" + std::to_string(order), order, moments, options);
        if (!incorporate_matrix(state, moment_matrix, checked_order))
        {
            return state.result;
        }
        if (order > 0)
        {
            MatrixAnalysis interval_matrix =
                analyze_matrix("L_interval_" + std::to_string(order - 1), order - 1, interval_localizing_sequence(moments, lower, upper), options);
            if (!incorporate_matrix(state, interval_matrix, checked_order))
            {
                return state.result;
            }
        }
    }

    return finish_classification(std::move(state), checked_order, MomentSupport::bounded_interval);
}

// The moments the check tests: those of y = (x - c) / s,
// with c = m1 / m0 and s the power of two just above the spread, divided by
// the power of two just above m0. Realizability is invariant under this
// affine change of coordinate, and the moment matrices change by a
// congruence, which keeps the signs of their eigenvalues; the tolerance,
// relative to the standardized matrices, then no longer depends on where the
// distribution sits or in which unit it is given.
struct TestedMoments
{
    MomentSequence moments;
    double centre = 0.0;
    int scale_exponent = 0;

    // A support bound in the tested coordinate. Its rounding, like that of
    // building the localizing sequences, is relative to the terms involved
    // and lies within the roundoff allowance of the eigen solve.
    [[nodiscard]] double bound(const double value) const
    {
        return std::ldexp(value - centre, -scale_exponent);
    }
};

[[nodiscard]] TestedMoments standardized_moments(const MomentSet& moments)
{
    const std::size_t count = moments.size();
    const double mass = moments[0];
    TestedMoments tested;
    tested.centre = count > 1 ? moments[1] / mass : 0.0;

    // Dividing by the power of two of m0 first is exact and keeps the centring
    // sums finite for moments near the overflow limit.
    int mass_exponent = 0;
    static_cast<void>(std::frexp(mass, &mass_exponent));
    std::vector<double> scaled(count);
    for (std::size_t k = 0; k < count; ++k)
    {
        scaled[k] = std::ldexp(moments[k], -mass_exponent);
    }

    // mu_k = sum over j of C(k, j) (-c)^(k-j) m_j. Each term passes through
    // at most 2k + 1 roundings here, each a factor 1 + delta with |delta| <=
    // epsilon / 2, which costs less than (k + 1) epsilon of the sum of the
    // absolute terms; the input moments, taken as exact to their last digit,
    // add one epsilon. Where the terms cancel, this limits the resolution.
    std::vector<double> central(count, 0.0);
    std::vector<double> central_uncertainty(count, 0.0);
    std::vector<double> binomial_row(count, 0.0);
    std::vector<double> centre_powers(count, 1.0);
    binomial_row[0] = 1.0;
    for (std::size_t k = 0; k < count; ++k)
    {
        if (k > 0)
        {
            for (std::size_t j = k; j > 0; --j)
            {
                binomial_row[j] += binomial_row[j - 1];
            }
            centre_powers[k] = centre_powers[k - 1] * -tested.centre;
        }
        double magnitude = 0.0;
        for (std::size_t j = 0; j <= k; ++j)
        {
            const double term = binomial_row[j] * centre_powers[k - j] * scaled[j];
            central[k] += term;
            magnitude += std::abs(term);
        }
        central_uncertainty[k] = static_cast<double>(k + 2) * std::numeric_limits<double>::epsilon() * magnitude;
    }

    // Scale by the spread, taken from the first even central moment that is
    // resolved, of either sign; without one, by the distance from the origin,
    // so that a point mass is not stretched.
    double spread = tested.centre != 0.0 ? std::abs(tested.centre) : 1.0;
    for (std::size_t k = 2; k < count; k += 2)
    {
        if (std::abs(central[k]) > 2.0 * central_uncertainty[k])
        {
            spread = std::pow(std::abs(central[k]) / scaled[0], 1.0 / static_cast<double>(k));
            break;
        }
    }
    static_cast<void>(std::frexp(spread, &tested.scale_exponent));

    tested.moments.values.resize(count);
    tested.moments.uncertainties.resize(count);
    for (std::size_t k = 0; k < count; ++k)
    {
        const long long shift = -static_cast<long long>(tested.scale_exponent) * static_cast<long long>(k);
        const int clamped = static_cast<int>(std::clamp(shift, -100000LL, 100000LL));
        tested.moments.values[k] = std::ldexp(central[k], clamped);
        tested.moments.uncertainties[k] = std::ldexp(central_uncertainty[k], clamped);
    }
    tested.moments.information_scale = tested.moments.values[0];
    return tested;
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
    {
        const TestedMoments tested = standardized_moments(moments);
        return with_support_context(check_hamburger(tested.moments, options), support);
    }
    case MomentSupport::nonnegative:
    {
        const TestedMoments tested = standardized_moments(moments);
        return with_support_context(check_stieltjes(tested.moments, tested.bound(0.0), options), support);
    }
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
    if (const char* error = interval_error(interval); error != nullptr)
    {
        return make_result(RealizabilityStatus::invalid_input, moments.empty() ? 0 : moments.size() - 1, "bounded_interval", 0, 0.0,
                           options.tolerance, error);
    }
    if (validation.status != RealizabilityStatus::strictly_realizable)
    {
        return with_support_context(std::move(validation), MomentSupport::bounded_interval);
    }

    const TestedMoments tested = standardized_moments(moments);
    return with_support_context(check_hausdorff(tested.moments, tested.bound(interval.lower), tested.bound(interval.upper), options),
                                MomentSupport::bounded_interval);
}

const char* internal::realizability_policy_error(const RealizabilityPolicy& policy)
{
    if (policy.check_mode == RealizabilityCheckMode::None)
    {
        return nullptr;
    }
    if (policy.check_mode != RealizabilityCheckMode::EveryStage && policy.check_mode != RealizabilityCheckMode::FinalOnly)
    {
        return "Unknown realizability check mode.";
    }
    if (const char* error = tolerance_error(policy.options); error != nullptr)
    {
        return error;
    }
    switch (policy.support)
    {
    case MomentSupport::real_line:
    case MomentSupport::nonnegative:
        return nullptr;
    case MomentSupport::bounded_interval:
        return policy.interval.has_value() ? interval_error(*policy.interval) : "A bounded-interval realizability policy requires an interval.";
    }
    return "Unknown moment support.";
}

} // namespace qmomkit
