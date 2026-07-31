#include "example_workflow.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/wheeler.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

// The diagnostic demonstrates the physical aspect-ratio support that the
// unbounded configuration intentionally omits.
constexpr std::size_t chi_dimension = 1;
constexpr double chi_lower = 0.0;
constexpr double chi_upper = 1.0;

struct Matrix
{
    std::size_t rows = 0;
    std::size_t cols = 0;
    std::vector<double> values;

    Matrix() = default;

    Matrix(const std::size_t row_count, const std::size_t col_count) : rows(row_count), cols(col_count), values(row_count * col_count, 0.0) {}

    double& operator()(const std::size_t row, const std::size_t col)
    {
        return values.at(row * cols + col);
    }

    double operator()(const std::size_t row, const std::size_t col) const
    {
        return values.at(row * cols + col);
    }
};

struct EigenResult
{
    std::vector<double> eigenvalues;
    Matrix eigenvectors;
};

struct MatrixDiagnostic
{
    std::string family;
    std::size_t k = 0;
    std::size_t l = 0;
    Matrix matrix;
    double determinant = 0.0;
    std::vector<double> eigenvalues;
    double minimum_eigenvalue = 0.0;
    double infinity_norm = 0.0;
    double effective_tolerance = 0.0;
    double signed_margin = 0.0;
    std::string classification;
};

struct ConditionalSolveDiagnostic
{
    std::vector<double> parent_nodes;
    std::vector<double> parent_weights;
    std::vector<double> parent_global_weights;
    Matrix vandermonde;
    Matrix rhs;
    Matrix solved;
    double condition_estimate = 0.0;
    double residual_norm = 0.0;
    double maximum_input_scale = 0.0;
    double maximum_output_scale = 0.0;
};

struct StageMoments
{
    std::string label;
    double time = 0.0;
    std::vector<double> values;
    std::vector<double> branch_moments;
};

struct SourceReplay
{
    std::optional<qmomkit::CqmomMomentLayout> layout;
    std::vector<StageMoments> stages;
    double max_difference_from_public_source_update = 0.0;
};

struct OffendingNode
{
    std::size_t terminal_node = 0;
    qmomkit::ConditioningPath path;
    std::vector<double> state;
    double terminal_weight = 0.0;
    double parent_global_weight = 0.0;
    double chi_conditional_weight = 0.0;
    std::vector<double> chi_nodes;
    std::vector<double> chi_conditional_weights;
    std::vector<double> chi_global_weights;
};

struct BoundedAcceptedBranch
{
    qmomkit::ConditioningPath parent_path;
    std::vector<double> chi_nodes;
    std::vector<double> conditional_weights;
    std::vector<double> global_weights;
    std::vector<double> canonical_moments;
    std::vector<double> retained_prefix_absolute_residuals;
    std::vector<double> retained_prefix_relative_residuals;
    double parent_global_weight = 0.0;
};

struct BoundaryAdvance
{
    bool ok = false;
    std::string message;
    std::optional<qmomkit::CqmomState> pre_source_state;
    std::optional<qmomkit::CqmomQuadrature> pre_source_quadrature;
    std::optional<qmomkit::CqmomState> boundary_state;
    std::optional<qmomkit::CqmomQuadrature> boundary_quadrature;
    std::optional<qmomkit::CqmomState> final_state;
    std::optional<qmomkit::CqmomQuadrature> final_quadrature;
    qmomkit::TransportDiagnostics source_diagnostics;
    qmomkit::TransportStage failed_stage = qmomkit::TransportStage::none;
    double source_start_time = 0.0;
    double boundary_time = 0.0;
    std::string stage_label;
};

struct DiagnosticReport
{
    std::size_t macro_step = 0;
    double time = 0.0;
    std::string operator_stage_label;
    OffendingNode offending;
    qmomkit::CqmomDiagnostics unbounded_branch;
    std::optional<qmomkit::CqmomDiagnostics> bounded_branch;
    std::optional<qmomkit::CqmomNodeCountTransition> bounded_transition;
    std::optional<BoundedAcceptedBranch> bounded_accepted_branch;
    ConditionalSolveDiagnostic conditional_solve;
    SourceReplay source_replay;
    std::vector<double> branch_moments;
    std::vector<double> reconstructed_moments;
    qmomkit::Quadrature1D wheeler;
    double configured_rmin = 0.0;
    std::vector<double> recurrence_a;
    std::vector<double> recurrence_b;
    Matrix sigma;
    std::vector<MatrixDiagnostic> matrices;
};

std::string path_string(const qmomkit::ConditioningPath& path)
{
    std::ostringstream text;
    text << '[';
    for (std::size_t level = 0; level < path.size(); ++level)
    {
        if (level > 0)
        {
            text << ',';
        }
        text << path.node_index_at_level(level);
    }
    text << ']';
    return text.str();
}

bool path_has_prefix(const qmomkit::ConditioningPath& path, const std::vector<std::size_t>& prefix)
{
    if (path.size() < prefix.size())
    {
        return false;
    }
    for (std::size_t level = 0; level < prefix.size(); ++level)
    {
        if (path.node_index_at_level(level) != prefix[level])
        {
            return false;
        }
    }
    return true;
}

bool same_path(const qmomkit::ConditioningPath& lhs, const qmomkit::ConditioningPath& rhs)
{
    if (lhs.size() != rhs.size())
    {
        return false;
    }
    for (std::size_t level = 0; level < lhs.size(); ++level)
    {
        if (lhs.node_index_at_level(level) != rhs.node_index_at_level(level))
        {
            return false;
        }
    }
    return true;
}

std::string vector_string(const std::vector<double>& values)
{
    std::ostringstream text;
    text << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (i > 0)
        {
            text << ';';
        }
        text << values[i];
    }
    return text.str();
}

std::string vector_string(const std::vector<std::size_t>& values)
{
    std::ostringstream text;
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (i > 0)
        {
            text << ';';
        }
        text << values[i];
    }
    return text.str();
}

std::string matrix_string(const Matrix& matrix)
{
    std::ostringstream text;
    text << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (std::size_t row = 0; row < matrix.rows; ++row)
    {
        if (row > 0)
        {
            text << " / ";
        }
        text << '[';
        for (std::size_t col = 0; col < matrix.cols; ++col)
        {
            if (col > 0)
            {
                text << "; ";
            }
            text << matrix(row, col);
        }
        text << ']';
    }
    return text.str();
}

std::string transport_scheme_string(const qmomkit::TransportScheme scheme)
{
    switch (scheme)
    {
    case qmomkit::TransportScheme::strang:
        return "strang";
    case qmomkit::TransportScheme::cd:
        return "cd";
    case qmomkit::TransportScheme::dc:
        return "dc";
    }
    return "unknown";
}

std::string transition_kind_string(const qmomkit::CqmomNodeCountTransitionKind kind)
{
    switch (kind)
    {
    case qmomkit::CqmomNodeCountTransitionKind::reduced:
        return "reduced";
    case qmomkit::CqmomNodeCountTransitionKind::grew:
        return "grew";
    case qmomkit::CqmomNodeCountTransitionKind::created:
        return "created";
    case qmomkit::CqmomNodeCountTransitionKind::removed:
        return "removed";
    }
    return "unknown";
}

double infinity_norm(const Matrix& matrix)
{
    double norm = 0.0;
    for (std::size_t row = 0; row < matrix.rows; ++row)
    {
        double row_sum = 0.0;
        for (std::size_t col = 0; col < matrix.cols; ++col)
        {
            row_sum += std::abs(matrix(row, col));
        }
        norm = std::max(norm, row_sum);
    }
    return norm;
}

double determinant(Matrix matrix)
{
    if (matrix.rows != matrix.cols)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double value = 1.0;
    for (std::size_t pivot_col = 0; pivot_col < matrix.rows; ++pivot_col)
    {
        std::size_t pivot_row = pivot_col;
        double pivot_abs = std::abs(matrix(pivot_col, pivot_col));
        for (std::size_t row = pivot_col + 1; row < matrix.rows; ++row)
        {
            const double candidate = std::abs(matrix(row, pivot_col));
            if (candidate > pivot_abs)
            {
                pivot_abs = candidate;
                pivot_row = row;
            }
        }
        if (pivot_abs == 0.0)
        {
            return 0.0;
        }
        if (pivot_row != pivot_col)
        {
            for (std::size_t col = 0; col < matrix.cols; ++col)
            {
                std::swap(matrix(pivot_col, col), matrix(pivot_row, col));
            }
            value = -value;
        }
        const double pivot = matrix(pivot_col, pivot_col);
        value *= pivot;
        for (std::size_t row = pivot_col + 1; row < matrix.rows; ++row)
        {
            const double factor = matrix(row, pivot_col) / pivot;
            for (std::size_t col = pivot_col + 1; col < matrix.cols; ++col)
            {
                matrix(row, col) -= factor * matrix(pivot_col, col);
            }
        }
    }
    return value;
}

EigenResult symmetric_jacobi_eigen(Matrix matrix)
{
    if (matrix.rows != matrix.cols)
    {
        throw std::invalid_argument("Eigenvalue diagnostics require a square matrix.");
    }
    const std::size_t n = matrix.rows;
    Matrix eigenvectors(n, n);
    for (std::size_t i = 0; i < n; ++i)
    {
        eigenvectors(i, i) = 1.0;
    }
    constexpr std::size_t max_sweeps_per_node = 80;
    constexpr double tolerance = 1.0e-14;
    const std::size_t max_sweeps = max_sweeps_per_node * std::max<std::size_t>(n, 1);
    for (std::size_t sweep = 0; sweep < max_sweeps; ++sweep)
    {
        std::size_t p = 0;
        std::size_t q = 0;
        double largest = 0.0;
        for (std::size_t row = 0; row < n; ++row)
        {
            for (std::size_t col = row + 1; col < n; ++col)
            {
                const double candidate = std::abs(matrix(row, col));
                if (candidate > largest)
                {
                    largest = candidate;
                    p = row;
                    q = col;
                }
            }
        }
        if (largest <= tolerance)
        {
            break;
        }
        const double tau = (matrix(q, q) - matrix(p, p)) / (2.0 * matrix(p, q));
        const double t = (tau >= 0.0 ? 1.0 : -1.0) / (std::abs(tau) + std::sqrt(1.0 + tau * tau));
        const double c = 1.0 / std::sqrt(1.0 + t * t);
        const double s = t * c;
        const double app = matrix(p, p);
        const double aqq = matrix(q, q);
        const double apq = matrix(p, q);
        matrix(p, p) = app - t * apq;
        matrix(q, q) = aqq + t * apq;
        matrix(p, q) = 0.0;
        matrix(q, p) = 0.0;
        for (std::size_t k = 0; k < n; ++k)
        {
            if (k != p && k != q)
            {
                const double akp = matrix(k, p);
                const double akq = matrix(k, q);
                matrix(k, p) = c * akp - s * akq;
                matrix(p, k) = matrix(k, p);
                matrix(k, q) = s * akp + c * akq;
                matrix(q, k) = matrix(k, q);
            }
            const double vkp = eigenvectors(k, p);
            const double vkq = eigenvectors(k, q);
            eigenvectors(k, p) = c * vkp - s * vkq;
            eigenvectors(k, q) = s * vkp + c * vkq;
        }
    }
    std::vector<double> eigenvalues(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        eigenvalues[i] = matrix(i, i);
    }
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](const std::size_t lhs, const std::size_t rhs) { return eigenvalues[lhs] < eigenvalues[rhs]; });
    std::vector<double> sorted_values(n, 0.0);
    Matrix sorted_vectors(n, n);
    for (std::size_t sorted = 0; sorted < n; ++sorted)
    {
        sorted_values[sorted] = eigenvalues[order[sorted]];
        for (std::size_t row = 0; row < n; ++row)
        {
            sorted_vectors(row, sorted) = eigenvectors(row, order[sorted]);
        }
    }
    return EigenResult{std::move(sorted_values), std::move(sorted_vectors)};
}

std::string classify_matrix(const double minimum_eigenvalue, const double effective_tolerance, const double norm)
{
    const double diagonalization_roundoff = 1000.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, norm);
    if (minimum_eigenvalue < -effective_tolerance)
    {
        return "not_positive_semidefinite";
    }
    if (std::abs(minimum_eigenvalue) <= diagonalization_roundoff)
    {
        return "positive_semidefinite";
    }
    if (minimum_eigenvalue < effective_tolerance)
    {
        return "numerically_indeterminate";
    }
    return "positive_definite";
}

MatrixDiagnostic analyze_matrix(std::string family, const std::size_t k, const std::size_t l, Matrix matrix)
{
    MatrixDiagnostic diagnostic;
    diagnostic.family = std::move(family);
    diagnostic.k = k;
    diagnostic.l = l;
    diagnostic.matrix = std::move(matrix);
    diagnostic.determinant = determinant(diagnostic.matrix);
    diagnostic.infinity_norm = infinity_norm(diagnostic.matrix);
    diagnostic.effective_tolerance = 1.0e-12 * std::max(1.0, diagnostic.infinity_norm);
    const EigenResult eigen = symmetric_jacobi_eigen(diagnostic.matrix);
    diagnostic.eigenvalues = eigen.eigenvalues;
    diagnostic.minimum_eigenvalue = diagnostic.eigenvalues.empty() ? 0.0 : diagnostic.eigenvalues.front();
    diagnostic.signed_margin = diagnostic.minimum_eigenvalue;
    diagnostic.classification = classify_matrix(diagnostic.minimum_eigenvalue, diagnostic.effective_tolerance, diagnostic.infinity_norm);
    return diagnostic;
}

Matrix hankel_matrix(const std::vector<double>& moments, const std::size_t k, const std::size_t l)
{
    Matrix matrix(l + 1, l + 1);
    for (std::size_t row = 0; row <= l; ++row)
    {
        for (std::size_t col = 0; col <= l; ++col)
        {
            matrix(row, col) = moments.at(k + row + col);
        }
    }
    return matrix;
}

Matrix localizing_matrix(const std::vector<double>& moments, const std::size_t l, const bool lower)
{
    Matrix matrix(l + 1, l + 1);
    for (std::size_t row = 0; row <= l; ++row)
    {
        for (std::size_t col = 0; col <= l; ++col)
        {
            const std::size_t index = row + col;
            matrix(row, col) = lower ? moments.at(index + 1) - chi_lower * moments.at(index) : chi_upper * moments.at(index) - moments.at(index + 1);
        }
    }
    return matrix;
}

std::vector<MatrixDiagnostic> build_matrix_diagnostics(const std::vector<double>& moments)
{
    std::vector<MatrixDiagnostic> diagnostics;
    for (std::size_t k = 0; k <= 1; ++k)
    {
        for (std::size_t l = 0; k + 2 * l < moments.size(); ++l)
        {
            diagnostics.push_back(analyze_matrix("classical_hankel_hadamard", k, l, hankel_matrix(moments, k, l)));
        }
    }
    for (std::size_t l = 0; 2 * l + 1 < moments.size(); ++l)
    {
        diagnostics.push_back(analyze_matrix("lower_localizing", 0, l, localizing_matrix(moments, l, true)));
        diagnostics.push_back(analyze_matrix("upper_localizing", 0, l, localizing_matrix(moments, l, false)));
    }
    return diagnostics;
}

std::vector<double> solve_linear_system(Matrix a, std::vector<double> rhs)
{
    if (a.rows != a.cols || rhs.size() != a.rows)
    {
        throw std::invalid_argument("Linear solve dimensions are inconsistent.");
    }
    const std::size_t n = a.rows;
    for (std::size_t pivot_col = 0; pivot_col < n; ++pivot_col)
    {
        std::size_t pivot_row = pivot_col;
        double pivot_abs = std::abs(a(pivot_col, pivot_col));
        for (std::size_t row = pivot_col + 1; row < n; ++row)
        {
            const double candidate = std::abs(a(row, pivot_col));
            if (candidate > pivot_abs)
            {
                pivot_abs = candidate;
                pivot_row = row;
            }
        }
        if (pivot_abs == 0.0)
        {
            throw std::runtime_error("Weighted Vandermonde diagnostic solve is singular.");
        }
        if (pivot_row != pivot_col)
        {
            for (std::size_t col = 0; col < n; ++col)
            {
                std::swap(a(pivot_col, col), a(pivot_row, col));
            }
            std::swap(rhs[pivot_col], rhs[pivot_row]);
        }
        for (std::size_t row = pivot_col + 1; row < n; ++row)
        {
            const double factor = a(row, pivot_col) / a(pivot_col, pivot_col);
            a(row, pivot_col) = 0.0;
            for (std::size_t col = pivot_col + 1; col < n; ++col)
            {
                a(row, col) -= factor * a(pivot_col, col);
            }
            rhs[row] -= factor * rhs[pivot_col];
        }
    }
    std::vector<double> solution(n, 0.0);
    for (std::size_t row = n; row-- > 0;)
    {
        double value = rhs[row];
        for (std::size_t col = row + 1; col < n; ++col)
        {
            value -= a(row, col) * solution[col];
        }
        solution[row] = value / a(row, row);
    }
    return solution;
}

double residual_infinity_norm(const Matrix& a, const Matrix& x, const Matrix& b)
{
    double norm = 0.0;
    for (std::size_t row = 0; row < a.rows; ++row)
    {
        for (std::size_t col = 0; col < x.cols; ++col)
        {
            double ax = 0.0;
            for (std::size_t k = 0; k < a.cols; ++k)
            {
                ax += a(row, k) * x(k, col);
            }
            norm = std::max(norm, std::abs(ax - b(row, col)));
        }
    }
    return norm;
}

double inverse_infinity_norm(const Matrix& a)
{
    Matrix inverse(a.rows, a.cols);
    for (std::size_t col = 0; col < a.cols; ++col)
    {
        std::vector<double> rhs(a.rows, 0.0);
        rhs[col] = 1.0;
        const std::vector<double> solution = solve_linear_system(a, rhs);
        for (std::size_t row = 0; row < a.rows; ++row)
        {
            inverse(row, col) = solution[row];
        }
    }
    return infinity_norm(inverse);
}

double max_abs_value(const Matrix& matrix)
{
    double scale = 0.0;
    for (const double value : matrix.values)
    {
        scale = std::max(scale, std::abs(value));
    }
    return scale;
}

std::optional<double> moment_value(const qmomkit::CqmomMomentLayout& layout, const std::vector<double>& values, const qmomkit::MultiIndex& index)
{
    const auto found = std::find(layout.required_indices.begin(), layout.required_indices.end(), index);
    if (found == layout.required_indices.end())
    {
        return std::nullopt;
    }
    const std::size_t position = static_cast<std::size_t>(std::distance(layout.required_indices.begin(), found));
    return values.at(position);
}

ConditionalSolveDiagnostic conditional_chi_solve(const qmomkit::CqmomQuadrature& quadrature, const qmomkit::CqmomMomentLayout& layout,
                                                 const std::vector<double>& values)
{
    // Reproduce the branch-local weighted Vandermonde conditionalization for
    // chi after diameter. This mirrors the production mathematical system from
    // public state and quadrature data, but it is not an internal solver trace.
    std::map<std::size_t, std::pair<double, double>> parent_by_path;
    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        if (quadrature.terminal_paths()[node].empty())
        {
            continue;
        }
        const std::size_t parent_index = quadrature.terminal_paths()[node].node_index_at_level(0);
        auto& entry = parent_by_path[parent_index];
        entry.first = quadrature.abscissas()[0][node];
        entry.second += quadrature.weights()[node];
    }

    ConditionalSolveDiagnostic diagnostic;
    for (const auto& [unused_path_index, node_and_weight] : parent_by_path)
    {
        static_cast<void>(unused_path_index);
        diagnostic.parent_nodes.push_back(node_and_weight.first);
        diagnostic.parent_weights.push_back(node_and_weight.second);
        diagnostic.parent_global_weights.push_back(node_and_weight.second);
    }
    const std::size_t parent_count = diagnostic.parent_nodes.size();
    diagnostic.vandermonde = Matrix(parent_count, parent_count);
    diagnostic.rhs = Matrix(parent_count, 6);
    diagnostic.solved = Matrix(parent_count, 6);
    for (std::size_t row = 0; row < parent_count; ++row)
    {
        for (std::size_t col = 0; col < parent_count; ++col)
        {
            diagnostic.vandermonde(row, col) = diagnostic.parent_weights[col] * std::pow(diagnostic.parent_nodes[col], static_cast<double>(row));
        }
        for (std::size_t chi_order = 0; chi_order < 6; ++chi_order)
        {
            qmomkit::MultiIndex index{row, chi_order, 0, 0, 0};
            const std::optional<double> value = moment_value(layout, values, index);
            if (!value.has_value())
            {
                throw std::runtime_error("The diagnostic CQMOM layout does not contain a required d/chi moment.");
            }
            diagnostic.rhs(row, chi_order) = *value;
        }
    }
    for (std::size_t chi_order = 0; chi_order < 6; ++chi_order)
    {
        std::vector<double> rhs(parent_count, 0.0);
        for (std::size_t row = 0; row < parent_count; ++row)
        {
            rhs[row] = diagnostic.rhs(row, chi_order);
        }
        const std::vector<double> solution = solve_linear_system(diagnostic.vandermonde, rhs);
        for (std::size_t row = 0; row < parent_count; ++row)
        {
            diagnostic.solved(row, chi_order) = solution[row];
        }
    }
    diagnostic.condition_estimate = infinity_norm(diagnostic.vandermonde) * inverse_infinity_norm(diagnostic.vandermonde);
    diagnostic.residual_norm = residual_infinity_norm(diagnostic.vandermonde, diagnostic.solved, diagnostic.rhs);
    diagnostic.maximum_input_scale = std::max(max_abs_value(diagnostic.vandermonde), max_abs_value(diagnostic.rhs));
    diagnostic.maximum_output_scale = max_abs_value(diagnostic.solved);
    return diagnostic;
}

std::vector<double> branch_moment_row(const ConditionalSolveDiagnostic& solve, const std::size_t parent_node_index)
{
    std::vector<double> moments(6, 0.0);
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        moments[order] = solve.solved(parent_node_index, order);
    }
    return moments;
}

std::optional<OffendingNode> find_first_out_of_support_chi(const qmomkit::CqmomQuadrature& quadrature)
{
    // The diagnostic stops on the first accepted post-source chi node outside
    // [0,1]. Waiting for a later continuous callback failure would conflate the
    // reconstruction event with downstream particle dynamics.
    constexpr double tolerance = 1.0e-12;
    for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
    {
        if (quadrature.terminal_paths()[node].size() <= 1)
        {
            continue;
        }
        const double chi = quadrature.abscissas()[chi_dimension][node];
        if (!std::isfinite(chi) || (chi >= chi_lower - tolerance && chi <= chi_upper + tolerance))
        {
            continue;
        }

        OffendingNode offending;
        offending.terminal_node = node;
        offending.path = quadrature.terminal_paths()[node];
        offending.terminal_weight = quadrature.weights()[node];
        for (std::size_t dimension = 0; dimension < quadrature.number_of_dimensions(); ++dimension)
        {
            offending.state.push_back(quadrature.abscissas()[dimension][node]);
        }

        const std::size_t parent_index = offending.path.node_index_at_level(0);
        std::map<std::size_t, std::pair<double, double>> chi_by_path;
        for (std::size_t terminal = 0; terminal < quadrature.weights().size(); ++terminal)
        {
            if (quadrature.terminal_paths()[terminal].size() <= 1)
            {
                continue;
            }
            if (quadrature.terminal_paths()[terminal].node_index_at_level(0) != parent_index)
            {
                continue;
            }
            offending.parent_global_weight += quadrature.weights()[terminal];
            const std::size_t chi_index = quadrature.terminal_paths()[terminal].node_index_at_level(1);
            auto& entry = chi_by_path[chi_index];
            entry.first = quadrature.abscissas()[chi_dimension][terminal];
            entry.second += quadrature.weights()[terminal];
        }
        for (const auto& [chi_index, node_and_weight] : chi_by_path)
        {
            static_cast<void>(chi_index);
            offending.chi_nodes.push_back(node_and_weight.first);
            offending.chi_global_weights.push_back(node_and_weight.second);
            offending.chi_conditional_weights.push_back(node_and_weight.second / offending.parent_global_weight);
        }
        const std::size_t offending_chi_index = offending.path.node_index_at_level(1);
        offending.chi_conditional_weight = chi_by_path.at(offending_chi_index).second / offending.parent_global_weight;
        return offending;
    }
    return std::nullopt;
}

std::optional<qmomkit::CqmomDiagnostics> find_branch_diagnostic(const std::vector<qmomkit::CqmomDiagnostics>& diagnostics,
                                                                const qmomkit::ConditioningPath& parent_path)
{
    for (const qmomkit::CqmomDiagnostics& diagnostic : diagnostics)
    {
        if (diagnostic.conditioning_level == 1 && diagnostic.phase_space_dimension == chi_dimension && same_path(diagnostic.parent_path, parent_path))
        {
            return diagnostic;
        }
    }
    return std::nullopt;
}

std::optional<qmomkit::CqmomNodeCountTransition> find_branch_transition(const std::vector<qmomkit::CqmomNodeCountTransition>& transitions,
                                                                        const qmomkit::ConditioningPath& parent_path)
{
    for (const qmomkit::CqmomNodeCountTransition& transition : transitions)
    {
        if (transition.stage == qmomkit::TransportStage::post_source_inversion && transition.conditioning_level == 1 &&
            transition.phase_space_dimension == chi_dimension && same_path(transition.parent_path, parent_path))
        {
            return transition;
        }
    }
    return std::nullopt;
}

qmomkit::CqmomMomentLayout layout_from_state(const qmomkit::CqmomState& state)
{
    return qmomkit::CqmomMomentLayout{state.hierarchy(), state.moment_indices()};
}

qmomkit::CqmomQuadrature reconstruct_stage_quadrature(const cqmom_non_spherical_particles::PreparedCqmomExample& example,
                                                      const qmomkit::CqmomMomentLayout& layout, const std::vector<double>& values)
{
    qmomkit::CqmomMomentData data{layout, values};
    qmomkit::CqmomLeanCase stage_case{std::move(data),
                                      example.cqmom_case.conditional_options(),
                                      example.cqmom_case.qmom_options(),
                                      example.cqmom_case.adaptive_branch_counts(),
                                      example.cqmom_case.diagnostic_level(),
                                      example.cqmom_case.vandermonde_solver(),
                                      example.cqmom_case.coordinate_support()};
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(stage_case);
    if (!result.ok())
    {
        throw std::runtime_error("Diagnostic stage CQMOM reconstruction failed: " + result.message);
    }
    return result.quadrature;
}

SourceReplay replay_source_stages(const cqmom_non_spherical_particles::PreparedCqmomExample& example, const qmomkit::CqmomState& pre_source_state,
                                  const qmomkit::ConditioningPath& parent_path, const double source_start_time, const double time_step)
{
    // QMOMKit's current moment-source update path is SSP-RK3. The replay uses
    // the same callback, state, stage times, and coefficients, then verifies the
    // final candidate against the public source-update API.
    SourceReplay replay;
    replay.layout = layout_from_state(pre_source_state);
    const std::vector<double> y0 = pre_source_state.moment_values();
    const auto rhs = [&](const double stage_time, const std::vector<double>& values)
    {
        const qmomkit::CqmomSourceContext context{stage_time, replay.layout->hierarchy, replay.layout->required_indices, values};
        qmomkit::CqmomMomentRate rate = example.source(stage_time, context);
        if (rate.values.size() != values.size())
        {
            throw std::runtime_error("Diagnostic source replay received a rate vector with the wrong size.");
        }
        return rate.values;
    };
    const auto branch_moments_for = [&](const std::vector<double>& values)
    {
        const qmomkit::CqmomQuadrature stage_quadrature = reconstruct_stage_quadrature(example, *replay.layout, values);
        const ConditionalSolveDiagnostic solve = conditional_chi_solve(stage_quadrature, *replay.layout, values);
        return branch_moment_row(solve, parent_path.node_index_at_level(0));
    };

    replay.stages.push_back({"pre_source", source_start_time, y0, branch_moments_for(y0)});
    const std::vector<double> k1 = rhs(source_start_time, y0);
    std::vector<double> stage1(y0.size(), 0.0);
    for (std::size_t i = 0; i < y0.size(); ++i)
    {
        stage1[i] = y0[i] + time_step * k1[i];
    }
    replay.stages.push_back({"stage1_candidate", source_start_time + time_step, stage1, branch_moments_for(stage1)});

    const std::vector<double> k2 = rhs(source_start_time + time_step, stage1);
    std::vector<double> embedded(y0.size(), 0.0);
    std::vector<double> stage2(y0.size(), 0.0);
    for (std::size_t i = 0; i < y0.size(); ++i)
    {
        embedded[i] = 0.5 * (y0[i] + (stage1[i] + time_step * k2[i]));
        stage2[i] = 0.75 * y0[i] + 0.25 * (stage1[i] + time_step * k2[i]);
    }
    replay.stages.push_back({"stage2_candidate", source_start_time + 0.5 * time_step, stage2, branch_moments_for(stage2)});
    replay.stages.push_back({"embedded_comparison_candidate", source_start_time + time_step, embedded, branch_moments_for(embedded)});

    const std::vector<double> k3 = rhs(source_start_time + 0.5 * time_step, stage2);
    std::vector<double> final(y0.size(), 0.0);
    for (std::size_t i = 0; i < y0.size(); ++i)
    {
        final[i] = (y0[i] + 2.0 * (stage2[i] + time_step * k3[i])) / 3.0;
    }
    replay.stages.push_back({"accepted_post_source_pre_inversion", source_start_time + time_step, final, branch_moments_for(final)});

    qmomkit::CqmomMomentData input_data{*replay.layout, y0};
    const qmomkit::CqmomMomentSourceUpdateResult public_update = qmomkit::advance_cqmom_moments_with_source(
        input_data, source_start_time, time_step, example.source, example.transport_options.source_options);
    if (!public_update.ok())
    {
        throw std::runtime_error("Public source update failed during diagnostic replay verification: " + public_update.message);
    }
    for (std::size_t i = 0; i < final.size(); ++i)
    {
        replay.max_difference_from_public_source_update =
            std::max(replay.max_difference_from_public_source_update, std::abs(final[i] - public_update.moments.values[i]));
    }
    if (replay.max_difference_from_public_source_update > 1.0e-11)
    {
        throw std::runtime_error("Diagnostic source replay does not match the public source update.");
    }
    return replay;
}

struct RecurrenceDiagnostics
{
    std::vector<double> a;
    std::vector<double> b;
    Matrix sigma;
};

RecurrenceDiagnostics recurrence_diagnostics(const std::vector<double>& moments, const std::size_t n)
{
    RecurrenceDiagnostics diagnostics{std::vector<double>(n, 0.0), std::vector<double>(n, 0.0), Matrix(2 * n + 1, 2 * n + 1)};
    for (std::size_t i = 1; i <= 2 * n; ++i)
    {
        diagnostics.sigma(1, i) = moments.at(i - 1);
    }
    diagnostics.a[0] = moments.at(1) / moments.at(0);
    diagnostics.b[0] = 0.0;
    for (std::size_t k = 2; k <= n; ++k)
    {
        for (std::size_t l = k; l <= 2 * n - k + 1; ++l)
        {
            diagnostics.sigma(k, l) = diagnostics.sigma(k - 1, l + 1) - diagnostics.a[k - 2] * diagnostics.sigma(k - 1, l) -
                                      diagnostics.b[k - 2] * diagnostics.sigma(k - 2, l);
        }
        diagnostics.a[k - 1] = diagnostics.sigma(k, k + 1) / diagnostics.sigma(k, k) - diagnostics.sigma(k - 1, k) / diagnostics.sigma(k - 1, k - 1);
        diagnostics.b[k - 1] = diagnostics.sigma(k, k) / diagnostics.sigma(k - 1, k - 1);
    }
    return diagnostics;
}

std::vector<double> reconstruct_1d_moments(const qmomkit::Quadrature1D& quadrature, const std::size_t count)
{
    std::vector<double> moments(count, 0.0);
    for (std::size_t order = 0; order < count; ++order)
    {
        for (std::size_t node = 0; node < quadrature.weights.size(); ++node)
        {
            moments[order] += quadrature.weights[node] * std::pow(quadrature.abscissas[node], static_cast<double>(order));
        }
    }
    return moments;
}

std::vector<double> reconstruct_1d_moments(const std::vector<double>& nodes, const std::vector<double>& weights, const std::size_t count)
{
    std::vector<double> moments(count, 0.0);
    for (std::size_t order = 0; order < count; ++order)
    {
        for (std::size_t node = 0; node < nodes.size(); ++node)
        {
            moments[order] += weights[node] * std::pow(nodes[node], static_cast<double>(order));
        }
    }
    return moments;
}

BoundedAcceptedBranch actual_bounded_branch_from_quadrature(const qmomkit::CqmomQuadrature& quadrature, const qmomkit::ConditioningPath& parent_path,
                                                            const std::vector<double>& rejected_candidate_moments, const std::size_t accepted_count)
{
    // The actual bounded closure is read from the accepted bounded quadrature.
    // This avoids presenting a diagnostic-only lower-count solve as production
    // behavior and preserves the bounded workflow's canonicalization contract.
    BoundedAcceptedBranch accepted;
    accepted.parent_path = parent_path;
    std::map<std::size_t, std::pair<double, double>> chi_by_path;

    for (std::size_t terminal = 0; terminal < quadrature.weights().size(); ++terminal)
    {
        const qmomkit::ConditioningPath& terminal_path = quadrature.terminal_paths()[terminal];
        if (terminal_path.size() <= 1 || terminal_path.node_index_at_level(0) != parent_path.node_index_at_level(0))
        {
            continue;
        }
        accepted.parent_global_weight += quadrature.weights()[terminal];
        const std::size_t chi_index = terminal_path.node_index_at_level(1);
        auto& entry = chi_by_path[chi_index];
        entry.first = quadrature.abscissas()[chi_dimension][terminal];
        entry.second += quadrature.weights()[terminal];
    }

    if (chi_by_path.empty() || accepted.parent_global_weight == 0.0)
    {
        throw std::runtime_error("Could not recover the accepted bounded chi branch from the bounded post-source quadrature.");
    }

    for (const auto& [unused_chi_index, node_and_weight] : chi_by_path)
    {
        static_cast<void>(unused_chi_index);
        accepted.chi_nodes.push_back(node_and_weight.first);
        accepted.global_weights.push_back(node_and_weight.second);
        accepted.conditional_weights.push_back(node_and_weight.second / accepted.parent_global_weight);
    }
    if (accepted.chi_nodes.size() != accepted_count)
    {
        throw std::runtime_error("The accepted bounded chi branch count does not match the bounded branch diagnostic.");
    }

    accepted.canonical_moments = reconstruct_1d_moments(accepted.chi_nodes, accepted.conditional_weights, rejected_candidate_moments.size());
    const std::size_t retained_moment_count = std::min(rejected_candidate_moments.size(), 2 * accepted_count);
    for (std::size_t order = 0; order < retained_moment_count; ++order)
    {
        const double absolute = std::abs(accepted.canonical_moments[order] - rejected_candidate_moments[order]);
        const double relative = absolute / std::max(std::abs(rejected_candidate_moments[order]), std::numeric_limits<double>::min());
        accepted.retained_prefix_absolute_residuals.push_back(absolute);
        accepted.retained_prefix_relative_residuals.push_back(relative);
    }
    return accepted;
}

BoundaryAdvance advance_to_source_boundary(const cqmom_non_spherical_particles::PreparedCqmomExample& example, const qmomkit::CqmomState& state,
                                           const double time, const double time_step)
{
    // Advance exactly to the source boundary for the configured splitting
    // scheme. The unbounded diagnostic inspects this boundary before any
    // following continuous stage can advect the accepted unsupported chi node.
    BoundaryAdvance advanced;
    advanced.source_start_time = time;
    const auto fail = [&](const qmomkit::CqmomTransportResult& result)
    {
        advanced.ok = false;
        advanced.message = result.message;
        advanced.failed_stage = result.diagnostics.failed_stage;
        if (result.quadrature.has_value())
        {
            advanced.final_quadrature = *result.quadrature;
        }
        return advanced;
    };

    if (example.transport_options.scheme == qmomkit::TransportScheme::strang)
    {
        const double half_step = 0.5 * time_step;
        qmomkit::CqmomTransportResult first = qmomkit::advance_cqmom_continuous(state, time, half_step, example.velocity, example.transport_options);
        if (!first.ok() || !first.quadrature.has_value())
        {
            return fail(first);
        }
        advanced.pre_source_state = first.state;
        advanced.pre_source_quadrature = *first.quadrature;
        advanced.source_start_time = time + half_step;
        qmomkit::CqmomTransportResult source = qmomkit::advance_cqmom_discontinuous(*advanced.pre_source_state, advanced.source_start_time, time_step,
                                                                                    example.source, example.transport_options);
        if (!source.ok() || !source.quadrature.has_value())
        {
            return fail(source);
        }
        advanced.source_diagnostics = source.diagnostics;
        advanced.boundary_state = source.state;
        advanced.boundary_quadrature = *source.quadrature;
        advanced.boundary_time = source.final_time;
        advanced.stage_label = "post_source_inversion_before_second_continuous";
        qmomkit::CqmomTransportResult second =
            qmomkit::advance_cqmom_continuous(*advanced.boundary_state, time + half_step, half_step, example.velocity, example.transport_options);
        if (!second.ok() || !second.quadrature.has_value())
        {
            return fail(second);
        }
        advanced.final_state = second.state;
        advanced.final_quadrature = *second.quadrature;
        advanced.ok = true;
        return advanced;
    }
    if (example.transport_options.scheme == qmomkit::TransportScheme::cd)
    {
        qmomkit::CqmomTransportResult continuous =
            qmomkit::advance_cqmom_continuous(state, time, time_step, example.velocity, example.transport_options);
        if (!continuous.ok() || !continuous.quadrature.has_value())
        {
            return fail(continuous);
        }
        advanced.pre_source_state = continuous.state;
        advanced.pre_source_quadrature = *continuous.quadrature;
        qmomkit::CqmomTransportResult source =
            qmomkit::advance_cqmom_discontinuous(*advanced.pre_source_state, time, time_step, example.source, example.transport_options);
        if (!source.ok() || !source.quadrature.has_value())
        {
            return fail(source);
        }
        advanced.source_diagnostics = source.diagnostics;
        advanced.boundary_state = source.state;
        advanced.boundary_quadrature = *source.quadrature;
        advanced.final_state = source.state;
        advanced.final_quadrature = *source.quadrature;
        advanced.boundary_time = source.final_time;
        advanced.stage_label = "post_source_inversion";
        advanced.ok = true;
        return advanced;
    }

    qmomkit::CqmomTransportResult source = qmomkit::advance_cqmom_discontinuous(state, time, time_step, example.source, example.transport_options);
    if (!source.ok() || !source.quadrature.has_value())
    {
        return fail(source);
    }
    advanced.pre_source_state = state;
    advanced.source_diagnostics = source.diagnostics;
    advanced.boundary_state = source.state;
    advanced.boundary_quadrature = *source.quadrature;
    advanced.boundary_time = source.final_time;
    advanced.stage_label = "post_source_inversion_before_continuous";
    qmomkit::CqmomTransportResult continuous =
        qmomkit::advance_cqmom_continuous(*advanced.boundary_state, time, time_step, example.velocity, example.transport_options);
    if (!continuous.ok() || !continuous.quadrature.has_value())
    {
        return fail(continuous);
    }
    advanced.final_state = continuous.state;
    advanced.final_quadrature = *continuous.quadrature;
    advanced.ok = true;
    return advanced;
}

void write_failure_context(const std::filesystem::path& directory, const qmomkit::CqmomQuadrature& last_quadrature,
                           const cqmom_non_spherical_particles::PreparedCqmomExample& example, const double time, const std::string& message)
{
    cqmom_non_spherical_particles::ensure_output_directory(directory);
    std::ofstream report(directory / "report.md");
    report << std::setprecision(std::numeric_limits<double>::max_digits10);
    double min_d = std::numeric_limits<double>::infinity();
    double max_abs_velocity = 0.0;
    double max_rate = 0.0;
    for (std::size_t node = 0; node < last_quadrature.weights().size(); ++node)
    {
        std::vector<double> state;
        for (std::size_t dimension = 0; dimension < last_quadrature.number_of_dimensions(); ++dimension)
        {
            state.push_back(last_quadrature.abscissas()[dimension][node]);
        }
        min_d = std::min(min_d, state[0]);
        max_abs_velocity = std::max({max_abs_velocity, std::abs(state[2]), std::abs(state[3]), std::abs(state[4])});
        const std::vector<double> rate = example.velocity(time, state);
        for (const double component : rate)
        {
            max_rate = std::max(max_rate, std::abs(component));
        }
    }
    report << "# Unbounded Diagnostic Failure Context\n\n";
    report << "The diagnostic stopped before reproducing an accepted out-of-support chi node.\n\n";
    report << "- time: " << time << '\n';
    report << "- message: " << message << '\n';
    report << "- smallest diameter: " << min_d << '\n';
    report << "- largest absolute velocity coordinate: " << max_abs_velocity << '\n';
    report << "- largest phase-space rate component: " << max_rate << '\n';
    report << "\nThe last accepted quadrature contained " << last_quadrature.weights().size() << " terminal nodes.\n";
}

void write_moment_table(std::ostream& output, const std::vector<double>& moments, const double branch_weight)
{
    output << "| order | conditional moment | branch-weighted moment |\n";
    output << "|---:|---:|---:|\n";
    for (std::size_t order = 0; order < moments.size(); ++order)
    {
        output << "| " << order << " | " << moments[order] << " | " << branch_weight * moments[order] << " |\n";
    }
}

void write_matrix_table(std::ostream& output, const std::vector<MatrixDiagnostic>& matrices, const std::string& family)
{
    output << "| matrix | entries | determinant | eigenvalues | min eigenvalue | infinity norm | effective tolerance | signed margin | "
              "classification |\n";
    output << "|---|---|---:|---|---:|---:|---:|---:|---|\n";
    for (const MatrixDiagnostic& matrix : matrices)
    {
        if (matrix.family != family)
        {
            continue;
        }
        output << "| " << family << "(k=" << matrix.k << ", l=" << matrix.l << ") | `" << matrix_string(matrix.matrix) << "` | " << matrix.determinant
               << " | `" << vector_string(matrix.eigenvalues) << "` | " << matrix.minimum_eigenvalue << " | " << matrix.infinity_norm << " | "
               << matrix.effective_tolerance << " | " << matrix.signed_margin << " | " << matrix.classification << " |\n";
    }
}

void write_interval_matrix_table(std::ostream& output, const std::vector<MatrixDiagnostic>& matrices, const std::string& family,
                                 const std::string& label)
{
    output << "| matrix | entries | determinant | eigenvalues | min eigenvalue | infinity norm | effective tolerance | signed margin | "
              "classification |\n";
    output << "|---|---|---:|---|---:|---:|---:|---:|---|\n";
    for (const MatrixDiagnostic& matrix : matrices)
    {
        if (matrix.family != family)
        {
            continue;
        }
        output << "| " << label << "(l=" << matrix.l << ") | `" << matrix_string(matrix.matrix) << "` | " << matrix.determinant << " | `"
               << vector_string(matrix.eigenvalues) << "` | " << matrix.minimum_eigenvalue << " | " << matrix.infinity_norm << " | "
               << matrix.effective_tolerance << " | " << matrix.signed_margin << " | " << matrix.classification << " |\n";
    }
}

void write_report_file(const std::filesystem::path& directory, const DiagnosticReport& report,
                       const cqmom_non_spherical_particles::PreparedCqmomExample& unbounded,
                       const cqmom_non_spherical_particles::PreparedCqmomExample& bounded)
{
    cqmom_non_spherical_particles::ensure_output_directory(directory);
    const auto precision = std::numeric_limits<double>::max_digits10;
    std::ofstream markdown(directory / "report.md");
    markdown << std::setprecision(precision);

    markdown << "# Unbounded conditional-support diagnostic\n\n";
    markdown << "This diagnostic intentionally omits the configured chi support [0,1], accepts the first unconstrained post-source reconstruction "
                "with chi outside that interval, writes this report, and stops before the next continuous stage.\n\n";

    markdown << "## Run configuration\n\n";
    markdown << "| quantity | bounded workflow | unbounded diagnostic |\n";
    markdown << "|---|---|---|\n";
    markdown << "| configuration | `" << bounded.configuration_file.string() << "` | `" << unbounded.configuration_file.string() << "` |\n";
    markdown << "| sampled initialization | seed " << cqmom_non_spherical_particles::initialization_seed << ", samples "
             << cqmom_non_spherical_particles::initialization_sample_count << " | same seed and sample count |\n";
    markdown << "| node counts | `" << vector_string(unbounded.cqmom_case.hierarchy().node_counts()) << "` | same |\n";
    markdown << "| conditioning order | `" << vector_string(unbounded.cqmom_case.hierarchy().conditioning_order()) << "` | same |\n";
    markdown << "| transport scheme | `" << transport_scheme_string(unbounded.transport_options.scheme) << "` | same |\n";
    markdown << "| macro time step | " << cqmom_non_spherical_particles::cqmom_time_step << " | same |\n";
    markdown << "| source model | binary fragmentation callback from shared setup | same callback construction |\n";
    markdown << "| chi support declaration | `[0,1]` | intentionally absent |\n\n";
    markdown
        << "The two CQMOM workflows are prepared by the same example-local helper, so the sampled realization, physical parameters, node counts, "
           "conditioning order, integration settings, tolerances, and source model match except for the chi support declaration and output "
           "directory.\n\n";

    markdown << "## Detection event\n\n";
    markdown << "| quantity | value |\n";
    markdown << "|---|---:|\n";
    markdown << "| macro step | " << report.macro_step << " |\n";
    markdown << "| physical time | " << report.time << " |\n";
    markdown << "| operator stage | " << report.operator_stage_label << " |\n";
    markdown << "| coordinate index | " << chi_dimension << " |\n";
    markdown << "| expected support interval | [" << chi_lower << ", " << chi_upper << "] |\n";
    markdown << "| terminal node | " << report.offending.terminal_node << " |\n";
    markdown << "| complete terminal path | `" << path_string(report.offending.path) << "` |\n";
    markdown << "| parent branch path | `" << path_string(report.unbounded_branch.parent_path) << "` |\n";
    markdown << "| requested count | " << report.unbounded_branch.requested_node_count << " |\n";
    markdown << "| accepted count | " << report.unbounded_branch.actual_node_count << " |\n";
    markdown << "| parent/global branch weight | " << report.offending.parent_global_weight << " |\n";
    markdown << "| offending terminal-node weight | " << report.offending.terminal_weight << " |\n";
    markdown << "| offending conditional chi weight | " << report.offending.chi_conditional_weight << " |\n\n";
    markdown << "Complete offending state `(d, chi, u, v, w)`: `" << vector_string(report.offending.state) << "`\n\n";
    markdown << "Unbounded chi nodes: `" << vector_string(report.offending.chi_nodes) << "`\n\n";
    markdown << "Unbounded conditional chi weights: `" << vector_string(report.offending.chi_conditional_weights) << "`\n\n";
    markdown << "Unbounded corresponding global weights: `" << vector_string(report.offending.chi_global_weights) << "`\n\n";

    markdown << "## Physical and hierarchy context\n\n";
    markdown
        << "The phase-space ordering is `(d, chi, u, v, w)`. The coordinate name `chi` is example metadata, and this model communicates its physical "
           "aspect-ratio interval through the generic coordinate-support declaration. The diagnostic inspects the accepted post-source quadrature "
           "immediately after moment-space fragmentation because the next continuous stage could obscure the closure event that produced the "
           "out-of-support value.\n\n";

    markdown << "## Conditional moment sequence\n\n";
    markdown << "The unbounded workflow supplies all six branch-local chi moments to the three-node inversion:\n\n";
    write_moment_table(markdown, report.branch_moments, report.offending.parent_global_weight);
    markdown << '\n';

    markdown << "## Source-stage evolution\n\n";
    markdown
        << "The source replay uses the same pre-source moment state, the same source callback, and the SSP-RK3 stage times and coefficients used by "
           "`advance_cqmom_moments_with_source(...)`. QMOMKit currently exposes this discontinuous source update through SSP-RK3 options; this "
           "example-local diagnostic should be extended before being used with any future selectable source scheme. The replay final candidate "
           "agrees with the public source update to "
        << report.source_replay.max_difference_from_public_source_update << ".\n\n";
    markdown << "| stage | time | m0 | m1 | m2 | m3 | m4 | m5 |\n";
    markdown << "|---|---:|---:|---:|---:|---:|---:|---:|\n";
    for (const StageMoments& stage : report.source_replay.stages)
    {
        markdown << "| " << stage.label << " | " << stage.time;
        for (const double moment : stage.branch_moments)
        {
            markdown << " | " << moment;
        }
        markdown << " |\n";
    }
    markdown << '\n';

    markdown << "## Classical Hankel-Hadamard determinants\n\n";
    markdown << "Classical determinants use `H(k,l) = [m_(k+i+j)]` for `i,j = 0,...,l` and `Delta(k,l) = det(H(k,l))`. "
                "For six moments the relevant classical checks are only the k = 0 and k = 1 Hankel families below.\n\n";
    markdown << "| determinant | value |\n";
    markdown << "|---|---:|\n";
    for (const std::size_t l : {std::size_t{0}, std::size_t{1}, std::size_t{2}})
    {
        for (const std::size_t k : {std::size_t{0}, std::size_t{1}})
        {
            for (const MatrixDiagnostic& matrix : report.matrices)
            {
                if (matrix.family == "classical_hankel_hadamard" && matrix.k == k && matrix.l == l)
                {
                    markdown << "| Delta(" << k << "," << l << ") | " << matrix.determinant << " |\n";
                }
            }
        }
    }
    markdown << '\n';

    markdown << "## Interval localizing matrices\n\n";
    markdown << "For interval `[a,b]`, `Delta_lower(l) = det([m_(i+j+1) - a*m_(i+j)])` and "
                "`Delta_upper(l) = det([b*m_(i+j) - m_(i+j+1)])`. For `[0,1]`, the upper matrix is `[m_(i+j) - m_(i+j+1)]`, "
                "which detects the missing upper chi bound that the classical Stieltjes conditions do not enforce.\n\n";
    markdown << "### Lower interval-localizing determinants\n\n";
    markdown << "| determinant | value |\n";
    markdown << "|---|---:|\n";
    for (const MatrixDiagnostic& matrix : report.matrices)
    {
        if (matrix.family == "lower_localizing")
        {
            markdown << "| Delta_lower(" << matrix.l << ") | " << matrix.determinant << " |\n";
        }
    }
    markdown << "\n### Upper interval-localizing determinants\n\n";
    markdown << "| determinant | value |\n";
    markdown << "|---|---:|\n";
    for (const MatrixDiagnostic& matrix : report.matrices)
    {
        if (matrix.family == "upper_localizing")
        {
            markdown << "| Delta_upper(" << matrix.l << ") | " << matrix.determinant << " |\n";
        }
    }
    markdown << '\n';

    markdown << "## Eigenvalue and tolerance analysis\n\n";
    markdown << "Determinant signs are explanatory. Positive-semidefinite acceptance is based on the minimum eigenvalue compared with "
                "`1e-12 * max(1, infinity_norm)`, matching the tolerance scaling used by QMOMKit realizability checks.\n\n";
    markdown << "### Classical Hankel-Hadamard matrices\n\n";
    write_matrix_table(markdown, report.matrices, "classical_hankel_hadamard");
    markdown << "\n### Lower interval-localizing matrices\n\n";
    write_interval_matrix_table(markdown, report.matrices, "lower_localizing", "Delta_lower");
    markdown << "\n### Upper interval-localizing matrices\n\n";
    write_interval_matrix_table(markdown, report.matrices, "upper_localizing", "Delta_upper");
    markdown << '\n';

    markdown << "## Weighted Vandermonde reconstruction\n\n";
    markdown
        << "The diagnostic reproduces the production weighted Vandermonde system and verifies numerical agreement; it is not a direct internal trace "
           "of private solver storage.\n\n";
    markdown << "- parent marginal nodes: `" << vector_string(report.conditional_solve.parent_nodes) << "`\n";
    markdown << "- parent local weights: `" << vector_string(report.conditional_solve.parent_weights) << "`\n";
    markdown << "- parent global weights: `" << vector_string(report.conditional_solve.parent_global_weights) << "`\n";
    markdown << "- accumulated offending branch weight: " << report.offending.parent_global_weight << '\n';
    markdown << "- weighted Vandermonde matrix: `" << matrix_string(report.conditional_solve.vandermonde) << "`\n";
    markdown << "- RHS values by row and chi order: `" << matrix_string(report.conditional_solve.rhs) << "`\n";
    markdown << "- solved conditional moments by parent row and chi order: `" << matrix_string(report.conditional_solve.solved) << "`\n";
    markdown << "- condition estimate: " << report.conditional_solve.condition_estimate << '\n';
    markdown << "- residual infinity norm: " << report.conditional_solve.residual_norm << '\n';
    markdown << "- maximum RHS/input magnitude: " << report.conditional_solve.maximum_input_scale << '\n';
    markdown << "- maximum conditional-solution magnitude: " << report.conditional_solve.maximum_output_scale << "\n\n";

    markdown << "## Wheeler reconstruction\n\n";
    markdown << "- supplied conditional moments: `" << vector_string(report.branch_moments) << "`\n";
    markdown << "- recurrence a coefficients: `" << vector_string(report.recurrence_a) << "`\n";
    markdown << "- recurrence b coefficients: `" << vector_string(report.recurrence_b) << "`\n";
    std::vector<double> sigma_diagonal;
    for (std::size_t i = 1; i < report.sigma.rows; ++i)
    {
        sigma_diagonal.push_back(report.sigma(i, i));
    }
    markdown << "- sigma diagonal: `" << vector_string(sigma_diagonal) << "`\n";
    markdown << "- returned nodes: `" << vector_string(report.wheeler.abscissas) << "`\n";
    markdown << "- conditional weights: `" << vector_string(report.wheeler.weights) << "`\n";
    const auto minmax_weight = std::minmax_element(report.wheeler.weights.begin(), report.wheeler.weights.end());
    markdown << "- local weight ratio: " << (*minmax_weight.first / *minmax_weight.second) << '\n';
    markdown << "- configured rmin: " << report.configured_rmin << '\n';
    markdown << "- corresponding global weights: `" << vector_string(report.offending.chi_global_weights) << "`\n\n";
    markdown << "| order | supplied moment | reconstructed moment | absolute residual | relative residual |\n";
    markdown << "|---:|---:|---:|---:|---:|\n";
    for (std::size_t order = 0; order < report.branch_moments.size(); ++order)
    {
        const double absolute = std::abs(report.reconstructed_moments[order] - report.branch_moments[order]);
        const double relative = absolute / std::max(std::abs(report.branch_moments[order]), std::numeric_limits<double>::min());
        markdown << "| " << order << " | " << report.branch_moments[order] << " | " << report.reconstructed_moments[order] << " | " << absolute
                 << " | " << relative << " |\n";
    }
    markdown
        << "\nWheeler accurately reconstructs a Stieltjes-realizable sequence that is infeasible on `[0,1]`. Wheeler did not create the unsupported "
           "sequence; the unsupported node is the result of supplying an unbounded conditional moment sequence for a bounded physical "
           "coordinate.\n\n";

    markdown << "## Bounded workflow comparison\n\n";
    if (report.bounded_branch.has_value() && report.bounded_accepted_branch.has_value())
    {
        const qmomkit::CqmomDiagnostics& bounded_branch = *report.bounded_branch;
        const BoundedAcceptedBranch& accepted = *report.bounded_accepted_branch;
        markdown << "| quantity | value |\n";
        markdown << "|---|---:|\n";
        markdown << "| actual bounded branch path | `" << path_string(accepted.parent_path) << "` |\n";
        markdown << "| actual rejected candidate count | " << bounded_branch.bounded_candidate_node_count << " |\n";
        markdown << "| rejected candidate attempt type | " << bounded_branch.bounded_candidate_attempt_type << " |\n";
        markdown << "| actual support margin | " << bounded_branch.bounded_minimum_margin << " |\n";
        markdown << "| failed support matrix | " << bounded_branch.bounded_failed_matrix << " |\n";
        markdown << "| actual accepted node count | " << accepted.chi_nodes.size() << " |\n";
        if (report.bounded_transition.has_value())
        {
            markdown << "| active count before | "
                     << (report.bounded_transition->active_node_count_before.has_value()
                             ? std::to_string(*report.bounded_transition->active_node_count_before)
                             : std::string("none"))
                     << " |\n";
            markdown << "| active count after | "
                     << (report.bounded_transition->active_node_count_after.has_value()
                             ? std::to_string(*report.bounded_transition->active_node_count_after)
                             : std::string("none"))
                     << " |\n";
            markdown << "| active-count transition | " << transition_kind_string(report.bounded_transition->kind) << " |\n";
        }
        else if (bounded_branch.bounded_previous_active_node_count.has_value())
        {
            markdown << "| active count before | " << *bounded_branch.bounded_previous_active_node_count << " |\n";
            markdown << "| active count after | " << accepted.chi_nodes.size() << " |\n";
        }
        markdown << "| actual accepted chi nodes | `" << vector_string(accepted.chi_nodes) << "` |\n";
        markdown << "| actual accepted conditional weights | `" << vector_string(accepted.conditional_weights) << "` |\n";
        markdown << "| actual accepted global weights | `" << vector_string(accepted.global_weights) << "` |\n";
        markdown << "| parent/global branch weight | " << accepted.parent_global_weight << " |\n\n";
        markdown << "Actual rejected six-moment candidate: `" << vector_string(bounded_branch.bounded_candidate_moments) << "`\n\n";
        markdown << "Actual canonicalized moments associated with the accepted bounded quadrature: `" << vector_string(accepted.canonical_moments)
                 << "`\n\n";
        markdown << "| retained prefix order | candidate moment | canonicalized moment | absolute residual | relative residual |\n";
        markdown << "|---:|---:|---:|---:|---:|\n";
        for (std::size_t order = 0; order < accepted.retained_prefix_absolute_residuals.size(); ++order)
        {
            markdown << "| " << order << " | " << bounded_branch.bounded_candidate_moments[order] << " | " << accepted.canonical_moments[order]
                     << " | " << accepted.retained_prefix_absolute_residuals[order] << " | " << accepted.retained_prefix_relative_residuals[order]
                     << " |\n";
        }
        markdown << '\n';
    }
    else
    {
        markdown << "No matching bounded branch diagnostic and accepted bounded quadrature were available at this boundary.\n\n";
    }
    markdown
        << "The bounded workflow rejects the support-infeasible six-moment three-node candidate, accepts a support-valid lower-count branch, "
           "rebuilds descendants, and canonicalizes the accepted state from that quadrature. This is not clipping and not projection of all six "
           "moments onto a bounded moment cone; it is adaptive reduction of represented branch order. Higher-order source-candidate information is "
           "therefore replaced by the closure implied by the accepted lower-count quadrature.\n\n";

    markdown << "## Interpretation\n\n";
    markdown << "The unbounded workflow uses all six branch moments and ordinary Wheeler reconstruction returns three positive nodes, one above one, "
                "with moment residuals at roundoff. The bounded workflow keeps the physical chi support active, so the same source-boundary event is "
                "handled by support-aware adaptive reconstruction rather than by altering the particle model or clipping node coordinates.\n\n";

    markdown << "## Conclusion\n\n";
    markdown
        << "The first support-loss event is a conditional finite-moment closure issue: classical Hankel-Hadamard positivity is insufficient for the "
           "bounded interval `[0,1]`, and the upper interval-localizing matrix is the first failing condition. The bounded workflow resolves the "
           "event through adaptive branch-order reduction while preserving ordinary user-facing simulation output.\n";
}

DiagnosticReport build_report(const cqmom_non_spherical_particles::PreparedCqmomExample& unbounded,
                              const std::optional<qmomkit::CqmomDiagnostics>& bounded_branch, const BoundaryAdvance& unbounded_boundary,
                              const BoundaryAdvance& bounded_boundary, OffendingNode offending, const std::size_t macro_step)
{
    if (!unbounded_boundary.pre_source_state.has_value() || !unbounded_boundary.boundary_quadrature.has_value())
    {
        throw std::runtime_error("Diagnostic boundary state is incomplete.");
    }
    const qmomkit::ConditioningPath parent_path{{offending.path.node_index_at_level(0)}};
    const std::optional<qmomkit::CqmomDiagnostics> unbounded_branch =
        find_branch_diagnostic(unbounded_boundary.source_diagnostics.post_source_cqmom_diagnostics, parent_path);
    if (!unbounded_branch.has_value())
    {
        throw std::runtime_error("Could not find the offending unbounded branch diagnostic.");
    }

    SourceReplay replay;
    try
    {
        replay = replay_source_stages(unbounded, *unbounded_boundary.pre_source_state, parent_path, unbounded_boundary.source_start_time,
                                      cqmom_non_spherical_particles::cqmom_time_step);
    }
    catch (const std::exception& error)
    {
        throw std::runtime_error(std::string("source replay failed: ") + error.what());
    }
    ConditionalSolveDiagnostic solve;
    try
    {
        solve = conditional_chi_solve(*unbounded_boundary.boundary_quadrature, *replay.layout, replay.stages.back().values);
    }
    catch (const std::exception& error)
    {
        throw std::runtime_error(std::string("conditional chi solve failed: ") + error.what());
    }
    const std::vector<double> branch_moments = branch_moment_row(solve, offending.path.node_index_at_level(0));

    qmomkit::WheelerOptions wheeler_options;
    wheeler_options.adaptive = false;
    wheeler_options.rmin = unbounded.cqmom_case.qmom_options().rmin;
    wheeler_options.eabs = unbounded.cqmom_case.qmom_options().eabs;
    wheeler_options.cutoff = unbounded.cqmom_case.qmom_options().cutoff;
    const qmomkit::Quadrature1D wheeler = qmomkit::compute_quadrature(qmomkit::MomentSet(branch_moments), 3, wheeler_options);
    if (!wheeler.ok())
    {
        throw std::runtime_error("Diagnostic Wheeler verification failed: " + wheeler.message);
    }
    const RecurrenceDiagnostics recurrence = recurrence_diagnostics(branch_moments, 3);

    DiagnosticReport report;
    report.macro_step = macro_step;
    report.time = unbounded_boundary.boundary_time;
    report.operator_stage_label = unbounded_boundary.stage_label;
    report.offending = std::move(offending);
    report.unbounded_branch = *unbounded_branch;
    report.bounded_branch = bounded_branch;
    if (bounded_branch.has_value() && bounded_boundary.boundary_quadrature.has_value())
    {
        // The comparison must use the bounded workflow's accepted post-source
        // quadrature, not a fresh lower-count solve on the rejected candidate.
        report.bounded_transition = find_branch_transition(bounded_boundary.source_diagnostics.cqmom_node_count_transitions, parent_path);
        report.bounded_accepted_branch =
            actual_bounded_branch_from_quadrature(*bounded_boundary.boundary_quadrature, parent_path, bounded_branch->bounded_candidate_moments,
                                                  bounded_branch->bounded_accepted_fallback_node_count);
    }
    report.conditional_solve = solve;
    report.source_replay = std::move(replay);
    report.branch_moments = branch_moments;
    report.wheeler = wheeler;
    report.configured_rmin = unbounded.cqmom_case.qmom_options().rmin;
    report.reconstructed_moments = reconstruct_1d_moments(wheeler, branch_moments.size());
    report.recurrence_a = recurrence.a;
    report.recurrence_b = recurrence.b;
    report.sigma = recurrence.sigma;
    report.matrices = build_matrix_diagnostics(branch_moments);
    return report;
}

} // namespace

int main()
{
    const std::filesystem::path unbounded_configuration_file = "configuration_unbounded.qmom";
    const std::filesystem::path bounded_configuration_file = "configuration.qmom";

    try
    {
        std::optional<cqmom_non_spherical_particles::PreparedCqmomExample> unbounded_prepared =
            cqmom_non_spherical_particles::prepare_cqmom_example(unbounded_configuration_file);
        std::optional<cqmom_non_spherical_particles::PreparedCqmomExample> bounded_prepared =
            cqmom_non_spherical_particles::prepare_cqmom_example(bounded_configuration_file);
        if (!unbounded_prepared.has_value() || !bounded_prepared.has_value())
        {
            return 1;
        }
        cqmom_non_spherical_particles::PreparedCqmomExample unbounded = std::move(*unbounded_prepared);
        cqmom_non_spherical_particles::PreparedCqmomExample bounded = std::move(*bounded_prepared);

        std::cout << "5D CQMOM unbounded chi diagnostic\n";
        std::cout << "configuration: " << unbounded_configuration_file.string() << '\n';
        cqmom_non_spherical_particles::print_transport_settings(unbounded);

        qmomkit::CqmomState unbounded_state = std::move(unbounded.state);
        qmomkit::CqmomState bounded_state = std::move(bounded.state);
        qmomkit::CqmomQuadrature last_unbounded_quadrature = unbounded.initial_quadrature;
        double time = cqmom_non_spherical_particles::cqmom_initial_time;
        std::size_t macro_step = 0;

        while (time < cqmom_non_spherical_particles::cqmom_final_time - 1.0e-12)
        {
            ++macro_step;
            const double current_step =
                std::min(cqmom_non_spherical_particles::cqmom_time_step, cqmom_non_spherical_particles::cqmom_final_time - time);
            BoundaryAdvance unbounded_boundary = advance_to_source_boundary(unbounded, unbounded_state, time, current_step);
            BoundaryAdvance bounded_boundary = advance_to_source_boundary(bounded, bounded_state, time, current_step);
            if (!unbounded_boundary.ok)
            {
                write_failure_context(unbounded.output_options.directory, last_unbounded_quadrature, unbounded, time, unbounded_boundary.message);
                std::cerr << "Unbounded diagnostic transport failed before an out-of-support chi node was accepted: " << unbounded_boundary.message
                          << '\n';
                return 1;
            }
            if (!bounded_boundary.ok)
            {
                write_failure_context(unbounded.output_options.directory, last_unbounded_quadrature, unbounded, time, bounded_boundary.message);
                std::cerr << "Bounded comparison transport failed before the diagnostic event: " << bounded_boundary.message << '\n';
                return 1;
            }
            const std::optional<OffendingNode> offending = find_first_out_of_support_chi(*unbounded_boundary.boundary_quadrature);
            if (offending.has_value())
            {
                std::optional<DiagnosticReport> report;
                try
                {
                    const qmomkit::ConditioningPath parent_path{{offending->path.node_index_at_level(0)}};
                    const std::optional<qmomkit::CqmomDiagnostics> bounded_branch =
                        find_branch_diagnostic(bounded_boundary.source_diagnostics.post_source_cqmom_diagnostics, parent_path);
                    report = build_report(unbounded, bounded_branch, unbounded_boundary, bounded_boundary, *offending, macro_step);
                    write_report_file(unbounded.output_options.directory, *report, unbounded, bounded);
                }
                catch (const std::exception& error)
                {
                    throw std::runtime_error(std::string("Failed while building unbounded diagnostic report for path ") +
                                             path_string(offending->path) + ": " + error.what());
                }
                std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << "Detected accepted chi node outside [0,1] at time "
                          << report->time << ", path " << path_string(report->offending.path) << ", chi = " << report->offending.state[chi_dimension]
                          << ".\n";
                std::cout << "Diagnostic report written to " << (unbounded.output_options.directory / "report.md").string() << '\n';
                return 0;
            }

            unbounded_state = std::move(*unbounded_boundary.final_state);
            bounded_state = std::move(*bounded_boundary.final_state);
            last_unbounded_quadrature = *unbounded_boundary.final_quadrature;
            time += current_step;
        }

        std::cout << "No accepted out-of-support chi node was reproduced before final time.\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
