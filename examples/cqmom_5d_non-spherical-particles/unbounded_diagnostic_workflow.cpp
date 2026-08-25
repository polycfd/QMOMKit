#include "unbounded_diagnostic_workflow.hpp"

#include "example_workflow.hpp"
#include "qmomkit/moment_reconstruction.hpp"
#include "qmomkit/wheeler.hpp"
#include "sampling.hpp"
#include "unbounded_diagnostic_report.hpp"
#include "unbounded_diagnostic_types.hpp"
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

using cqmom_5d_non_spherical_particles::BoundaryAdvance;
using cqmom_5d_non_spherical_particles::BoundedAcceptedBranch;
using cqmom_5d_non_spherical_particles::ConditionalSolveDiagnostic;
using cqmom_5d_non_spherical_particles::DiagnosticReport;
using cqmom_5d_non_spherical_particles::EigenResult;
using cqmom_5d_non_spherical_particles::Matrix;
using cqmom_5d_non_spherical_particles::MatrixDiagnostic;
using cqmom_5d_non_spherical_particles::OffendingNode;
using cqmom_5d_non_spherical_particles::SourceReplay;
using cqmom_5d_non_spherical_particles::StageMoments;

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

std::string phase_space_evolution_scheme_string(const qmomkit::PhaseSpaceEvolutionScheme scheme)
{
    switch (scheme)
    {
    case qmomkit::PhaseSpaceEvolutionScheme::strang:
        return "strang";
    case qmomkit::PhaseSpaceEvolutionScheme::advsrc:
        return "advsrc";
    case qmomkit::PhaseSpaceEvolutionScheme::srcadv:
        return "srcadv";
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
    // [0,1]. Waiting for a later advection callback failure would conflate the
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
        if (transition.stage == qmomkit::PhaseSpaceEvolutionStage::post_source_inversion && transition.conditioning_level == 1 &&
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

SourceReplay record_phase_space_source_stages(const BoundaryAdvance& boundary, const qmomkit::ConditioningPath& parent_path)
{
    if (!boundary.pre_source_state.has_value() || !boundary.boundary_state.has_value() || !boundary.boundary_quadrature.has_value())
    {
        throw std::runtime_error("The phase-space source-stage record is incomplete.");
    }

    // The diagnostic records the states accepted by the same phase-space source
    // operation that evaluates the breakage callback at its SSP-RK stages.
    SourceReplay replay;
    replay.layout = layout_from_state(*boundary.boundary_state);
    const auto record = [&](const char* label, const double time, const qmomkit::CqmomState& state, const qmomkit::CqmomQuadrature& quadrature)
    {
        const qmomkit::CqmomMomentLayout layout = layout_from_state(state);
        const ConditionalSolveDiagnostic solve = conditional_chi_solve(quadrature, layout, state.moment_values());
        replay.stages.push_back({label, time, state.moment_values(), branch_moment_row(solve, parent_path.node_index_at_level(0))});
    };

    if (boundary.pre_source_quadrature.has_value())
    {
        record("pre_source", boundary.source_start_time, *boundary.pre_source_state, *boundary.pre_source_quadrature);
    }
    else
    {
        replay.stages.push_back({"pre_source_quadrature_unavailable", boundary.source_start_time, boundary.pre_source_state->moment_values(), {}});
    }
    record("accepted_post_source", boundary.boundary_time, *boundary.boundary_state, *boundary.boundary_quadrature);
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

BoundedAcceptedBranch actual_bounded_branch_from_quadrature(const qmomkit::CqmomQuadrature& quadrature, const qmomkit::ConditioningPath& parent_path,
                                                            const std::size_t accepted_count)
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

    return accepted;
}

BoundaryAdvance advance_to_source_boundary(const cqmom_5d_non_spherical_particles::PreparedCqmomExample& example, const qmomkit::CqmomState& state,
                                           const double time, const double time_step)
{
    // Advance exactly to the source boundary for the configured splitting
    // scheme. The unbounded diagnostic inspects this boundary before any
    // following advection stage can advect the accepted unsupported chi node.
    BoundaryAdvance advanced;
    advanced.source_start_time = time;
    const auto fail = [&](const qmomkit::CqmomPhaseSpaceEvolutionResult& result)
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

    if (example.phase_space_evolution_options.scheme == qmomkit::PhaseSpaceEvolutionScheme::strang)
    {
        const double half_step = 0.5 * time_step;
        qmomkit::CqmomPhaseSpaceEvolutionResult first =
            qmomkit::advance_cqmom_advection(state, time, half_step, example.velocity, example.phase_space_evolution_options);
        if (!first.ok() || !first.quadrature.has_value())
        {
            return fail(first);
        }
        advanced.pre_source_state = first.state;
        advanced.pre_source_quadrature = *first.quadrature;
        advanced.source_start_time = time + half_step;
        qmomkit::CqmomPhaseSpaceEvolutionResult source = qmomkit::advance_cqmom_moment_source(
            *advanced.pre_source_state, advanced.source_start_time, time_step, example.source, example.phase_space_evolution_options);
        if (!source.ok() || !source.quadrature.has_value())
        {
            return fail(source);
        }
        advanced.moment_source_diagnostics = source.diagnostics;
        advanced.boundary_state = source.state;
        advanced.boundary_quadrature = *source.quadrature;
        advanced.boundary_time = source.final_time;
        advanced.stage_label = "post_source_inversion_before_second_advection";
        qmomkit::CqmomPhaseSpaceEvolutionResult second = qmomkit::advance_cqmom_advection(*advanced.boundary_state, time + half_step, half_step,
                                                                                          example.velocity, example.phase_space_evolution_options);
        if (!second.ok() || !second.quadrature.has_value())
        {
            return fail(second);
        }
        advanced.final_state = second.state;
        advanced.final_quadrature = *second.quadrature;
        advanced.ok = true;
        return advanced;
    }
    if (example.phase_space_evolution_options.scheme == qmomkit::PhaseSpaceEvolutionScheme::advsrc)
    {
        qmomkit::CqmomPhaseSpaceEvolutionResult advection =
            qmomkit::advance_cqmom_advection(state, time, time_step, example.velocity, example.phase_space_evolution_options);
        if (!advection.ok() || !advection.quadrature.has_value())
        {
            return fail(advection);
        }
        advanced.pre_source_state = advection.state;
        advanced.pre_source_quadrature = *advection.quadrature;
        qmomkit::CqmomPhaseSpaceEvolutionResult source =
            qmomkit::advance_cqmom_moment_source(*advanced.pre_source_state, time, time_step, example.source, example.phase_space_evolution_options);
        if (!source.ok() || !source.quadrature.has_value())
        {
            return fail(source);
        }
        advanced.moment_source_diagnostics = source.diagnostics;
        advanced.boundary_state = source.state;
        advanced.boundary_quadrature = *source.quadrature;
        advanced.final_state = source.state;
        advanced.final_quadrature = *source.quadrature;
        advanced.boundary_time = source.final_time;
        advanced.stage_label = "post_source_inversion";
        advanced.ok = true;
        return advanced;
    }

    qmomkit::CqmomPhaseSpaceEvolutionResult source =
        qmomkit::advance_cqmom_moment_source(state, time, time_step, example.source, example.phase_space_evolution_options);
    if (!source.ok() || !source.quadrature.has_value())
    {
        return fail(source);
    }
    advanced.pre_source_state = state;
    advanced.moment_source_diagnostics = source.diagnostics;
    advanced.boundary_state = source.state;
    advanced.boundary_quadrature = *source.quadrature;
    advanced.boundary_time = source.final_time;
    advanced.stage_label = "post_source_inversion_before_advection";
    qmomkit::CqmomPhaseSpaceEvolutionResult advection =
        qmomkit::advance_cqmom_advection(*advanced.boundary_state, time, time_step, example.velocity, example.phase_space_evolution_options);
    if (!advection.ok() || !advection.quadrature.has_value())
    {
        return fail(advection);
    }
    advanced.final_state = advection.state;
    advanced.final_quadrature = *advection.quadrature;
    advanced.ok = true;
    return advanced;
}

DiagnosticReport build_report(const cqmom_5d_non_spherical_particles::PreparedCqmomExample& unbounded,
                              const std::optional<qmomkit::CqmomDiagnostics>& bounded_branch, const BoundaryAdvance& unbounded_boundary,
                              const BoundaryAdvance& bounded_boundary, OffendingNode offending, const std::size_t macro_step)
{
    if (!unbounded_boundary.pre_source_state.has_value() || !unbounded_boundary.boundary_quadrature.has_value())
    {
        throw std::runtime_error("Diagnostic boundary state is incomplete.");
    }
    const qmomkit::ConditioningPath parent_path{{offending.path.node_index_at_level(0)}};
    const std::optional<qmomkit::CqmomDiagnostics> unbounded_branch =
        find_branch_diagnostic(unbounded_boundary.moment_source_diagnostics.post_source_cqmom_diagnostics, parent_path);
    if (!unbounded_branch.has_value())
    {
        throw std::runtime_error("Could not find the offending unbounded branch diagnostic.");
    }

    SourceReplay replay;
    try
    {
        replay = record_phase_space_source_stages(unbounded_boundary, parent_path);
    }
    catch (const std::exception& error)
    {
        throw std::runtime_error(std::string("phase-space source-stage record failed: ") + error.what());
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
        report.bounded_transition = find_branch_transition(bounded_boundary.moment_source_diagnostics.cqmom_node_count_transitions, parent_path);
        report.bounded_accepted_branch = actual_bounded_branch_from_quadrature(*bounded_boundary.boundary_quadrature, parent_path,
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

int cqmom_5d_non_spherical_particles::run_unbounded_diagnostic()
{
    const std::filesystem::path unbounded_configuration_file = "configuration_unbounded.qmom";
    const std::filesystem::path bounded_configuration_file = "configuration.qmom";

    try
    {
        std::optional<cqmom_5d_non_spherical_particles::PreparedCqmomExample> unbounded_prepared =
            cqmom_5d_non_spherical_particles::prepare_cqmom_example(unbounded_configuration_file);
        std::optional<cqmom_5d_non_spherical_particles::PreparedCqmomExample> bounded_prepared =
            cqmom_5d_non_spherical_particles::prepare_cqmom_example(bounded_configuration_file);
        if (!unbounded_prepared.has_value() || !bounded_prepared.has_value())
        {
            return 1;
        }
        cqmom_5d_non_spherical_particles::PreparedCqmomExample unbounded = std::move(*unbounded_prepared);
        cqmom_5d_non_spherical_particles::PreparedCqmomExample bounded = std::move(*bounded_prepared);

        std::cout << "5D CQMOM unbounded chi diagnostic\n";
        std::cout << "configuration: " << unbounded_configuration_file.string() << '\n';
        cqmom_5d_non_spherical_particles::print_phase_space_evolution_settings(unbounded);

        qmomkit::CqmomState unbounded_state = std::move(unbounded.state);
        qmomkit::CqmomState bounded_state = std::move(bounded.state);
        qmomkit::CqmomQuadrature last_unbounded_quadrature = unbounded.initial_quadrature;
        double time = cqmom_5d_non_spherical_particles::cqmom_initial_time;
        std::size_t macro_step = 0;

        while (time < cqmom_5d_non_spherical_particles::cqmom_final_time - 1.0e-12)
        {
            ++macro_step;
            const double current_step =
                std::min(cqmom_5d_non_spherical_particles::cqmom_time_step, cqmom_5d_non_spherical_particles::cqmom_final_time - time);
            BoundaryAdvance unbounded_boundary = advance_to_source_boundary(unbounded, unbounded_state, time, current_step);
            BoundaryAdvance bounded_boundary = advance_to_source_boundary(bounded, bounded_state, time, current_step);
            if (!unbounded_boundary.ok)
            {
                cqmom_5d_non_spherical_particles::write_failure_context(unbounded.output_options.directory, last_unbounded_quadrature, unbounded,
                                                                        time, unbounded_boundary.message);
                std::cerr << "Unbounded diagnostic phase-space evolution failed before an out-of-support chi node was accepted: "
                          << unbounded_boundary.message << '\n';
                return 1;
            }
            if (!bounded_boundary.ok)
            {
                cqmom_5d_non_spherical_particles::write_failure_context(unbounded.output_options.directory, last_unbounded_quadrature, unbounded,
                                                                        time, bounded_boundary.message);
                std::cerr << "Bounded comparison phase-space evolution failed before the diagnostic event: " << bounded_boundary.message << '\n';
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
                        find_branch_diagnostic(bounded_boundary.moment_source_diagnostics.post_source_cqmom_diagnostics, parent_path);
                    report = build_report(unbounded, bounded_branch, unbounded_boundary, bounded_boundary, *offending, macro_step);
                    cqmom_5d_non_spherical_particles::write_report_file(unbounded.output_options.directory, *report, unbounded, bounded);
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
