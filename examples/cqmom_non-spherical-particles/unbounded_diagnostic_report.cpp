#include "unbounded_diagnostic_report.hpp"

#include "qmomkit/qmomkit.hpp"
#include "sampling.hpp"
#include "utilities.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace cqmom_non_spherical_particles
{

namespace
{

constexpr std::size_t chi_dimension = 1;
constexpr double chi_lower = 0.0;
constexpr double chi_upper = 1.0;

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

} // namespace

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
                "with chi outside that interval, writes this report, and stops before the next advection stage.\n\n";

    markdown << "## Run configuration\n\n";
    markdown << "| quantity | bounded workflow | unbounded diagnostic |\n";
    markdown << "|---|---|---|\n";
    markdown << "| configuration | `" << bounded.configuration_file.string() << "` | `" << unbounded.configuration_file.string() << "` |\n";
    markdown << "| sampled initialization | seed " << cqmom_non_spherical_particles::initialization_seed << ", samples "
             << cqmom_non_spherical_particles::initialization_sample_count << " | same seed and sample count |\n";
    markdown << "| node counts | `" << vector_string(unbounded.cqmom_case.hierarchy().node_counts()) << "` | same |\n";
    markdown << "| conditioning order | `" << vector_string(unbounded.cqmom_case.hierarchy().conditioning_order()) << "` | same |\n";
    markdown << "| phase-space evolution scheme | `" << phase_space_evolution_scheme_string(unbounded.phase_space_evolution_options.scheme)
             << "` | same |\n";
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
           "immediately after moment-space fragmentation because the next advection stage could obscure the closure event that produced the "
           "out-of-support value.\n\n";

    markdown << "## Conditional moment sequence\n\n";
    markdown << "The unbounded workflow supplies all six branch-local chi moments to the three-node inversion:\n\n";
    write_moment_table(markdown, report.branch_moments, report.offending.parent_global_weight);
    markdown << '\n';

    markdown << "## Moment-Source Evolution\n\n";
    markdown << "The reported pre-source and accepted post-source states come directly from the same "
                "`advance_cqmom_moment_source(...)` phase-space operation that evaluates the breakage callback. This preserves the active "
                "adaptive CQMOM topology and source-stage representation used by the diagnostic run.\n\n";
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

    markdown << "## Support-localizing matrices\n\n";
    markdown << "For interval `[a,b]`, `Delta_lower(l) = det([m_(i+j+1) - a*m_(i+j)])` and "
                "`Delta_upper(l) = det([b*m_(i+j) - m_(i+j+1)])`. For `[0,1]`, the upper matrix is `[m_(i+j) - m_(i+j+1)]`, "
                "which detects the missing upper chi bound that the classical Stieltjes conditions do not enforce.\n\n";
    markdown << "### Lower-support localizing determinants\n\n";
    markdown << "| determinant | value |\n";
    markdown << "|---|---:|\n";
    for (const MatrixDiagnostic& matrix : report.matrices)
    {
        if (matrix.family == "lower_localizing")
        {
            markdown << "| Delta_lower(" << matrix.l << ") | " << matrix.determinant << " |\n";
        }
    }
    markdown << "\n### Upper-support localizing determinants\n\n";
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
    markdown << "\n### Lower-support localizing matrices\n\n";
    write_interval_matrix_table(markdown, report.matrices, "lower_localizing", "Delta_lower");
    markdown << "\n### Upper-support localizing matrices\n\n";
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
           "bounded interval `[0,1]`, and the upper-support localizing matrix is the first failing condition. The bounded workflow resolves the "
           "event through adaptive branch-order reduction while preserving ordinary user-facing simulation output.\n";
}

} // namespace cqmom_non_spherical_particles
