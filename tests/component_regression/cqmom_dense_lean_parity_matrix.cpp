/*
 * WHAT IS TESTED
 *   That the dense and lean CQMOM moment storage schemes produce the same
 *   reconstruction, over the whole case matrix: phase-space dimension,
 *   conditioning order, node counts, adaptive branch counts, bounded coordinate
 *   support, every Vandermonde solver policy, and all three splitting schemes.
 *
 *   This complements cqmom_dense_lean_parity_3d.cpp, which examines a single 3D
 *   case in close detail. This test is broad rather than deep: it checks
 *   agreement across the matrix but does not examine the intermediate solves.
 *
 * HOW IT IS TESTED
 *   Each case is built once as a dense raw moment tensor and once as the lean
 *   required-moment set extracted from that same tensor, so both paths start
 *   from the identical numbers. Both are reconstructed and the accepted
 *   quadratures compared node by node. Evolution cases are advanced for the
 *   same steps under the same splitting scheme and their final quadratures
 *   compared.
 *
 * PASS CRITERIA
 *   Status, node count and the branches adaptive reconstruction kept agree
 *   exactly. Weights and abscissas agree to 5e-8 absolute or 5e-10 relative,
 *   the tolerance the detailed 3D parity test already uses for the same
 *   comparison: both paths run the same core but in a different order of
 *   operations, so they agree to round-off rather than exactly.
 *
 * WHEN IT IS TESTED
 *   Debug and Release. The checks are explicit and do not rely on assert().
 *
 * LIFETIME
 *   This test exists to record, while both storage schemes are still present,
 *   that removing the dense scheme loses no numerical behaviour. It is expected
 *   to be deleted along with the dense storage scheme, because what it compares
 *   will no longer exist.
 */

#include "../verification/verification_helpers.hpp"
#include "cqmom_case_matrix.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/phase_space_evolution.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{

constexpr double parity_absolute_tolerance = 5.0e-8;
constexpr double parity_relative_tolerance = 5.0e-10;

[[nodiscard]] qmomkit::CqmomDenseCase dense_case_of(const qmomkit::case_matrix::StaticCase& case_)
{
    return qmomkit::CqmomDenseCase{qmomkit::case_matrix::hierarchy_of(case_),
                                   qmomkit::case_matrix::dense_moments_of(case_),
                                   {},
                                   {},
                                   case_.adaptive_branch_counts,
                                   qmomkit::CqmomDiagnosticLevel::basic,
                                   case_.vandermonde_solver,
                                   case_.coordinate_support};
}

// Reports which case and which node disagreed. The executable can be run over
// the whole matrix at once, where the CTest case name is not available to
// identify the failure.
void check_agrees(const std::string& case_name, const char* quantity, const std::size_t node, const double lean_value, const double dense_value)
{
    const double difference = std::abs(lean_value - dense_value);
    const double relative = difference / std::max(std::abs(dense_value), 1.0);
    if (difference <= parity_absolute_tolerance || relative <= parity_relative_tolerance)
    {
        return;
    }
    std::ostringstream detail;
    detail << std::setprecision(17) << "case " << case_name << ", " << quantity << " at node " << node << ": lean " << lean_value << ", dense "
           << dense_value << ", absolute difference " << difference << ", relative difference " << relative;
    qmomkit::verification::fail_check("dense/lean parity comparison", __FILE__, __LINE__, detail.str());
}

void compare_quadratures(const std::string& case_name, const qmomkit::CqmomQuadrature& dense, const qmomkit::CqmomQuadrature& lean)
{
    QMOMKIT_VERIFY_CHECK(dense.number_of_dimensions() == lean.number_of_dimensions());
    QMOMKIT_VERIFY_CHECK(dense.weights().size() == lean.weights().size());
    QMOMKIT_VERIFY_CHECK(dense.total_nodes() == lean.total_nodes());
    QMOMKIT_VERIFY_CHECK(dense.has_rectangular_storage() == lean.has_rectangular_storage());
    QMOMKIT_VERIFY_CHECK(dense.adaptive_branch_counts_enabled() == lean.adaptive_branch_counts_enabled());

    for (std::size_t node = 0; node < dense.weights().size(); ++node)
    {
        QMOMKIT_VERIFY_CHECK(dense.terminal_paths()[node].to_conditioning_path() == lean.terminal_paths()[node].to_conditioning_path());
        check_agrees(case_name, "weight", node, lean.weights()[node], dense.weights()[node]);
        for (std::size_t dimension = 0; dimension < dense.number_of_dimensions(); ++dimension)
        {
            check_agrees(case_name, "abscissa", node, lean.abscissas()[dimension][node], dense.abscissas()[dimension][node]);
        }
    }
}

void run_static_case(const qmomkit::case_matrix::StaticCase& case_)
{
    const qmomkit::CqmomResult dense = qmomkit::compute_cqmom(dense_case_of(case_));
    const qmomkit::CqmomResult lean = qmomkit::compute_cqmom(qmomkit::case_matrix::lean_case_of(case_));

    QMOMKIT_VERIFY_CHECK(dense.status == lean.status);
    QMOMKIT_VERIFY_CHECK(dense.ok() == lean.ok());
    QMOMKIT_VERIFY_CHECK(dense.ok());
    QMOMKIT_VERIFY_CHECK(dense.diagnostics.size() == lean.diagnostics.size());
    compare_quadratures(case_.name, dense.quadrature, lean.quadrature);
}

void run_evolution_case(const qmomkit::case_matrix::EvolutionCase& case_)
{
    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::case_matrix::matrix_velocity();
    const qmomkit::CqmomMomentSource source = qmomkit::case_matrix::matrix_source();

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = case_.scheme;

    qmomkit::CqmomState dense_state = qmomkit::CqmomCase(dense_case_of(case_.reconstruction)).initial_state();
    qmomkit::CqmomState lean_state = qmomkit::CqmomCase(qmomkit::case_matrix::lean_case_of(case_.reconstruction)).initial_state();

    QMOMKIT_VERIFY_CHECK(dense_state.moment_storage() == qmomkit::CqmomMomentStorage::dense);
    QMOMKIT_VERIFY_CHECK(lean_state.moment_storage() == qmomkit::CqmomMomentStorage::lean);

    double time = 0.0;
    for (std::size_t step = 0; step < case_.step_count; ++step)
    {
        const qmomkit::CqmomPhaseSpaceEvolutionResult dense =
            qmomkit::advance_cqmom_phase_space(dense_state, time, case_.time_step, velocity, source, options);
        const qmomkit::CqmomPhaseSpaceEvolutionResult lean =
            qmomkit::advance_cqmom_phase_space(lean_state, time, case_.time_step, velocity, source, options);

        QMOMKIT_VERIFY_CHECK(dense.status == lean.status);
        QMOMKIT_VERIFY_CHECK(dense.ok() == lean.ok());
        QMOMKIT_VERIFY_CHECK(dense.ok());
        QMOMKIT_VERIFY_CHECK(dense.diagnostics.failed_stage == lean.diagnostics.failed_stage);
        QMOMKIT_VERIFY_CHECK(dense.final_time == lean.final_time);

        dense_state = dense.state;
        lean_state = lean.state;
        time = dense.final_time;
    }

    compare_quadratures(case_.reconstruction.name, dense_state.quadrature(), lean_state.quadrature());
}

} // namespace

int main(const int argc, char** argv)
{
    const std::string selected_case = argc == 2 ? argv[1] : "";
    const auto selected = [&selected_case](const std::string& name) { return selected_case.empty() || selected_case == name; };
    bool matched = selected_case.empty();

    for (const qmomkit::case_matrix::StaticCase& case_ : qmomkit::case_matrix::static_cases())
    {
        if (selected(case_.name))
        {
            matched = true;
            run_static_case(case_);
        }
    }

    for (const qmomkit::case_matrix::EvolutionCase& case_ : qmomkit::case_matrix::evolution_cases())
    {
        if (selected(case_.reconstruction.name))
        {
            matched = true;
            run_evolution_case(case_);
        }
    }

    if (!matched)
    {
        std::cerr << "Unknown dense/lean parity case '" << selected_case << "'.\n";
        return 2;
    }

    return 0;
}
