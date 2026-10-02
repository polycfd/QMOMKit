/*
 * WHAT IS TESTED
 *   That CQMOM reconstruction and phase-space evolution produce exactly the
 *   results recorded here, across a matrix spanning phase-space dimension,
 *   conditioning order, node counts, adaptive branch counts, bounded coordinate
 *   support, every Vandermonde solver policy, and all three splitting schemes.
 *
 * HOW IT IS TESTED
 *   Three checks per static case.
 *
 *   1. Analytical recovery. Every case is built from an explicit list of
 *      quadrature nodes, and its raw moments are computed from those nodes.
 *      That node list is the only quadrature reproducing those moments for the
 *      configured node counts, so reconstruction must return the nodes the case
 *      was built from. This is an analytical check: it does not depend on any
 *      recorded value, and would catch a wrong result even if the recorded
 *      values were themselves wrong.
 *
 *   2. Prepared-path agreement. The same case is reconstructed a second time
 *      through CqmomPreparedCase, which reuses a reconstruction plan built once from
 *      the hierarchy instead of deriving it per call. Both routes must return
 *      the identical quadrature, down to the last bit.
 *
 *   3. Recorded values. The accepted weights and abscissas are compared against
 *      the values recorded below, which pins the result far more tightly than
 *      recovery alone: it also fixes node ordering, conditioning-path ordering,
 *      and which branches adaptive reconstruction kept.
 *
 *   In addition, an initial state is created both from the case and from
 *   CqmomPreparedCase. The two must be identical, and the state's quadrature
 *   must be the one reconstruction returned.
 *
 *   Evolution cases carry the recorded check and the initial-state check: the
 *   two initial states are advanced side by side and must stay identical after
 *   every step. Recovery does not apply, because once the state has been
 *   advanced it is no longer the node list the case was built from.
 *
 * PASS CRITERIA
 *   Recovery: every node the case was built from is matched by an accepted
 *   node, with weight and coordinates agreeing to 1e-12 absolute.
 *   Prepared-path and initial-state agreement: bit-exact, because both routes
 *   run in the same build and any difference is a real one.
 *   Recorded values: status and node count exactly; every weight and abscissa
 *   within 1e-12 relative, with a 1e-15 absolute floor.
 *
 *   The recorded comparison cannot be exact. The values were produced by one
 *   compiler, and another reproduces the same computation with different last
 *   bits, so an exact comparison reports a toolchain difference as a
 *   regression. The budget is set from measurement, far above that noise and
 *   far below any change of algorithm or ordering; the rationale and the
 *   measured figures are at check_recorded_value below. A failure reports the
 *   absolute difference, the number of representable values separating the two
 *   numbers, and the budget it exceeded.
 *
 *   This bounds how tightly a refactor can be pinned by this file alone. To
 *   show that a refactor changed nothing at all, run this suite before and
 *   after the change with one compiler and compare: within a single toolchain
 *   the results are reproducible bit for bit.
 *
 * WHEN IT IS TESTED
 *   Debug and Release. The checks are explicit and do not rely on assert().
 */

#include "../verification/verification_helpers.hpp"
#include "cqmom_case_matrix.hpp"
#include "qmomkit/cqmom.hpp"
#include "qmomkit/phase_space_evolution.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace
{

constexpr double node_recovery_tolerance = 1.0e-12;

struct ReferenceQuadrature
{
    const char* name;
    int status;
    std::vector<double> weights;
    // Abscissas in dimension-major order: all nodes of dimension 0, then all
    // nodes of dimension 1, and so on.
    std::vector<double> abscissas;
};

#include "cqmom_reference_values.inc"

[[nodiscard]] const ReferenceQuadrature& reference_for(const std::string& name)
{
    for (const ReferenceQuadrature& reference : reference_quadratures())
    {
        if (name == reference.name)
        {
            return reference;
        }
    }
    std::cerr << "No recorded reference for case '" << name << "'.\n";
    std::exit(EXIT_FAILURE);
}

// How many representable floating-point values lie between the two numbers, so
// that a failure says whether the result moved by a rounding step or by
// something that matters.
[[nodiscard]] long double representable_values_between(const double actual, const double expected)
{
    if (actual == expected)
    {
        return 0.0L;
    }
    if (!std::isfinite(actual) || !std::isfinite(expected))
    {
        return std::numeric_limits<long double>::infinity();
    }
    std::int64_t actual_bits = 0;
    std::int64_t expected_bits = 0;
    std::memcpy(&actual_bits, &actual, sizeof(actual_bits));
    std::memcpy(&expected_bits, &expected, sizeof(expected_bits));
    if ((actual_bits < 0) != (expected_bits < 0))
    {
        return std::numeric_limits<long double>::infinity();
    }
    return static_cast<long double>(std::llabs(actual_bits - expected_bits));
}

// The recorded values were produced by one compiler. Another compiler evaluates
// the same algorithm with different instruction selection and fused multiply-add
// contraction, so the last bits differ even when nothing about the computation
// has changed. Comparing exactly would therefore report a toolchain difference
// as a regression, so the recorded check is relative.
//
// The budget is chosen from measurement: across the whole matrix, GCC differs
// from the recording toolchain by at most 27 representable values, which is
// about 6e-15 relative for the magnitudes recorded here. A change of algorithm,
// node ordering or branch selection moves results by many orders more than that,
// so this budget separates the two.
constexpr double recorded_relative_tolerance = 1.0e-12;
// Guards a future recording close to zero, where the relative term alone would
// be unusable. Every value recorded today is far enough from zero that the
// relative term governs.
constexpr double recorded_absolute_floor = 1.0e-15;

void check_recorded_value(const std::string& case_name, const char* quantity, const std::size_t position, const double actual, const double expected)
{
    const double difference = std::abs(actual - expected);
    const double allowed = std::max(recorded_absolute_floor, recorded_relative_tolerance * std::abs(expected));
    if (std::isfinite(actual) && std::isfinite(expected) && difference <= allowed)
    {
        return;
    }
    std::ostringstream detail;
    detail << std::setprecision(17) << "case " << case_name << ", " << quantity << '[' << position << "]: produced " << actual << ", recorded "
           << expected << ", differing by " << difference << " (" << static_cast<double>(representable_values_between(actual, expected))
           << " representable values), which exceeds the permitted " << allowed;
    qmomkit::verification::fail_check("recorded value comparison", __FILE__, __LINE__, detail.str());
}

void compare_against_reference(const std::string& case_name, const qmomkit::CqmomQuadrature& quadrature, const int status)
{
    const ReferenceQuadrature& reference = reference_for(case_name);
    QMOMKIT_VERIFY_CHECK(status == reference.status);
    QMOMKIT_VERIFY_CHECK(quadrature.weights().size() == reference.weights.size());

    for (std::size_t node = 0; node < reference.weights.size(); ++node)
    {
        check_recorded_value(case_name, "weight", node, quadrature.weights()[node], reference.weights[node]);
    }

    const std::size_t dimensions = quadrature.number_of_dimensions();
    QMOMKIT_VERIFY_CHECK(dimensions * reference.weights.size() == reference.abscissas.size());
    std::size_t position = 0;
    for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
    {
        for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
        {
            check_recorded_value(case_name, "abscissa", position, quadrature.abscissas()[dimension][node], reference.abscissas[position]);
            ++position;
        }
    }
}

/*
 * The nodes a case was built from are the only quadrature reproducing its
 * moments, so each of them must appear among the accepted nodes. They are
 * matched by nearest position rather than by index, because node ordering is
 * already fixed by the recorded comparison above.
 *
 * Each accepted node may be claimed only once. Without that, a reconstruction
 * that returned the same node twice and dropped another could still match every
 * node the case was built from, provided the dropped one happened to lie close
 * to the duplicate.
 */
void check_recovers_generating_nodes(const qmomkit::case_matrix::StaticCase& case_, const qmomkit::CqmomQuadrature& quadrature)
{
    const std::size_t dimensions = quadrature.number_of_dimensions();
    QMOMKIT_VERIFY_CHECK(quadrature.weights().size() == case_.nodes.size());
    std::vector<bool> claimed(quadrature.weights().size(), false);

    for (const qmomkit::case_matrix::QuadratureNode& expected : case_.nodes)
    {
        std::size_t closest = 0;
        double closest_distance = std::numeric_limits<double>::infinity();
        for (std::size_t node = 0; node < quadrature.weights().size(); ++node)
        {
            if (claimed[node])
            {
                continue;
            }
            double distance = 0.0;
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
            {
                const double offset = quadrature.abscissas()[dimension][node] - expected.coordinate[dimension];
                distance += offset * offset;
            }
            if (distance < closest_distance)
            {
                closest_distance = distance;
                closest = node;
            }
        }

        QMOMKIT_VERIFY_CHECK(closest_distance < std::numeric_limits<double>::infinity());
        claimed[closest] = true;
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
        {
            QMOMKIT_VERIFY_NEAR(quadrature.abscissas()[dimension][closest], expected.coordinate[dimension], node_recovery_tolerance);
        }
        QMOMKIT_VERIFY_NEAR(quadrature.weights()[closest], expected.weight, node_recovery_tolerance);
    }
}

void check_quadratures_identical(const qmomkit::CqmomQuadrature& actual, const qmomkit::CqmomQuadrature& expected)
{
    QMOMKIT_VERIFY_CHECK(actual.weights().size() == expected.weights().size());
    QMOMKIT_VERIFY_CHECK(actual.number_of_dimensions() == expected.number_of_dimensions());

    for (std::size_t node = 0; node < actual.weights().size(); ++node)
    {
        QMOMKIT_VERIFY_CHECK(actual.weights()[node] == expected.weights()[node]);
        for (std::size_t dimension = 0; dimension < actual.number_of_dimensions(); ++dimension)
        {
            QMOMKIT_VERIFY_CHECK(actual.abscissas()[dimension][node] == expected.abscissas()[dimension][node]);
        }
    }
}

/*
 * Reconstructing through CqmomPreparedCase reuses a reconstruction plan built once from
 * the hierarchy, where the explicit case derives one per call. The two routes
 * must agree exactly: the plan is bookkeeping about which moments are needed
 * and in what order, and does not enter the arithmetic.
 */
void check_prepared_path_agrees(const qmomkit::CqmomResult& direct, const qmomkit::CqmomResult& prepared)
{
    QMOMKIT_VERIFY_CHECK(direct.status == prepared.status);
    check_quadratures_identical(direct.quadrature, prepared.quadrature);
}

/*
 * A state created from the case builds its reconstruction plan in that call,
 * where a state created through CqmomPreparedCase takes the plan the prepared
 * case already holds. The plan does not enter the arithmetic, so the two states
 * must agree exactly: moments, active branch counts, and accepted quadrature.
 */
void check_states_identical(const qmomkit::CqmomState& actual, const qmomkit::CqmomState& expected)
{
    QMOMKIT_VERIFY_CHECK(actual.moment_values() == expected.moment_values());

    const std::vector<qmomkit::CqmomActiveBranchCount>& actual_counts = actual.active_branch_counts();
    const std::vector<qmomkit::CqmomActiveBranchCount>& expected_counts = expected.active_branch_counts();
    QMOMKIT_VERIFY_CHECK(actual_counts.size() == expected_counts.size());
    for (std::size_t branch = 0; branch < actual_counts.size(); ++branch)
    {
        QMOMKIT_VERIFY_CHECK(actual_counts[branch].conditioning_level == expected_counts[branch].conditioning_level);
        QMOMKIT_VERIFY_CHECK(actual_counts[branch].parent_path == expected_counts[branch].parent_path);
        QMOMKIT_VERIFY_CHECK(actual_counts[branch].phase_space_dimension == expected_counts[branch].phase_space_dimension);
        QMOMKIT_VERIFY_CHECK(actual_counts[branch].configured_node_count == expected_counts[branch].configured_node_count);
        QMOMKIT_VERIFY_CHECK(actual_counts[branch].active_node_count == expected_counts[branch].active_node_count);
    }

    check_quadratures_identical(actual.quadrature(), expected.quadrature());
}

void run_static_case(const qmomkit::case_matrix::StaticCase& case_)
{
    const qmomkit::CqmomCase explicit_case = qmomkit::case_matrix::lean_case_of(case_);
    const qmomkit::CqmomResult result = qmomkit::compute_cqmom(explicit_case);
    QMOMKIT_VERIFY_CHECK(result.ok());

    check_recovers_generating_nodes(case_, result.quadrature);
    const qmomkit::CqmomPreparedCase cqmom_prepared_case(explicit_case);
    check_prepared_path_agrees(result, qmomkit::compute_cqmom(cqmom_prepared_case));

    const qmomkit::CqmomState initial = explicit_case.initial_state();
    check_quadratures_identical(initial.quadrature(), result.quadrature);
    check_states_identical(initial, cqmom_prepared_case.initial_state());

    compare_against_reference(case_.name, result.quadrature, static_cast<int>(result.status));
}

void run_evolution_case(const qmomkit::case_matrix::EvolutionCase& case_)
{
    const qmomkit::PhaseSpaceVelocity velocity = qmomkit::case_matrix::matrix_velocity();
    const qmomkit::CqmomMomentSource source = qmomkit::case_matrix::matrix_source();

    qmomkit::CqmomPhaseSpaceEvolutionOptions options;
    options.scheme = case_.scheme;

    const qmomkit::CqmomCase explicit_case = qmomkit::case_matrix::lean_case_of(case_.reconstruction);
    qmomkit::CqmomState state = qmomkit::CqmomPreparedCase(explicit_case).initial_state();
    qmomkit::CqmomState case_state = explicit_case.initial_state();
    check_states_identical(case_state, state);

    double time = 0.0;
    int status = 0;
    for (std::size_t step = 0; step < case_.step_count; ++step)
    {
        const qmomkit::CqmomPhaseSpaceEvolutionResult evolved =
            qmomkit::advance_cqmom_phase_space(state, time, case_.time_step, velocity, source, options);
        status = static_cast<int>(evolved.status);
        QMOMKIT_VERIFY_CHECK(evolved.ok());

        const qmomkit::CqmomPhaseSpaceEvolutionResult case_evolved =
            qmomkit::advance_cqmom_phase_space(case_state, time, case_.time_step, velocity, source, options);
        QMOMKIT_VERIFY_CHECK(case_evolved.status == evolved.status);
        QMOMKIT_VERIFY_CHECK(case_evolved.final_time == evolved.final_time);
        check_states_identical(case_evolved.state, evolved.state);

        state = evolved.state;
        case_state = case_evolved.state;
        time = evolved.final_time;
    }

    compare_against_reference(case_.reconstruction.name, state.quadrature(), status);
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
        std::cerr << "Unknown CQMOM reference-value case '" << selected_case << "'.\n";
        return 2;
    }

    return 0;
}
